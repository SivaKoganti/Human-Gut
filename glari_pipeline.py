"""GLARI: Gut-Liver Axis Risk & Intervention simulation pipeline.

Integrates a five-state ODE model of the gut-liver-immune axis for a patient
genotype, with and without a multi-target intervention, and reports Day-30
outcomes.

State vector y = [p, f, delta, lam, D]
    p      mucosal tight-junction integrity (percolation occupancy, 0-1)
    f      portal endotoxin/viral load (EU/mL)
    delta  LSEC fenestral diameter (nm)
    lam    microbial network resilience (0-1)
    D      pathobiont dysbiosis index
"""

import argparse
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

# ==============================================================================
# 1. ClinVar Genomic Risk Profiler (Monte Carlo Penetrance)
# ==============================================================================
CLINVAR_DB = pd.DataFrame([
    {"gene": "CDH1",   "variant_id": "VCV000012423", "p_impact_mean": -0.15, "p_impact_std": 0.03,  "gc_weight": 2.5},
    {"gene": "TFF1",   "variant_id": "VCV000018092", "p_impact_mean": -0.10, "p_impact_std": 0.02,  "gc_weight": 1.8},
    {"gene": "NOD2",   "variant_id": "VCV000004077", "p_impact_mean": -0.08, "p_impact_std": 0.01,  "gc_weight": 1.2},
    {"gene": "TLR4",   "variant_id": "VCV000014302", "p_impact_mean": -0.05, "p_impact_std": 0.01,  "gc_weight": 1.4},
    {"gene": "CLEC4M", "variant_id": "VCV000441201", "p_impact_mean": -0.02, "p_impact_std": 0.005, "gc_weight": 1.1},
])

DEFAULT_VARIANTS = ("VCV000012423", "VCV000004077")  # CDH1 + NOD2


def sample_patient_genomics_mc(variant_ids, db, rng, num_samples=1):
    """Sample integrity penalties and multiplicative GC risk weights for a genotype."""
    unknown = set(variant_ids) - set(db["variant_id"])
    if unknown:
        raise ValueError(f"Unknown ClinVar variant IDs: {sorted(unknown)}")

    p_penalties = np.zeros(num_samples)
    gc_weights = np.ones(num_samples)
    for _, row in db[db["variant_id"].isin(variant_ids)].iterrows():
        p_penalties += rng.normal(row["p_impact_mean"], row["p_impact_std"], num_samples)
        gc_weights *= row["gc_weight"]
    return p_penalties, gc_weights


# ==============================================================================
# 2. Model Parameters
# ==============================================================================
@dataclass(frozen=True)
class ModelParams:
    p_c: float = 0.50          # Percolation critical threshold
    delta_max: float = 105.0   # Baseline LSEC fenestra (nm)
    delta_min: float = 25.0    # Capillarization floor (nm)
    sigma: float = 0.08        # Transition width
    gc_stage: float = 2.0
    ebv_active: bool = True
    gc_risk_factor: float = 1.0
    therapy_start: float = 15.0  # Day the combined intervention begins
    # PLACEHOLDERS: chosen only to keep the system bounded, not fitted to data.
    # Replace with literature-derived values before reporting results.
    D_max: float = 10.0        # Dysbiosis carrying capacity (logistic cap on D)
    k_kupffer: float = 5.0     # Kupffer-cell first-order endotoxin clearance (1/day)

    @property
    def ebv_multiplier(self):
        return 2.2 if self.ebv_active else 1.0


def initial_state(p_genetic_penalty):
    p_0 = float(np.clip(0.65 + p_genetic_penalty, 0.10, 0.95))
    return [p_0, 0.40, 85.0, 0.68, 2.20]


# ==============================================================================
# 3. Stiff ODE System (Radau Integration)
# ==============================================================================
def gut_liver_rhs(t, y, params, apply_therapy):
    p = np.clip(y[0], 0.01, 1.0)
    f = max(0.0, y[1])
    delta = np.clip(y[2], params.delta_min, params.delta_max)
    # lam is not clipped: its loss term is proportional to lam, so it stays >= 0
    # on its own, and clipping would decouple the rate from the stored state.
    lam = max(0.0, y[3])
    D = max(0.1, y[4])

    phi_gut = 5.0 * (1.0 + 0.35 * params.gc_stage * params.gc_risk_factor) * params.ebv_multiplier

    P_percolation = (phi_gut * (1.0 + 0.5 * D)) / (1.0 + np.exp((p - params.p_c) / params.sigma))
    R_sieve = 2.5 * f * ((delta / params.delta_max) ** 2) / (1.0 + 0.3 * f)
    R_network = lam * f
    R_kupffer = params.k_kupffer * f

    dp_dt = -0.15 * D * (f / (1.0 + f)) * p + 0.02 * (1.0 - p)  # loss scales with remaining junctions
    df_dt = P_percolation - R_sieve - R_network - R_kupffer
    ddelta_dt = -0.70 * f * (delta - params.delta_min) + 0.1 * (params.delta_max - delta)
    dlam_dt = -0.25 * (D / (1.0 + D)) * f * lam + 0.04 * (0.92 - lam)
    dD_dt = 0.20 * (f / (1.0 + f)) * D * (1.0 - D / params.D_max) - 0.05 * (D - 1.0)

    if apply_therapy:
        dp_dt += 0.40 * (0.85 - p)          # Oral Butyrate / SCFA tight-junction rescue
        ddelta_dt += 3.00 * (95.0 - delta)  # sGC stimulator fenestration restoration
        dlam_dt += 0.30 * (0.88 - lam)      # Synbiotic probiotic network stabilization
        dD_dt -= 0.35 * (D - 1.0)           # Pathobiont suppression

    return [dp_dt, df_dt, ddelta_dt, dlam_dt, dD_dt]


def simulate(y0, params, days, treat):
    """Integrate the model; therapy (if any) switches on at params.therapy_start.

    The intervention is a discontinuity in the RHS, so the integration is split
    at the switch time rather than letting the stiff solver step across it.
    """
    switch = params.therapy_start if treat and params.therapy_start < days else days
    segments = [(0.0, switch, False)]
    if switch < days:
        segments.append((switch, days, True))

    state = y0
    for t0, t1, therapy_on in segments:
        sol = solve_ivp(gut_liver_rhs, (t0, t1), state, method="Radau",
                        args=(params, therapy_on), rtol=1e-6, atol=1e-8)
        if not sol.success:
            raise RuntimeError(f"Solver failed at day {sol.t[-1]:.2f}: {sol.message}")
        state = sol.y[:, -1]
    return state


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for reproducible genomic sampling")
    parser.add_argument("--days", type=float, default=30.0, help="Simulation horizon in days")
    parser.add_argument("--variants", nargs="*", default=list(DEFAULT_VARIANTS),
                        help="ClinVar variant IDs carried by the patient")
    parser.add_argument("--no-ebv", action="store_true", help="Simulate EBV-negative disease")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.days <= 0:
        sys.exit("--days must be positive")

    rng = np.random.default_rng(args.seed)
    try:
        p_penalties, gc_weights = sample_patient_genomics_mc(args.variants, CLINVAR_DB, rng)
    except ValueError as exc:
        sys.exit(str(exc))

    params = ModelParams(gc_risk_factor=float(gc_weights[0]), ebv_active=not args.no_ebv)
    y0 = initial_state(p_penalties[0])

    try:
        untreated = simulate(y0, params, args.days, treat=False)
        treated = simulate(y0, params, args.days, treat=True)
    except RuntimeError as exc:
        sys.exit(f"Simulation did not complete: {exc}")

    print("GLARI Simulation Pipeline Complete.")
    print(f"Untreated Final Load (Day {args.days:g}): {untreated[1]:.2f} EU/mL | Integrity p = {untreated[0]:.2f}")
    print(f"Treated Final Load (Day {args.days:g}):   {treated[1]:.2f} EU/mL | Integrity p = {treated[0]:.2f}")


if __name__ == "__main__":
    main()
