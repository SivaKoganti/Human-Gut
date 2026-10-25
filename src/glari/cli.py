"""Command-line entry point: `glari` or `python -m glari`."""

import argparse
import sys

import numpy as np

from .model import (
    CLINVAR_DB,
    DEFAULT_VARIANTS,
    ModelParams,
    initial_state,
    sample_patient_genomics_mc,
    simulate,
)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="glari",
        description="Gut-Liver Axis Risk & Intervention simulation pipeline.",
    )
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
