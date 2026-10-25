# Multi-Scale Gut–Liver–Immune Axis Computational Model

This repository contains the multi-scale mathematical modeling pipeline, AASLD abstract, and computational framework for investigating the gut–liver–immune axis, EBV-associated gastric carcinoma (EBVaGC), ClinVar genomic risk profiles (`CDH1`, `NOD2`), and LSEC capillarization.

## Co-Authors
- **Siva Koganti** (Department of Medicine, University of Florida) - Co-Corresponding Author
- **Karuna Rasineni** - Co-Corresponding Author

## Model Overview
- **Percolation Mechanics**: Site-bond percolation ($p_c \approx 0.50$) modeling mucosal barrier breach.
- **LSEC Sieve Filtration**: Hong Wang Fourier restriction estimates governing fenestral diameter ($\delta \to 25\text{ nm}$).
- **Kinetic Dispersive Waves**: Yu Deng fluid dynamics modeling portal flow of LPS, $EBER1/2$ non-coding RNAs, and viral entities.
- **Network Resilience**: June Huh matroid Hodge theory log-concavity ($\lambda_H$).

## Contents
- `glari_pipeline.py`: Executable Python engine with stiff ODE solver (`Radau`) and Monte Carlo variant sampling.
- `AASLD_Abstract.md`: Official AASLD Liver Meeting abstract submission template.
- `requirements.txt`: Python package dependencies.

## Usage
```bash
python -m venv .venv
.venv/Scripts/activate        # Windows; use `source .venv/bin/activate` on macOS/Linux
pip install -r requirements.txt
python glari_pipeline.py --seed 42
```

Options: `--seed` (reproducible genomic sampling), `--days` (horizon, default 30), `--variants` (ClinVar IDs), `--no-ebv`.

## Status
The original equations diverged within 2 days for every genotype (unbounded dysbiosis feedback, saturating clearance, unbounded integrity loss). Three structural fixes keep the system bounded:
1. Logistic cap on dysbiosis growth: `0.20 * f/(1+f) * D * (1 - D/D_max)`.
2. Integrity loss scales with remaining junctions (`* p`); `lam` is no longer clipped in the rate equations.
3. Kupffer-cell clearance term `k_kupffer * f` on portal load.

`D_max` and `k_kupffer` in `ModelParams` are **placeholders** chosen only to keep the system bounded. They are not fitted to data and must be replaced with literature-derived values.

Known gap: a healthy, untreated patient (`--variants --no-ebv`) still declines (Day-30 `p = 0.02`), because the initial state is not a steady state. Gut influx scale and `D_0` still need calibration. The Results figures in `AASLD_Abstract.md` are not reproduced by this code.
