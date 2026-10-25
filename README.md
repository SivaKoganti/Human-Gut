# Multi-Scale Gut–Liver–Immune Axis Computational Model

[![CI](https://github.com/SivaKoganti/Human-Gut/actions/workflows/ci.yml/badge.svg)](https://github.com/SivaKoganti/Human-Gut/actions/workflows/ci.yml) [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Project site: https://sivakoganti.github.io/Human-Gut/

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
- `src/glari/model.py`: genotype sampling, parameters, and the five-state stiff ODE system (`Radau`).
- `src/glari/cli.py`: `glari` command-line entry point.
- `tests/`: pytest suite (reproducibility, physical bounds, therapy effect, CLI validation).
- `docs/index.html`: project landing page, deployed to GitHub Pages.
- `AASLD_Abstract.md`: AASLD Liver Meeting abstract draft.
- `CITATION.cff`, `LICENSE` (MIT), `SECURITY.md`.

## Run locally
Requires Python 3.10+.
```bash
git clone https://github.com/SivaKoganti/Human-Gut.git
cd Human-Gut
python -m venv .venv
.venv\Scripts\activate        # macOS/Linux: source .venv/bin/activate
pip install -e ".[dev]"
glari --seed 42
pytest
```

Options: `--seed` (reproducible genomic sampling), `--days` (horizon, default 30), `--variants` (ClinVar IDs; pass none for a variant-free patient), `--no-ebv`. `python -m glari` works too.

Use it from Python:
```python
import numpy as np
from glari import CLINVAR_DB, ModelParams, initial_state, sample_patient_genomics_mc, simulate

penalties, weights = sample_patient_genomics_mc(["VCV000012423"], CLINVAR_DB, np.random.default_rng(42))
params = ModelParams(gc_risk_factor=float(weights[0]))
p, f, delta, lam, D = simulate(initial_state(penalties[0]), params, days=30, treat=True)
```

## Citing and license
Please cite using `CITATION.cff` (GitHub's "Cite this repository" button). Released under the MIT License. Not for clinical use.

## Status
The original equations diverged within 2 days for every genotype (unbounded dysbiosis feedback, saturating clearance, unbounded integrity loss). Three structural fixes keep the system bounded:
1. Logistic cap on dysbiosis growth: `0.20 * f/(1+f) * D * (1 - D/D_max)`.
2. Integrity loss scales with remaining junctions (`* p`); `lam` is no longer clipped in the rate equations.
3. Kupffer-cell clearance term `k_kupffer * f` on portal load.

`D_max` and `k_kupffer` in `ModelParams` are **placeholders** chosen only to keep the system bounded. They are not fitted to data and must be replaced with literature-derived values.

Known gap: a healthy, untreated patient (`--variants --no-ebv`) still declines (Day-30 `p = 0.02`), because the initial state is not a steady state. Gut influx scale and `D_0` still need calibration. The Results figures in `AASLD_Abstract.md` are not reproduced by this code.
