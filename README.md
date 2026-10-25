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
With the current parameter set the model diverges within the first 2 days for every genotype tested (including no variants and `--no-ebv`) (dysbiosis `D` and portal load `f` grow without bound) and the solver stops. The pipeline now reports this and exits with an error rather than printing partial values. The Results figures in `AASLD_Abstract.md` are not yet reproduced by this code.
