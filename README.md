# Individual-level reconstruction of anti-Spike antibody kinetics under hybrid immunity

Code accompanying the manuscript *"Individual-level reconstruction of longitudinal
anti-Spike binding antibody kinetics under hybrid immunity."*

The framework fits a piecewise-exponential nonlinear mixed-effects (NLME) model to
sparse, irregular longitudinal anti-Spike antibody measurements, and reconstructs
continuous individual-level antibody trajectories across successive immune events
(second dose, booster dose, and breakthrough infection).

## Scripts

Run `01` first; it produces the fitted model object used by all other scripts.
Each script is run once per vaccination-group dataset.

| Script | Description |
|--------|-------------|
| `01_fit_nlme_model.R` | Fits the piecewise-exponential NLME model by SAEM (the core framework). Outputs the fitted model object and the long-format observations. |
| `02_simulation_recovery.R` | Simulation-recovery analysis of individual random effects (identifiability). |
| `03_bi_date_perturbation.R` | Sensitivity of post-infection estimates to uncertainty in the breakthrough-infection date. |
| `04_taub_fixed_refit.R` | Refits a reduced model with the post-infection time-to-peak fixed, and compares fit criteria (OFV/AIC/BIC) with the full model. |
| `05_reduced_sampling_density.R` | Reconstruction error as a function of reduced within-person sampling density. |
| `06_heldout_prediction.R` | Internal validation by leave-one-time-point-out held-out prediction. |

Scripts `02`, `03`, `05`, and `06` re-estimate only the individual random effects by
maximum a posteriori (MAP) estimation with the population model held fixed; they do not
refit the population model (no SAEM). Script `04` performs one SAEM refit of the reduced
model.

## Input files

Input filenames are written as `sample1` / `sample2` placeholders at the top of each
script; replace them with your own files.

- **`sample1`** — empirical antibody data, wide format, one row per participant, with columns:
  - `BI` — breakthrough-infection date (blank/NA if none)
  - `1stVaccine`, `2ndVaccine`, `3rdVaccine` — dose dates
  - `Sab1Val`..`Sab9Val` — antibody measurements (BAU/mL)
  - `Sab1Date`..`Sab9Date` — measurement dates
- **`sample2`** — fitted model object (`.rds`) produced by `01_fit_nlme_model.R`.

> Note for `04_taub_fixed_refit.R`: the post-infection time-to-peak is fixed via `fix()`,
> which requires a literal value in the model block. The script prints the full-model
> estimate on load; set the `fix()` value in the reduced model to that value before running.

## Dependencies

- R 4.5.2
- nlmixr2 5.0.0
- rxode2 5.0.1
- dplyr, tidyr, readr

## Data availability

Individual-level participant data are not publicly available due to the terms of the
cohort consent. These scripts are provided so that the modeling framework and analyses
are reproducible by other groups on comparable data.
