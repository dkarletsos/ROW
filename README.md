# ROW: regularised outcome-weighted estimation of the average treatment effect

`ROW` estimates the average treatment effect (ATE) for comparisons of two or more treatments in
observational data with limited covariate overlap.

- Per-arm outcome models are averaged over all patients and corrected with residuals weighted by
  `w = 1 / (e_A(X) + 1/lambda)`, where `e_A(X)` is the generalised propensity score of the treatment received.
- `lambda` is chosen **without looking at outcomes**: the smallest value (up to 8K²) at which every weighted
  arm matches the full population on the covariates and their squares (max SMD ≤ 0.1).
- Standard errors come from an M-estimation sandwich that accounts for estimating the propensity score
  and the outcome models.
- `lambda → 0` gives outcome regression; `lambda = Inf` gives AIPW.

## Installation

```r
# install.packages("remotes")
remotes::install_github("dkarletsos/ROW")
```

Dependencies: base R and `nnet`.

## Example: education and wages (NCDS)

```r
library(ROW)
data(ncds)
f <- ~ white + maemp + factor(scht) + qmab + qmab2 + qvab + qvab2 +
       paed_u + maed_u + agepa + agema + sib_u
fit <- row_ate(ncds, treatment = "Dmult", outcome = "wage", ps_formula = f, reference = "None")
fit
#>        contrast estimate    se lower upper
#>  >=A/eq vs None    0.279 0.024 0.233 0.325
#>    O/eq vs None    0.113 0.024 0.065 0.160

row_path(ncds, "Dmult", "wage", f, reference = "None")   # estimates along lambda (sensitivity analysis)
```

Main functions: `row_ate()` (estimate, sandwich SE, balance rule), `row_path()` (estimates across lambda),
`row_weights()` (the weights themselves).

## Reference

Karletsos D. *Overlap Weighting, Trimmed IPTW and a Regularised Outcome-Weighted Estimator (ROW) for Multi-Group
Real-World Comparisons: Estimand Choice, Bias and Precision Under Positivity Violations.* PHUSE EU Connect 2026,
paper AS03 (`PHUSE2026_final/PAP_AS03.pdf`).

---

## Reproducing the PHUSE 2026 results

Simulation code is in `code/` (excluded from the package build).

### Paired comparison of estimators (`code/v5`, results in `results_v5`)

    cd code/v5
    Rscript paired_sim.R 2000 2 results_paired        # OW, outcome regression, AIPW, truncated AIPW (same datasets for all methods)
    Rscript balcap_sim.R 2000 2 results_balcap        # ROW with the balance rule (caps 4K², 8K², 16K²)
    Rscript comparators_sim.R 2000 2 results_comp     # augmented OW, entropy balancing, Su et al. shrinkage, TMLE
    Rscript arow_sim.R 2000 2 results_arow            # ROW with fixed lambda; AIPW truncated 5/95
    Rscript analyse_comp.R                            # worst-case regret and head-to-head comparisons
    Rscript row_sandwich_coverage.R 1000 2 results_rowsand
    Rscript ncds_analysis.R                           # real-data illustration

All scripts use the seed 20261115 + scenario number, so every method sees identical datasets.

### Original study (OW vs trimmed IPTW)

    Rscript code/ow_vs_iptw_sim.R 2000 2 results   # 24 scenarios x 2,000 reps; base R + nnet
    Rscript code/summarise_sim.R results
    python3 code/make_figures.py results figures

Raw replicate files (`raw_sNN.rds`) are not included; rerun the simulations to regenerate them.

### Held-out validation, stress test and extended sandwich coverage (`code/v6`, results in `results_v6`)

Run from `code/v6` (scripts source `ow_vs_iptw_sim.R` and `true_estimands.csv` from the same folder; `cover_ext.R`
sources the package code in `R/row.R`).

    Rscript heldout_sim.R 800 2 results_heldout      # independent design: 8 confounders, binary outcome, n = 1,000; cap fixed at 8K^2
    Rscript analyse_heldout.R
    Rscript stress_sim.R 600 2 results_stress        # misspecified GPS and outcome models, t3 errors, very poor overlap
    Rscript analyse_stress.R
    for d in main heldout stress; do Rscript cover_ext.R $d 500 2; done   # sandwich coverage, incl. misspecified outcome models
    Rscript analyse_cover.R
    Rscript supp_tables.R                            # per-scenario bias / SE / RMSE / coverage (supplementary Tables S2-S3)
    python3 make_paired_figure.py ../../results_v5 ../../figures   # Figure 1
    python3 make_fig2.py                             # Figure 2: worst-case regret across the three designs
