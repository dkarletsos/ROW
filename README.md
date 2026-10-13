# Overlap weighting vs trimmed IPTW: PHUSE EU Connect 2026 simulation

Reproduces every number in the paper and slides.

    Rscript code/ow_vs_iptw_sim.R 2000 2 results   # 24 scenarios x 2,000 reps (~25 min, 2 cores); base R + nnet
    Rscript code/summarise_sim.R results           # performance_by_contrast.csv, performance_by_scenario.csv
    python3 code/make_figures.py results figures   # matplotlib, pandas
    python3 code/make_framework.py figures

Each scenario uses its own seed (20261115 + scenario number), so any single scenario can be rerun on its own.
Raw replicate files (raw_sNN.rds) are not included; rerun the simulation to regenerate them.

## ROW (Regularised Overlap Weights) extension
    cd code && Rscript row_extension.R 2000 2 ../results_row   # same seeds/datasets as main study (~31 min)
    Rscript summarise_row.R ../results_row ../results
    python3 make_row_figure.py ../results_row ../figures
ROW: h_lambda(X) = 1 / (sum_k 1/e_k(X) + lambda), w = h_lambda / e_A; lambda = c * K^2, c in {0.5, 2, 8, 32}.

## DR-ROW / A-ROW (augmented ROW) extension
    cd code && Rscript drrow_extension.R 2000 2 ../results_drrow   # ~58 min on 2 cores; same datasets as main study
    Rscript summarise_drrow.R ../results_drrow ../results ../results_row
    python3 make_drrow_figure.py ../results_drrow ../results ../figures
DR-ROW (ATE): mu_k = mean(m_k(X)) + Hajek mean of residuals in arm k with ROW weights h_lambda/e_k.
A-ROW (own estimand): plug-in averaged over the h_lambda-tilted population + same residual correction.
Outcome models: lin (X1..X5) and quad (+ X1^2, X2^2, X3^2).

## Outcome-regression (g-computation) comparator
    cd code && Rscript gcomp_extension.R 2000 2 ../results_gcomp   # ~8 min
    Rscript summarise_gcomp.R ../results_gcomp ../results
regret_table.csv: worst-case / mean RMSE relative to the best ATE estimator per scenario x outcome model (96 cases).
