############################################################################
# Summarise the paired simulation: performance vs ATE (and OW vs ATO),
# paired RMSE ratios ROW / OW with Monte Carlo significance of the MSE difference.
# Usage: Rscript summarise_paired.R [rawdir] [truths.csv] [outdir]
############################################################################
args  <- commandArgs(trailingOnly = TRUE)
RAW   <- if (length(args) >= 1) args[1] else "results_paired"
TRUTH <- if (length(args) >= 2) args[2] else "true_estimands.csv"
OUT   <- if (length(args) >= 3) args[3] else RAW

raw <- do.call(rbind, lapply(list.files(RAW, "^raw_s\\d+\\.rds$", full.names = TRUE), readRDS))
tr  <- read.csv(TRUTH)
key <- function(d) paste(d$K, d$overlap, d$positivity, d$arm)
m   <- match(key(raw), key(tr))
raw$ATE <- ifelse(raw$effect == "hom", tr$ATE_hom[m], tr$ATE_het[m])
raw$ATO <- ifelse(raw$effect == "hom", tr$ATO_hom[m], tr$ATO_het[m])
raw$err <- raw$est - raw$ATE
raw$cov <- abs(raw$err) <= qnorm(0.975) * raw$se
raw$cov_ATO <- abs(raw$est - raw$ATO) <= qnorm(0.975) * raw$se

grp <- c("scenario", "n", "K", "overlap", "positivity", "method", "om", "effect")
# per arm, then averaged over arms
pa <- aggregate(cbind(err, cov, cov_ATO, se, est, ATO, ATE) ~ scenario + n + K + overlap + positivity + method + om + effect + arm,
                data = raw, FUN = mean)
pa_sd  <- aggregate(est ~ scenario + n + K + overlap + positivity + method + om + effect + arm, data = raw, FUN = sd)
pa_mse <- aggregate(I(err^2) ~ scenario + n + K + overlap + positivity + method + om + effect + arm, data = raw, FUN = mean)
pa_cnt <- aggregate(est ~ scenario + n + K + overlap + positivity + method + om + effect + arm, data = raw, FUN = length)
pa$emp_se <- pa_sd$est; pa$rmse <- sqrt(pa_mse[[ncol(pa_mse)]]); pa$nrep <- pa_cnt$est
pa$abs_bias <- abs(pa$err); pa$abs_bias_ATO <- abs(pa$est - pa$ATO); pa$se_ratio <- pa$se / pa$emp_se
perf <- aggregate(cbind(abs_bias, emp_se, rmse, se_ratio, cov, abs_bias_ATO, cov_ATO, nrep) ~
                    scenario + n + K + overlap + positivity + method + om + effect, data = pa, FUN = mean)
perf <- perf[order(perf$scenario, perf$effect, perf$method, perf$om), ]
write.csv(perf, file.path(OUT, "paired_performance.csv"), row.names = FALSE)

# Paired comparisons against OW (same replicate, same arm)
ow <- raw[raw$method == "OW", c("scenario", "effect", "arm", "rep", "err")]
names(ow)[5] <- "err_ow"
cmp <- merge(raw[raw$method != "OW", ], ow, by = c("scenario", "effect", "arm", "rep"))
cmp$d <- cmp$err^2 - cmp$err_ow^2                       # paired squared-error difference
pc <- do.call(rbind, lapply(split(cmp, list(cmp$scenario, cmp$method, cmp$om, cmp$effect), drop = TRUE), function(x) {
  # average over arms within replicate first, then MC SE across replicates
  dr <- tapply(x$d, x$rep, mean)
  data.frame(x[1, c("scenario", "n", "K", "overlap", "positivity", "method", "om", "effect")],
             rmse = sqrt(mean(x$err^2)), rmse_ow = sqrt(mean(x$err_ow^2)),
             ratio = sqrt(mean(x$err^2) / mean(x$err_ow^2)),
             z = mean(dr) / (sd(dr) / sqrt(length(dr))), row.names = NULL)
}))
pc <- pc[order(pc$method, pc$om, pc$effect, pc$scenario), ]
write.csv(pc, file.path(OUT, "paired_vs_ow.csv"), row.names = FALSE)

# Regret (worst-case RMSE ratio vs best ATE estimator in each scenario x om x effect), excluding OW
ate <- perf[perf$method != "OW", ]
best <- aggregate(rmse ~ scenario + om + effect, data = ate, FUN = min); names(best)[4] <- "best"
ate <- merge(ate, best); ate$regret <- ate$rmse / ate$best
ate$correct <- with(ate, (om == "quad") | (om == "lin" & effect == "hom"))
reg <- do.call(rbind, lapply(split(ate, ate$method), function(x) data.frame(method = x$method[1],
  reg_cor_max = max(x$regret[x$correct]), reg_mis_max = max(x$regret[!x$correct]), reg_max = max(x$regret),
  cov_cor_min = min(x$cov[x$correct]), cov_cor_max = max(x$cov[x$correct]),
  cov_mis_min = min(x$cov[!x$correct]), cov_mis_max = max(x$cov[!x$correct]))))
write.csv(reg, file.path(OUT, "paired_regret.csv"), row.names = FALSE)
print(reg)
s <- pc[pc$method == "ROW", ]
print(do.call(rbind, lapply(split(s, list(s$effect, s$om), drop = TRUE), function(x) data.frame(
  effect = x$effect[1], om = x$om[1], wins = sum(x$ratio < 1), sig_wins = sum(x$z < -1.96),
  sig_losses = sum(x$z > 1.96), med = median(x$ratio), min = min(x$ratio), max = max(x$ratio)))))
