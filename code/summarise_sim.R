############################################################################
# Summarise simulation output: bias, empirical SE, RMSE, coverage, ESS, SMD
# Usage: Rscript summarise_sim.R results
############################################################################
args <- commandArgs(trailingOnly = TRUE)
OUTDIR <- if (length(args) >= 1) args[1] else "results"

raw <- do.call(rbind, lapply(list.files(OUTDIR, "^raw_s\\d+\\.rds$", full.names = TRUE), readRDS))
tr  <- read.csv(file.path(OUTDIR, "true_estimands.csv"))

# attach the true ATE and ATO for every row
key <- function(d) paste(d$K, d$overlap, d$positivity, d$arm)
m <- match(key(raw), key(tr))
raw$ATE <- ifelse(raw$effect == "hom", tr$ATE_hom[m], tr$ATE_het[m])
raw$ATO <- ifelse(raw$effect == "hom", tr$ATO_hom[m], tr$ATO_het[m])
# each method's own (nominal) target: IPTW variants -> ATE; OW -> ATO
raw$target <- ifelse(raw$method == "OW", raw$ATO, raw$ATE)

grp <- c("scenario", "n", "K", "overlap", "positivity", "method", "effect", "arm")
sp  <- split(raw, interaction(raw[, grp], drop = TRUE, lex.order = TRUE))
res <- do.call(rbind, lapply(sp, function(d) {
  lo <- d$est - qnorm(0.975) * d$se; hi <- d$est + qnorm(0.975) * d$se
  R  <- nrow(d)
  data.frame(d[1, grp], n_rep = R,
    true_ATE = d$ATE[1], true_ATO = d$ATO[1], target = d$target[1],
    mean_est   = mean(d$est),
    bias       = mean(d$est) - d$target[1],
    bias_vs_ATE = mean(d$est) - d$ATE[1],
    emp_se     = sd(d$est),
    rmse       = sqrt(mean((d$est - d$target[1])^2)),
    rmse_vs_ATE = sqrt(mean((d$est - d$ATE[1])^2)),
    mod_se     = mean(d$se),
    se_ratio   = mean(d$se) / sd(d$est),
    coverage   = mean(lo <= d$target[1] & hi >= d$target[1]),
    coverage_vs_ATE = mean(lo <= d$ATE[1] & hi >= d$ATE[1]),
    mcse_bias  = sd(d$est) / sqrt(R),
    mcse_cov   = sqrt(mean(lo <= d$target[1] & hi >= d$target[1]) *
                      (1 - mean(lo <= d$target[1] & hi >= d$target[1])) / R),
    ess_pct    = mean(d$ess), max_smd = mean(d$smd),
    pct_smd_gt_0.1 = 100 * mean(d$smd > 0.1),
    wmax_ratio = median(d$wmax_ratio))
}))
rownames(res) <- NULL
res <- res[order(res$scenario, res$effect, res$arm, res$method), ]
write.csv(res, file.path(OUTDIR, "performance_by_contrast.csv"), row.names = FALSE)

# Scenario-level summary averaged over the K-1 contrasts
agg_cols <- c("bias", "bias_vs_ATE", "emp_se", "rmse", "rmse_vs_ATE", "se_ratio",
              "coverage", "coverage_vs_ATE", "ess_pct", "max_smd", "pct_smd_gt_0.1")
res$abs_bias <- abs(res$bias); res$abs_bias_vs_ATE <- abs(res$bias_vs_ATE)
sp2 <- split(res, interaction(res[, grp[-8]], drop = TRUE, lex.order = TRUE))
agg <- do.call(rbind, lapply(sp2, function(d) data.frame(d[1, grp[-8]],
  mean_abs_bias = mean(d$abs_bias), mean_abs_bias_vs_ATE = mean(d$abs_bias_vs_ATE),
  t(colMeans(d[, setdiff(agg_cols, c("bias", "bias_vs_ATE"))])))))
rownames(agg) <- NULL
write.csv(agg, file.path(OUTDIR, "performance_by_scenario.csv"), row.names = FALSE)
cat("Rows:", nrow(raw), " contrasts summarised:", nrow(res), "\n")
