############################################################################
# Summarise ROW extension: performance vs own target and vs the ATE
# Usage: Rscript summarise_row.R results_row results
############################################################################
args <- commandArgs(trailingOnly = TRUE)
RD <- if (length(args) >= 1) args[1] else "results_row"
MD <- if (length(args) >= 2) args[2] else "results"

raw <- do.call(rbind, lapply(list.files(RD, "^raw_s\\d+\\.rds$", full.names = TRUE), readRDS))
tm  <- read.csv(file.path(MD, "true_estimands.csv"))
tr  <- read.csv(file.path(RD, "true_row_estimands.csv"))

k1 <- paste(raw$K, raw$overlap, raw$positivity, raw$arm)
m  <- match(k1, paste(tm$K, tm$overlap, tm$positivity, tm$arm))
raw$ATE <- ifelse(raw$effect == "hom", tm$ATE_hom[m], tm$ATE_het[m])
raw$ATO <- ifelse(raw$effect == "hom", tm$ATO_hom[m], tm$ATO_het[m])
mr <- match(paste(k1, raw$method), paste(tr$K, tr$overlap, tr$positivity, tr$arm, tr$method))
raw$ROWT <- ifelse(raw$effect == "hom", tr$T_hom[mr], tr$T_het[mr])
raw$target <- ifelse(raw$method == "OW", raw$ATO, ifelse(grepl("^ROW", raw$method), raw$ROWT, raw$ATE))

grp <- c("scenario", "n", "K", "overlap", "positivity", "method", "effect", "arm")
sp  <- split(raw, interaction(raw[, grp], drop = TRUE, lex.order = TRUE))
res <- do.call(rbind, lapply(sp, function(d) {
  z <- qnorm(0.975); lo <- d$est - z * d$se; hi <- d$est + z * d$se
  data.frame(d[1, grp], target = d$target[1], ATE = d$ATE[1],
    abs_bias = abs(mean(d$est) - d$target[1]), abs_bias_ATE = abs(mean(d$est) - d$ATE[1]),
    emp_se = sd(d$est),
    rmse = sqrt(mean((d$est - d$target[1])^2)), rmse_ATE = sqrt(mean((d$est - d$ATE[1])^2)),
    cov = mean(lo <= d$target[1] & hi >= d$target[1]), cov_ATE = mean(lo <= d$ATE[1] & hi >= d$ATE[1]),
    ess = mean(d$ess), smd = mean(d$smd), wmax = median(d$wmax_ratio))
}))
agg <- aggregate(cbind(target, abs_bias, abs_bias_ATE, emp_se, rmse, rmse_ATE, cov, cov_ATE, ess, smd, wmax) ~
                 scenario + n + K + overlap + positivity + method + effect, res, mean)
write.csv(agg, file.path(RD, "row_performance_by_scenario.csv"), row.names = FALSE)
cat("rows:", nrow(raw), "\n")
