# Summarise g-computation comparator against the true ATE
args <- commandArgs(trailingOnly = TRUE)
GD <- if (length(args) >= 1) args[1] else "results_gcomp"
MD <- if (length(args) >= 2) args[2] else "results"
raw <- do.call(rbind, lapply(list.files(GD, "^raw_s\\d+\\.rds$", full.names = TRUE), readRDS))
tm  <- read.csv(file.path(MD, "true_estimands.csv"))
m <- match(paste(raw$K, raw$overlap, raw$positivity, raw$arm), paste(tm$K, tm$overlap, tm$positivity, tm$arm))
raw$ATE <- ifelse(raw$effect == "hom", tm$ATE_hom[m], tm$ATE_het[m])
grp <- c("scenario", "n", "K", "overlap", "positivity", "method", "om", "effect", "arm")
sp <- split(raw, interaction(raw[, grp], drop = TRUE, lex.order = TRUE))
res <- do.call(rbind, lapply(sp, function(d) {
  lo <- d$est - 1.96 * d$se; hi <- d$est + 1.96 * d$se
  data.frame(d[1, grp], abs_bias = abs(mean(d$est) - d$ATE[1]), emp_se = sd(d$est),
             se_ratio = mean(d$se) / sd(d$est), rmse = sqrt(mean((d$est - d$ATE[1])^2)),
             cov = mean(lo <= d$ATE[1] & hi >= d$ATE[1]))
}))
agg <- aggregate(cbind(abs_bias, emp_se, se_ratio, rmse, cov) ~ scenario + n + K + overlap + positivity + method + om + effect, res, mean)
write.csv(agg, file.path(GD, "gcomp_performance_by_scenario.csv"), row.names = FALSE)
cat("rows:", nrow(raw), "\n")
