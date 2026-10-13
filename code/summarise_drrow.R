############################################################################
# Summarise DR-ROW / A-ROW extension
# Usage: Rscript summarise_drrow.R results_drrow results results_row
############################################################################
args <- commandArgs(trailingOnly = TRUE)
DD <- if (length(args) >= 1) args[1] else "results_drrow"
MD <- if (length(args) >= 2) args[2] else "results"
RD <- if (length(args) >= 3) args[3] else "results_row"

raw <- do.call(rbind, lapply(list.files(DD, "^raw_s\\d+\\.rds$", full.names = TRUE), readRDS))
tm  <- read.csv(file.path(MD, "true_estimands.csv"))
tr  <- read.csv(file.path(RD, "true_row_estimands.csv"))

k1 <- paste(raw$K, raw$overlap, raw$positivity, raw$arm)
m  <- match(k1, paste(tm$K, tm$overlap, tm$positivity, tm$arm))
raw$ATE <- ifelse(raw$effect == "hom", tm$ATE_hom[m], tm$ATE_het[m])
ATO     <- ifelse(raw$effect == "hom", tm$ATO_hom[m], tm$ATO_het[m])
cc      <- sub("^A-ROW-c|^DR-ROW-c", "", raw$method)
mr  <- match(paste(k1, paste0("ROW-c", cc)), paste(tr$K, tr$overlap, tr$positivity, tr$arm, tr$method))
ROWT <- ifelse(raw$effect == "hom", tr$T_hom[mr], tr$T_het[mr])
raw$target <- ifelse(grepl("^DR", raw$method) | cc == "Inf", raw$ATE, ifelse(cc == "0", ATO, ROWT))

grp <- c("scenario", "n", "K", "overlap", "positivity", "method", "om", "effect", "arm")
sp  <- split(raw, interaction(raw[, grp], drop = TRUE, lex.order = TRUE))
res <- do.call(rbind, lapply(sp, function(d) {
  z <- qnorm(0.975); lo <- d$est - z * d$se; hi <- d$est + z * d$se
  data.frame(d[1, grp], target = d$target[1], ATE = d$ATE[1],
    abs_bias = abs(mean(d$est) - d$target[1]), abs_bias_ATE = abs(mean(d$est) - d$ATE[1]),
    emp_se = sd(d$est), se_ratio = mean(d$se) / sd(d$est),
    rmse = sqrt(mean((d$est - d$target[1])^2)), rmse_ATE = sqrt(mean((d$est - d$ATE[1])^2)),
    cov = mean(lo <= d$target[1] & hi >= d$target[1]), cov_ATE = mean(lo <= d$ATE[1] & hi >= d$ATE[1]))
}))
agg <- aggregate(cbind(abs_bias, abs_bias_ATE, emp_se, se_ratio, rmse, rmse_ATE, cov, cov_ATE) ~
                 scenario + n + K + overlap + positivity + method + om + effect, res, mean)
write.csv(agg, file.path(DD, "drrow_performance_by_scenario.csv"), row.names = FALSE)
cat("rows:", nrow(raw), "\n")
