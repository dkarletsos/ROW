res <- list()
for (d in c("main", "heldout", "stress")) for (f in list.files(file.path("results_cover", d), full.names = TRUE)) {
  x <- readRDS(f); x$design <- d
  g <- do.call(rbind, lapply(split(x, list(x$method, x$effect, x$om, x$arm), drop = TRUE), function(z) data.frame(
    design = d, scenario = z$scenario[1], method = z$method[1], effect = z$effect[1], om = z$om[1], arm = z$arm[1],
    cov = mean(abs(z$est - z$truth) <= 1.96 * z$se), cov_naive = mean(abs(z$est - z$truth) <= 1.96 * z$se_naive),
    se_ratio = mean(z$se) / sd(z$est), bias_sd = abs(mean(z$est - z$truth)) / sd(z$est), nrep = nrow(z))))
  res[[f]] <- aggregate(cbind(cov, cov_naive, se_ratio, bias_sd) ~ design + scenario + method + effect + om, data = g, FUN = mean)
}
R <- do.call(rbind, res)
R$correct <- (R$effect == "cont_hom") | (R$design %in% c("main", "heldout") & R$effect == "cont_het" & R$om == "flex")
write.csv(R, "results_cover/coverage_by_cell.csv", row.names = FALSE)
S <- do.call(rbind, lapply(split(R, list(R$design, R$correct, R$method), drop = TRUE), function(z) data.frame(
  design = z$design[1], outcome_model = ifelse(z$correct[1], "correct", "misspecified"), method = z$method[1], cells = nrow(z),
  cov_min = min(z$cov), cov_med = median(z$cov), cov_max = max(z$cov), naive_med = median(z$cov_naive),
  se_ratio_min = min(z$se_ratio), se_ratio_med = median(z$se_ratio), se_ratio_max = max(z$se_ratio), bias_sd_med = median(z$bias_sd))))
S <- S[order(S$design, S$outcome_model, S$method), ]
write.csv(S, "results_cover/coverage_summary.csv", row.names = FALSE); print(S, digits = 3, row.names = FALSE)
# coverage among cells where |bias| < 0.25 SD (bias-free inference check)
B <- R[R$bias_sd < 0.25, ]
print(aggregate(cbind(cov, se_ratio) ~ method + design, data = B, FUN = function(v) round(c(min = min(v), med = median(v), n = length(v)), 3)))
