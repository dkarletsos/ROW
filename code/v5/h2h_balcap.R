# Paired comparison of adaptive-lambda ROW with truncated AIPW, fixed ROW (2K^2), and OW
tr <- read.csv("true_estimands.csv"); out <- list(); perf <- list(); lam <- list()
for (s in 1:24) {
  a <- readRDS(sprintf("results_paired/raw_s%02d.rds", s)); b <- readRDS(sprintf("results_arow/raw_s%02d.rds", s))
  c <- readRDS(sprintf("results_balcap/raw_s%02d.rds", s))
  x <- rbind(a[a$method %in% c("AIPW-T1/99", "OW", "GCOMP", "AIPW"), names(a)], b[b$method == "AROW-c2", names(a)],
             c[c$method %in% c("ROW-balcap4", "ROW-balcap8", "ROW-balcap16"), names(a)])
  lam[[s]] <- c[c$method != "OW", c("scenario", "n", "K", "overlap", "positivity", "method", "om", "effect", "lam")]
  k <- match(paste(x$K, x$overlap, x$positivity, x$arm), paste(tr$K, tr$overlap, tr$positivity, tr$arm))
  x$ate <- ifelse(x$effect == "hom", tr$ATE_hom[k], tr$ATE_het[k]); x$se2 <- (x$est - x$ate)^2
  x$cov <- abs(x$est - x$ate) <= 1.96 * x$se
  perf[[s]] <- aggregate(cbind(se2, cov) ~ scenario + n + K + overlap + positivity + method + om + effect, data = x, FUN = mean)
  for (refm in c("AIPW-T1/99", "OW", "AROW-c2")) {
    ref <- x[x$method == refm, c("effect", "om", "arm", "rep", "se2")]; names(ref)[5] <- "r2"
    for (m in c("ROW-balcap4", "ROW-balcap8", "ROW-balcap16", "AROW-c2")) { if (m == refm) next
      y <- merge(x[x$method == m, ], if (refm == "OW") ref[, -2] else ref, by = if (refm == "OW") c("effect", "arm", "rep") else c("effect", "om", "arm", "rep"))
      for (g in split(y, list(y$effect, y$om), drop = TRUE)) { dd <- tapply(g$se2 - g$r2, g$rep, mean)
        out[[length(out) + 1]] <- data.frame(scenario = s, ref = refm, method = m, effect = g$effect[1], om = g$om[1],
          ratio = sqrt(mean(g$se2) / mean(g$r2)), z = mean(dd) / (sd(dd) / sqrt(length(dd)))) } } }
}
o <- do.call(rbind, out); P <- do.call(rbind, perf); L <- do.call(rbind, lam)
write.csv(o, "results_balcap/h2h.csv", row.names = FALSE); write.csv(P, "results_balcap/perf.csv", row.names = FALSE)
cat("\n== head-to-head (ratio = RMSE method / RMSE ref) ==\n")
print(do.call(rbind, lapply(split(o, list(o$ref, o$method, o$effect, o$om), drop = TRUE), function(g) data.frame(
  ref = g$ref[1], method = g$method[1], effect = g$effect[1], om = g$om[1], better = sum(g$ratio < 1), sig_b = sum(g$z < -1.96),
  sig_w = sum(g$z > 1.96), med = round(median(g$ratio), 3), min = round(min(g$ratio), 3), max = round(max(g$ratio), 3)))), row.names = FALSE)
P$rmse <- sqrt(P$se2); P$correct <- P$om == "quad" | (P$om == "lin" & P$effect == "hom")
A <- P[P$method != "OW", ]; best <- aggregate(rmse ~ scenario + om + effect, data = A[A$method %in% c("GCOMP", "AIPW", "AIPW-T1/99", "AROW-c2", "ROW-balcap4", "ROW-balcap8", "ROW-balcap16"), ], FUN = min)
names(best)[4] <- "best"; A <- merge(A, best); A$reg <- A$rmse / A$best
cat("\n== regret (common best over all ATE estimators) and coverage (IF SE, selection ignored) ==\n")
print(do.call(rbind, lapply(split(A, A$method), function(g) data.frame(method = g$method[1], reg_cor = round(max(g$reg[g$correct]), 3),
  reg_mis = round(max(g$reg[!g$correct]), 3), reg_all = round(max(g$reg), 3), med_reg = round(median(g$reg), 3),
  cov_cor = paste(round(100 * range(g$cov[g$correct])), collapse = "-"), cov_mis = paste(round(100 * range(g$cov[!g$correct])), collapse = "-")))), row.names = FALSE)
cat("\n== chosen lambda/K^2 (share) ==\n"); L$correct <- L$om == "quad" | (L$om == "lin" & L$effect == "hom")
print(round(prop.table(table(paste(L$method, ifelse(L$correct, "correct", "misspec")), L$lam), 1), 2))
