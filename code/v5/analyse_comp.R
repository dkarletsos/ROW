# Merge: paired (OW, GCOMP, AIPW, AIPW-T1/99), balcap (ROW = balance rule cap 4K^2), comparators; per-scenario aggregation
tr <- read.csv("true_estimands.csv"); P <- list(); H <- list()
for (s in 1:24) {
  a <- readRDS(sprintf("results_paired/raw_s%02d.rds", s)); a <- a[a$method %in% c("OW", "GCOMP", "AIPW", "AIPW-T1/99"), ]
  b <- readRDS(sprintf("results_balcap/raw_s%02d.rds", s)); b <- b[b$method == "ROW-balcap8", names(a)]; b$method <- "ROW"
  c <- readRDS(sprintf("results_comp/raw_s%02d.rds", s)); c <- c[c$method != "OW", names(a)]
  e <- readRDS(sprintf("results_arow/raw_s%02d.rds", s)); e <- e[e$method == "AIPW-T5/95", names(a)]
  f1 <- readRDS(sprintf("results_balcap/raw_s%02d.rds", s)); f1 <- f1[f1$method %in% c("ROW-balcap4", "ROW-balcap16"), names(a)]
  f2 <- readRDS(sprintf("results_arow/raw_s%02d.rds", s)); f2 <- f2[f2$method == "AROW-c2", names(a)]; f2$method <- "ROW-fixed2"
  f3 <- readRDS(sprintf("results_balance/raw_s%02d.rds", s)); f3 <- f3[f3$method == "ROW-bal-sq-10", names(a)]; f3$method <- "ROW-nocap"
  f4 <- readRDS(sprintf("results_adaptive/raw_s%02d.rds", s)); f4 <- f4[f4$method == "ROW-mse", names(a)]
  x <- rbind(a, b, c, e, f1, f2, f3, f4)
  k <- match(paste(x$K, x$overlap, x$positivity, x$arm), paste(tr$K, tr$overlap, tr$positivity, tr$arm))
  x$ate <- ifelse(x$effect == "hom", tr$ATE_hom[k], tr$ATE_het[k]); x$se2 <- (x$est - x$ate)^2
  x$cov <- ifelse(is.na(x$se), NA, abs(x$est - x$ate) <= 1.96 * x$se)
  # methods without an outcome model (om = none) are replicated for both outcome models
  nm <- x[x$om == "none", ]; x <- rbind(x[x$om != "none", ], transform(nm, om = "lin"), transform(nm, om = "quad"))
  P[[s]] <- aggregate(cbind(se2, cov) ~ scenario + n + K + overlap + positivity + method + om + effect, data = x, FUN = mean, na.action = na.pass)
  ref <- x[x$method == "ROW", c("effect", "om", "arm", "rep", "se2")]; names(ref)[5] <- "r2"
  for (m in setdiff(unique(x$method), "ROW")) {
    y <- merge(x[x$method == m, ], ref, by = c("effect", "om", "arm", "rep"))
    for (g in split(y, list(y$effect, y$om), drop = TRUE)) { dd <- tapply(g$r2 - g$se2, g$rep, mean)
      H[[length(H) + 1]] <- data.frame(scenario = s, method = m, effect = g$effect[1], om = g$om[1],
        ratio_row = sqrt(mean(g$r2) / mean(g$se2)), z = mean(dd) / (sd(dd) / sqrt(length(dd)))) } }
}
P <- do.call(rbind, P); H <- do.call(rbind, H); P$rmse <- sqrt(P$se2)
P$correct <- P$om == "quad" | (P$om == "lin" & P$effect == "hom")
ATE <- c("GCOMP", "AIPW", "AIPW-T1/99", "AIPW-T5/95", "ROW", "SU", "TMLE", "AEB-main", "AEB", "EB-main", "EB")
best <- aggregate(rmse ~ scenario + om + effect, data = P[P$method %in% setdiff(ATE, c("EB", "AEB")), ], FUN = min); names(best)[4] <- "best"
P <- merge(P, best); P$reg <- P$rmse / P$best
write.csv(P, "results_comp/perf_all.csv", row.names = FALSE); write.csv(H, "results_comp/h2h_vs_row.csv", row.names = FALSE)
RS <- do.call(rbind, lapply(split(P, P$method), function(g) data.frame(method = g$method[1],
  reg_cor = max(g$reg[g$correct]), reg_mis = max(g$reg[!g$correct]), reg_all = max(g$reg), med = median(g$reg),
  cov_cor_min = min(g$cov[g$correct], na.rm = TRUE), cov_cor_max = max(g$cov[g$correct], na.rm = TRUE),
  cov_mis_min = min(g$cov[!g$correct], na.rm = TRUE), cov_mis_max = max(g$cov[!g$correct], na.rm = TRUE))))
write.csv(RS, "results_comp/regret_all.csv", row.names = FALSE)
cat("== worst-case regret (best among ATE estimators excluding oracle-feature EB/AEB) ==\n")
print(do.call(rbind, lapply(split(P, P$method), function(g) data.frame(method = g$method[1],
  reg_cor = round(max(g$reg[g$correct]), 2), reg_mis = round(max(g$reg[!g$correct]), 2), reg_all = round(max(g$reg), 2),
  med = round(median(g$reg), 3), cov_cor = if (all(is.na(g$cov))) "-" else paste(round(100 * range(g$cov[g$correct], na.rm = TRUE)), collapse = "-"),
  cov_mis = if (all(is.na(g$cov))) "-" else paste(round(100 * range(g$cov[!g$correct], na.rm = TRUE)), collapse = "-")))), row.names = FALSE)
cat("\n== ROW vs each method (ratio = RMSE ROW / RMSE method; <1 favours ROW) ==\n")
H$correct <- H$om == "quad" | (H$om == "lin" & H$effect == "hom")
print(do.call(rbind, lapply(split(H, list(H$method, H$correct), drop = TRUE), function(g) data.frame(method = g$method[1],
  model = ifelse(g$correct[1], "correct", "misspec"), row_better = sum(g$ratio_row < 1), sig_row_better = sum(g$z < -1.96),
  sig_row_worse = sum(g$z > 1.96), med = round(median(g$ratio_row), 3), min = round(min(g$ratio_row), 3), max = round(max(g$ratio_row), 3)))), row.names = FALSE)
