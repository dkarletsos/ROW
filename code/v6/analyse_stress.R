D <- "results_stress"; tr <- read.csv(file.path(D, "truth.csv"))
grid <- read.csv(file.path(D, "grid.csv"))
P <- list(); H <- list()
for (s in grid$scenario) {
  f <- file.path(D, sprintf("raw_s%02d.rds", s)); if (!file.exists(f)) next
  x <- readRDS(f); t <- tr[tr$scenario == s, ]
  k <- match(x$arm, t$arm); x$ate <- t[cbind(k, match(x$effect, names(t)))]; x$ate <- as.numeric(x$ate)
  x$se2 <- (x$est - x$ate)^2
  nm <- x[x$om == "none", ]; x <- rbind(x[x$om != "none", ], transform(nm, om = "lin"), transform(nm, om = "flex"))
  P[[s]] <- aggregate(se2 ~ scenario + method + om + effect, data = x, FUN = mean)
  ref <- x[x$method == "ROW", c("effect", "om", "arm", "rep", "se2")]; names(ref)[5] <- "r2"
  for (m in setdiff(unique(x$method), "ROW")) {
    y <- merge(x[x$method == m, ], ref, by = c("effect", "om", "arm", "rep"))
    for (g in split(y, list(y$effect, y$om), drop = TRUE)) { dd <- tapply(g$r2 - g$se2, g$rep, mean)
      H[[length(H) + 1]] <- data.frame(scenario = s, method = m, effect = g$effect[1], om = g$om[1],
        ratio = sqrt(mean(g$r2) / mean(g$se2)), z = mean(dd) / (sd(dd) / sqrt(length(dd)))) } }
}
P <- do.call(rbind, P); H <- do.call(rbind, H); P$rmse <- sqrt(P$se2)
P$correct <- FALSE
ATE <- c("GCOMP", "AIPW", "AIPW-T1/99", "AIPW-T5/95", "ROW", "SU", "TMLE", "AEB-main", "EB-main")
best <- aggregate(rmse ~ scenario + om + effect, data = P[P$method %in% ATE, ], FUN = min); names(best)[4] <- "best"
P <- merge(P, best); P$reg <- P$rmse / P$best
RS <- do.call(rbind, lapply(split(P, P$method), function(g) data.frame(method = g$method[1],
  worst_all = max(g$reg), worst_mis = max(g$reg[!g$correct]),
  worst_binary = max(g$reg[g$effect == "bin_het"]), median = median(g$reg))))
RS <- RS[order(RS$worst_all), ]
HS <- do.call(rbind, lapply(split(H, H$method), function(g) data.frame(method = g$method[1], n = nrow(g),
  row_better = sum(g$ratio < 1), row_sig_better = sum(g$z < -1.96), row_sig_worse = sum(g$z > 1.96), median_ratio = median(g$ratio),
  het_row_better = sum(g$ratio[g$effect != "cont_hom"] < 1), het_n = sum(g$effect != "cont_hom"),
  het_median_ratio = median(g$ratio[g$effect != "cont_hom"]))))
write.csv(P, file.path(D, "perf.csv"), row.names = FALSE); write.csv(RS, file.path(D, "regret.csv"), row.names = FALSE)
write.csv(HS, file.path(D, "h2h.csv"), row.names = FALSE)
cat("Scenarios analysed:", length(unique(P$scenario)), "\n"); print(RS, digits = 3, row.names = FALSE); print(HS, digits = 3, row.names = FALSE)
G <- read.csv(file.path(D, "grid.csv")); P <- merge(P, G[, c("scenario", "level", "psmis", "n", "K")])
cat("\n== worst-case regret by stress factor ==\n")
for (v in c("level", "psmis", "n", "effect")) { cat("--", v, "\n")
  print(round(tapply(P$reg, list(P$method, P[[v]]), max), 2)) }
