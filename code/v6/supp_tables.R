tr <- read.csv("true_estimands.csv"); out <- list()
for (s in 1:24) {
  a <- readRDS(sprintf("results_paired/raw_s%02d.rds", s)); a <- a[a$method == "OW", ]
  b <- readRDS(sprintf("results_balcap/raw_s%02d.rds", s)); b <- b[b$method == "ROW-balcap8", names(a)]; b$method <- "ROW"
  x <- rbind(a, b)
  k <- match(paste(x$K, x$overlap, x$positivity, x$arm), paste(tr$K, tr$overlap, tr$positivity, tr$arm))
  x$ate <- ifelse(x$effect == "hom", tr$ATE_hom[k], tr$ATE_het[k])
  x$cov <- ifelse(is.na(x$se), NA, abs(x$est - x$ate) <= 1.96 * x$se)
  g <- aggregate(cbind(bias = est - ate, se2 = (est - ate)^2, cov) ~ scenario + n + K + overlap + positivity + method + om + effect + arm, data = x, FUN = mean, na.action = na.pass)
  e <- aggregate(est ~ scenario + method + om + effect + arm, data = x, FUN = sd); g <- merge(g, e)
  out[[s]] <- aggregate(cbind(absbias = abs(bias), empse = est, rmse = sqrt(se2), cov) ~ scenario + n + K + overlap + positivity + method + om + effect, data = g, FUN = mean, na.action = na.pass)
}
P <- do.call(rbind, out); write.csv(P, "results_comp/supp_per_scenario.csv", row.names = FALSE); print(head(P)); print(nrow(P))
