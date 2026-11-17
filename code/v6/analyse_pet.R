## PET / AIPW-PET versus ROW and the other estimators, on matched replicates (identical datasets)
ATE9 <- c("GCOMP", "AIPW", "AIPW-T1/99", "AIPW-T5/95", "ROW", "SU", "TMLE", "AEB-main", "EB-main")
load_main <- function(s) {
  a <- readRDS(sprintf("results_paired/raw_s%02d.rds", s)); a <- a[a$method %in% c("OW", "GCOMP", "AIPW", "AIPW-T1/99"), ]
  b <- readRDS(sprintf("results_balcap/raw_s%02d.rds", s)); b <- b[b$method == "ROW-balcap8", names(a)]; b$method <- "ROW"
  c <- readRDS(sprintf("results_comp/raw_s%02d.rds", s)); c <- c[c$method != "OW", names(a)]
  e <- readRDS(sprintf("results_arow/raw_s%02d.rds", s)); e <- e[e$method == "AIPW-T5/95", names(a)]
  x <- rbind(a, b, c, e); tr <- read.csv("true_estimands.csv")
  k <- match(paste(x$K, x$overlap, x$positivity, x$arm), paste(tr$K, tr$overlap, tr$positivity, tr$arm))
  x$truth <- ifelse(x$effect == "hom", tr$ATE_hom[k], tr$ATE_het[k]); x$scenario <- s
  p <- readRDS(sprintf("results_pet/main/raw_s%02d.rds", s))
  list(x = x[, c("scenario", "method", "om", "effect", "arm", "rep", "est", "truth")], p = p)
}
load_other <- function(d, s) {
  D <- if (d == "heldout") "results_heldout" else "results_stress"; tr <- read.csv(file.path(D, "truth.csv"))
  x <- readRDS(file.path(D, sprintf("raw_s%02d.rds", s))); t <- tr[tr$scenario == s, ]
  x$truth <- as.numeric(t[cbind(match(x$arm, t$arm), match(x$effect, names(t)))])
  p <- readRDS(sprintf("results_pet/%s/raw_s%02d.rds", d, s))
  list(x = x[, c("scenario", "method", "om", "effect", "arm", "rep", "est", "truth")], p = p)
}
ALL <- list(); CHK <- list(); H <- list()
for (d in c("main", "heldout", "stress")) {
  ns <- c(main = 24, heldout = 8, stress = 16)[d]
  for (s in seq_len(ns)) {
    f <- sprintf("results_pet/%s/raw_s%02d.rds", d, s); if (!file.exists(f)) next
    L <- if (d == "main") load_main(s) else load_other(d, s); x <- L$x; p <- L$p
    reps <- unique(p$rep); x <- x[x$rep %in% reps, ]
    g <- p[p$method == "GCOMP-check", ]; g2 <- merge(g, x[x$method == "GCOMP", ], by = c("om", "effect", "arm", "rep"))
    CHK[[length(CHK) + 1]] <- data.frame(design = d, scenario = s, maxdiff = max(abs(g2$est.x - g2$est.y)), n = nrow(g2))
    tr <- unique(x[, c("effect", "arm", "truth")]); p <- merge(p[p$method != "GCOMP-check", ], tr, by = c("effect", "arm"))
    z <- rbind(x, p[, names(x)])
    nm <- z[z$om == "none", ]; oms <- setdiff(unique(z$om), "none")
    z <- rbind(z[z$om != "none", ], do.call(rbind, lapply(oms, function(o) transform(nm, om = o))))
    z$se2 <- (z$est - z$truth)^2; z$design <- d
    ALL[[length(ALL) + 1]] <- aggregate(se2 ~ design + scenario + method + om + effect, data = z, FUN = mean)
    ref <- z[z$method == "ROW", c("effect", "om", "arm", "rep", "se2")]; names(ref)[5] <- "r2"
    for (m in c("PET", "AIPW-PET")) { y <- merge(z[z$method == m, ], ref, by = c("effect", "om", "arm", "rep"))
      for (gg in split(y, list(y$effect, y$om), drop = TRUE)) { dd <- tapply(gg$r2 - gg$se2, gg$rep, mean)
        H[[length(H) + 1]] <- data.frame(design = d, scenario = s, method = m, effect = gg$effect[1], om = gg$om[1],
          ratio = sqrt(mean(gg$r2) / mean(gg$se2)), z = mean(dd) / (sd(dd) / sqrt(length(dd)))) } }
  }
}
P <- do.call(rbind, ALL); P$rmse <- sqrt(P$se2); CHK <- do.call(rbind, CHK); H <- do.call(rbind, H)
cat("Dataset identity check (max |GCOMP difference|):", signif(max(CHK$maxdiff), 3), "over", sum(CHK$n), "estimates\n")
for (refset in list(ATE9, c(ATE9, "PET", "AIPW-PET"))) {
  best <- aggregate(rmse ~ design + scenario + om + effect, data = P[P$method %in% refset, ], FUN = min); names(best)[5] <- "best"
  Q <- merge(P, best); Q$reg <- Q$rmse / Q$best
  W <- aggregate(reg ~ method + design, data = Q, FUN = max); W2 <- aggregate(reg ~ method, data = Q, FUN = max)
  Md <- aggregate(reg ~ method + design, data = Q, FUN = median)
  cat("\n== reference set:", if (length(refset) == 9) "nine ATE estimators" else "nine + PET + AIPW-PET", "==\n")
  w <- reshape(W, idvar = "method", timevar = "design", direction = "wide"); w$all <- W2$reg[match(w$method, W2$method)]
  m <- reshape(Md, idvar = "method", timevar = "design", direction = "wide"); w <- merge(w, m, by = "method", suffixes = c("", ".med"))
  print(w[order(w$all), ], digits = 3, row.names = FALSE)
  if (length(refset) == 9) write.csv(w, "results_pet/regret_pet_ref9.csv", row.names = FALSE) else write.csv(w, "results_pet/regret_pet_ref11.csv", row.names = FALSE)
}
HS <- do.call(rbind, lapply(split(H, list(H$method, H$design), drop = TRUE), function(g) data.frame(method = g$method[1], design = g$design[1],
  n = nrow(g), row_better = sum(g$ratio < 1), sig_better = sum(g$z < -1.96), sig_worse = sum(g$z > 1.96), median_ratio = round(median(g$ratio), 3))))
print(HS, row.names = FALSE); write.csv(HS, "results_pet/h2h_row_vs_pet.csv", row.names = FALSE); write.csv(CHK, "results_pet/identity_check.csv", row.names = FALSE)
best <- aggregate(rmse ~ design + scenario + om + effect, data = P[P$method %in% ATE9, ], FUN = min); names(best)[5] <- "best"
Q <- merge(P, best); Q$reg <- Q$rmse / Q$best; q <- Q[Q$design == "stress" & Q$method == "ROW", ]; w <- q[which.max(q$reg), ]; print(w)
print(Q[Q$design == "stress" & Q$scenario == w$scenario & Q$om == w$om & Q$effect == w$effect, c("method", "rmse", "reg")])
