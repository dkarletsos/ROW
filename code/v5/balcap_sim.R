############################################################################
# Outcome-blind, data-driven lambda for ROW: balance rule.
# For each dataset (before looking at outcomes): choose the smallest lambda on the grid such that,
# in every arm, the ROW-weighted covariate means match the full-population means:
#   max_k max_j |mean_w,k(f_j) - mean_all(f_j)| / sd_all(f_j) <= threshold
# If no lambda meets the threshold, take the lambda with the smallest maximum SMD.
# Features: "main" = X1..X5 ; "sq" = X1..X5 + X1^2, X2^2, X3^2.
# Same seeds as paired_sim.R -> identical datasets.
############################################################################
Sys.setenv(SIM_SOURCE_ONLY = "1")
source("ow_vs_iptw_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_balcap"
SCEN   <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else seq_len(nrow(grid))
dir.create(OUTDIR, showWarnings = FALSE)
CG <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, Inf)
RULES <- list(`ROW-balcap4` = list(f = "sq", t = 0.10, cap = 4), `ROW-balcap8` = list(f = "sq", t = 0.10, cap = 8),
              `ROW-balcap16` = list(f = "sq", t = 0.10, cap = 16))

wfun <- function(ea, cc, K) if (is.infinite(cc)) 1 / ea else 1 / (ea + 1 / (cc * K^2))
maxsmd <- function(Fm, A, w, K) {
  mu <- colMeans(Fm); sdv <- apply(Fm, 2, sd)
  max(sapply(1:K, function(k) { i <- A == k; max(abs(colSums(Fm[i, , drop = FALSE] * w[i]) / sum(w[i]) - mu) / sdv) }))
}
aug <- function(M, Y, A, w, K) {
  n <- length(Y); res <- Y - M[cbind(seq_len(n), A)]; corr <- numeric(K); ICc <- matrix(0, n, K)
  for (k in 1:K) { i <- A == k; corr[k] <- sum(w[i] * res[i]) / sum(w[i]); ICc[i, k] <- n * w[i] * (res[i] - corr[k]) / sum(w[i]) }
  mu <- colMeans(M) + corr
  t(sapply(2:K, function(k) { dA <- M[, k] - M[, 1]; c(est = mu[k] - mu[1], se = sqrt(sum(((dA - mean(dA)) + ICc[, k] - ICc[, 1])^2)) / n) }))
}

one_rep_b <- function(n, K, gamma, lambda) {
  d  <- gen_data(n, K, gamma, lambda)
  if (any(tabulate(d$A, K) < 5)) return(NULL)
  ps <- fit_gps(d, K, flag = lambda > 0)
  ea <- ps[cbind(seq_len(n), d$A)]
  X  <- as.matrix(d[, paste0("X", 1:5)])
  FM <- list(main = X, sq = cbind(X, X[, 1]^2, X[, 2]^2, X[, 3]^2))
  # ---- outcome-blind lambda selection ----
  sel <- sapply(RULES, function(r) {
    cg <- CG[CG <= r$cap]
    smd <- sapply(cg, function(cc) maxsmd(FM[[r$f]], d$A, wfun(ea, cc, K), K))
    ok <- which(smd <= r$t); if (length(ok)) cg[min(ok)] else cg[which.min(smd)]
  })
  Xl <- cbind(1, X); XD <- list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2))
  out <- list()
  for (y in c("hom", "het")) {
    Y <- d[[paste0("Y_", y)]]
    out[[length(out) + 1]] <- data.frame(method = "OW", om = "none", effect = y, arm = 2:K,
      est = estimate(Y, d$A, (1 / ea) / rowSums(1 / ps), K)[, "est"], se = NA, lam = NA)
    for (om in names(XD)) {
      M <- matrix(NA_real_, n, K)
      for (k in 1:K) { i <- d$A == k; b <- lm.fit(XD[[om]][i, , drop = FALSE], Y[i])$coefficients; b[is.na(b)] <- 0; M[, k] <- XD[[om]] %*% b }
      for (r in names(RULES)) {
        a <- aug(M, Y, d$A, wfun(ea, sel[[r]], K), K)
        out[[length(out) + 1]] <- data.frame(method = r, om = om, effect = y, arm = 2:K, est = a[, "est"], se = a[, "se"], lam = sel[[r]])
      }
    }
  }
  do.call(rbind, out)
}

RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
for (s in SCEN) {
  g <- grid[s, ]; f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s))
  if (file.exists(f)) next
  set.seed(20261115 + s)
  res <- mclapply(seq_len(N_REP), function(r) {
    o <- tryCatch(one_rep_b(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o
  }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res)
  res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL)
  saveRDS(res, f)
  cat(sprintf("Scenario %2d done %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
