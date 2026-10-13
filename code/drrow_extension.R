############################################################################
# Augmented / doubly robust ROW -- simulation extension
#
# ROW tilting:   h_l(X) = 1 / (sum_k 1/e_k(X) + l),  l = c * K^2
# Outcome models (per-arm OLS):
#   lin  = X1..X5            (correct for homogeneous effects; misspecified, no X1^2, for heterogeneous)
#   quad = X1..X5 + X1^2, X2^2, X3^2   (flexible; correct in all scenarios)
#
# A-ROW  (own estimand tau_l):
#   mu_k = sum_i h_i m_k(X_i) / sum_i h_i  +  Hajek_{A=k, w=h/e_k}( Y - m_k )
# DR-ROW (target: ATE):
#   mu_k = mean_i m_k(X_i)                 +  Hajek_{A=k, w=h/e_k}( Y - m_k )
#   l = Inf -> standard AIPW (Hajek form); smaller l -> more stable weights,
#   outcome model extrapolates into the poor-overlap tails
# SEs: influence-function based, nuisance models treated as fixed.
# Usage: Rscript drrow_extension.R [n_rep] [n_cores] [outdir]
############################################################################
Sys.setenv(SIM_SOURCE_ONLY = "1")
source("ow_vs_iptw_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_drrow"
dir.create(OUTDIR, showWarnings = FALSE)

C_GRID <- c(0, 0.5, 2, 8, 32, Inf)               # 0 = overlap weights, Inf = IPW

dr_estimates <- function(Y, A, X1, ps, K, c) {  # X1 = outcome design matrix
  n  <- length(Y); S <- rowSums(1 / ps)
  h  <- if (is.infinite(c)) rep(1, n) else 1 / (S + c * K^2)
  M  <- matrix(NA_real_, n, K)                    # m_k(X_i) for every patient and arm
  for (k in 1:K) {
    i <- A == k
    b <- lm.fit(X1[i, , drop = FALSE], Y[i])$coefficients
    b[is.na(b)] <- 0
    M[, k] <- X1 %*% b
  }
  res  <- Y - M[cbind(seq_len(n), A)]
  corr <- numeric(K); ICc <- matrix(0, n, K)      # weighted residual correction + its IF
  for (k in 1:K) {
    i <- A == k; w <- h[i] / ps[i, k]
    corr[k] <- sum(w * res[i]) / sum(w)
    ICc[i, k] <- n * w * (res[i] - corr[k]) / sum(w)
  }
  # DR-ROW (ATE): plug-in over everyone
  muA <- colMeans(M) + corr
  # A-ROW (own estimand): plug-in over the h-tilted population
  hn  <- h / sum(h)
  muR <- colSums(M * hn) + corr
  out <- list()
  for (k in 2:K) {
    dA <- M[, k] - M[, 1]
    ifA <- (dA - mean(dA)) + ICc[, k] - ICc[, 1]
    ifR <- n * hn * (dA - sum(hn * dA)) + ICc[, k] - ICc[, 1]
    out[[length(out) + 1]] <- data.frame(arm = k,
      est_dr = muA[k] - muA[1], se_dr = sqrt(sum(ifA^2)) / n,
      est_ar = muR[k] - muR[1], se_ar = sqrt(sum(ifR^2)) / n)
  }
  do.call(rbind, out)
}

one_rep_dr <- function(n, K, gamma, lambda) {
  d  <- gen_data(n, K, gamma, lambda)             # identical datasets to the main study
  if (any(tabulate(d$A, K) < 5)) return(NULL)
  ps <- fit_gps(d, K, flag = lambda > 0)
  Xl <- cbind(1, as.matrix(d[, paste0("X", 1:5)]))
  XD <- list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2))
  out <- list()
  for (om in names(XD)) for (c in C_GRID) for (y in c("hom", "het")) {
    r <- dr_estimates(d[[paste0("Y_", y)]], d$A, XD[[om]], ps, K, c)
    lab <- if (is.infinite(c)) "Inf" else format(c)
    out[[length(out) + 1]] <- rbind(
      data.frame(method = paste0("DR-ROW-c", lab), om = om, effect = y, arm = r$arm, est = r$est_dr, se = r$se_dr),
      data.frame(method = paste0("A-ROW-c", lab),  om = om, effect = y, arm = r$arm, est = r$est_ar, se = r$se_ar))
  }
  do.call(rbind, out)
}

RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
for (s in seq_len(nrow(grid))) {
  g <- grid[s, ]; f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s))
  if (file.exists(f)) next
  set.seed(20261115 + s)
  res <- mclapply(seq_len(N_REP), function(r) {
    o <- tryCatch(one_rep_dr(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o
  }, mc.cores = CORES)
  res <- do.call(rbind, res)
  res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL)
  saveRDS(res, f)
  cat(sprintf("Scenario %2d/%d done %.1f min\n", s, nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
