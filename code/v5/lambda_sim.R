############################################################################
# Paired comparison on identical datasets: OW vs ROW (augmented, lambda = 8K^2)
# vs outcome regression vs AIPW vs AIPW with 1/99-truncated weights.
#
# Every estimator is applied to the same simulated dataset in each replicate,
# so performance differences are paired (common random numbers).
# Target: ATE (pairwise contrasts arm k vs arm 1). OW is also scored vs ATO.
# Outcome models (per-arm OLS): lin = X1..X5 ; quad = lin + X1^2, X2^2, X3^2
#   (lin is misspecified for heterogeneous effects; quad is correct throughout)
# Usage: Rscript paired_sim.R [n_rep] [n_cores] [outdir]
############################################################################
Sys.setenv(SIM_SOURCE_ONLY = "1")
source("ow_vs_iptw_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 500L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_lambda"
SCEN   <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else seq_len(nrow(grid))
dir.create(OUTDIR, showWarnings = FALSE)

# Augmented estimator of the ATE with arbitrary weights for the residual correction
#   mu_k = mean_i m_k(X_i) + Hajek_{A=k, w}(Y - m_k)
aug_ate <- function(M, Y, A, w, K) {
  n <- length(Y); res <- Y - M[cbind(seq_len(n), A)]
  corr <- numeric(K); ICc <- matrix(0, n, K)
  for (k in 1:K) {
    i <- A == k; wi <- w[i]
    corr[k] <- sum(wi * res[i]) / sum(wi)
    ICc[i, k] <- n * wi * (res[i] - corr[k]) / sum(wi)
  }
  mu <- colMeans(M) + corr
  t(sapply(2:K, function(k) {
    dA <- M[, k] - M[, 1]
    ifk <- (dA - mean(dA)) + ICc[, k] - ICc[, 1]
    c(est = mu[k] - mu[1], se = sqrt(sum(ifk^2)) / n)
  }))
}

# Outcome regression (g-computation) with OLS-estimation term in the IF
gcomp_fit <- function(Y, A, X, K) {
  n <- length(Y); M <- matrix(NA_real_, n, K); ICr <- matrix(0, n, K)
  for (k in 1:K) {
    i <- A == k; Xi <- X[i, , drop = FALSE]
    f <- lm.fit(Xi, Y[i]); b <- f$coefficients; b[is.na(b)] <- 0
    M[, k] <- X %*% b
    XtXi <- tryCatch(solve(crossprod(Xi)), error = function(e) MASS::ginv(crossprod(Xi)))
    ICr[i, k] <- n * as.vector(Xi %*% (XtXi %*% colMeans(X))) * f$residuals
  }
  g <- t(sapply(2:K, function(k) {
    dA <- M[, k] - M[, 1]; ifk <- (dA - mean(dA)) + ICr[, k] - ICr[, 1]
    c(est = mean(dA), se = sqrt(sum(ifk^2)) / n)
  }))
  list(M = M, g = g)
}

one_rep_p <- function(n, K, gamma, lambda) {
  d  <- gen_data(n, K, gamma, lambda)
  if (any(tabulate(d$A, K) < 5)) return(NULL)
  ps <- fit_gps(d, K, flag = lambda > 0)
  idx <- cbind(seq_len(n), d$A); ea <- ps[idx]; S <- rowSums(1 / ps)
  marg <- tabulate(d$A, K)[d$A] / n
  W <- list(OW = (1 / ea) / S)
  for (cc in c(2, 32, 128, 512)) W[[paste0("ROW-c", cc)]] <- (1 / ea) / (S + cc * K^2)
  Xl <- cbind(1, as.matrix(d[, paste0("X", 1:5)]))
  XD <- list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2))
  out <- list()
  for (y in c("hom", "het")) {
    Y <- d[[paste0("Y_", y)]]
    r <- estimate(Y, d$A, W$OW, K)                  # OW: weighting only (no outcome model)
    out[[length(out) + 1]] <- data.frame(method = "OW", om = "none", effect = y, arm = 2:K,
                                         est = r[, "est"], se = r[, "se"])
    for (om in names(XD)) {
      gf <- gcomp_fit(Y, d$A, XD[[om]], K)
      for (m in setdiff(names(W), "OW")) {
        a <- aug_ate(gf$M, Y, d$A, W[[m]], K)
        out[[length(out) + 1]] <- data.frame(method = m, om = om, effect = y, arm = 2:K,
                                             est = a[, "est"], se = a[, "se"])
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
    o <- tryCatch(one_rep_p(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o
  }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res)
  res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL)
  saveRDS(res, f)
  cat(sprintf("Scenario %2d done %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
