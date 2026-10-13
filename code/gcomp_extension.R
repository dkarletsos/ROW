############################################################################
# Outcome-regression (g-computation) comparator for DR-ROW
#   mu_k = mean_i m_k(X_i), per-arm OLS (lin / quad), SE from influence function
#   (nuisance treated as fixed: IF_i = m_k(X_i) - m_1(X_i) - tau + residual terms)
# Same datasets as the main study (same seeds, data generated first).
############################################################################
Sys.setenv(SIM_SOURCE_ONLY = "1")
source("ow_vs_iptw_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_gcomp"
dir.create(OUTDIR, showWarnings = FALSE)

gcomp <- function(Y, A, X, K) {
  n <- length(Y); M <- matrix(NA_real_, n, K); ICr <- matrix(0, n, K)
  for (k in 1:K) {
    i <- A == k; Xi <- X[i, , drop = FALSE]
    f <- lm.fit(Xi, Y[i]); b <- f$coefficients; b[is.na(b)] <- 0
    M[, k] <- X %*% b
    # residual contribution of OLS estimation to mean(m_k): xbar' (X'X)^-1 x_i r_i
    XtXi <- tryCatch(solve(crossprod(Xi)), error = function(e) MASS::ginv(crossprod(Xi)))
    ICr[i, k] <- n * as.vector(Xi %*% (XtXi %*% colMeans(X))) * f$residuals
  }
  out <- list()
  for (k in 2:K) {
    dA <- M[, k] - M[, 1]
    ifk <- (dA - mean(dA)) + ICr[, k] - ICr[, 1]
    out[[length(out) + 1]] <- data.frame(arm = k, est = mean(dA), se = sqrt(sum(ifk^2)) / n)
  }
  do.call(rbind, out)
}

one_rep_g <- function(n, K, gamma, lambda) {
  d <- gen_data(n, K, gamma, lambda)
  if (any(tabulate(d$A, K) < 5)) return(NULL)
  Xl <- cbind(1, as.matrix(d[, paste0("X", 1:5)]))
  XD <- list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2))
  out <- list()
  for (om in names(XD)) for (y in c("hom", "het")) {
    r <- gcomp(d[[paste0("Y_", y)]], d$A, XD[[om]], K)
    out[[length(out) + 1]] <- data.frame(method = "GCOMP", om = om, effect = y, r)
  }
  do.call(rbind, out)
}

RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
for (s in seq_len(nrow(grid))) {
  g <- grid[s, ]; f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s))
  if (file.exists(f)) next
  set.seed(20261115 + s)
  res <- mclapply(seq_len(N_REP), function(r) {
    o <- tryCatch(one_rep_g(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o
  }, mc.cores = CORES)
  res <- do.call(rbind, res)
  res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL)
  saveRDS(res, f)
  cat(sprintf("Scenario %2d/%d done %.1f min\n", s, nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
