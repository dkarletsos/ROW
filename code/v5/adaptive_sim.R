############################################################################
# Data-driven lambda for ROW (arm-specific regularised weights w = 1/(e_A + 1/lambda)).
# Family runs from outcome regression (lambda -> 0) to AIPW (lambda = Inf).
# Rules (per contrast, using influence-function variances; AIPW = low-bias anchor):
#   MSE : minimise  max(0, (t_l - t_Inf)^2 - V(t_l - t_Inf)) + V(t_l)
#   LEP : smallest lambda such that every larger lambda agrees with AIPW within 2 SD of the difference
# Same seeds as paired_sim.R -> identical datasets.
# Usage: Rscript adaptive_sim.R [n_rep] [n_cores] [outdir] [scenarios]
############################################################################
Sys.setenv(SIM_SOURCE_ONLY = "1")
source("ow_vs_iptw_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_adaptive"
SCEN   <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else seq_len(nrow(grid))
dir.create(OUTDIR, showWarnings = FALSE)
CG <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, Inf)   # lambda / K^2

fit_M <- function(Y, A, X, K) {
  M <- matrix(NA_real_, length(Y), K)
  for (k in 1:K) { i <- A == k; b <- lm.fit(X[i, , drop = FALSE], Y[i])$coefficients; b[is.na(b)] <- 0; M[, k] <- X %*% b }
  M
}
# returns est (K-1) and IF matrix (n x K-1) for weights w
aug_if <- function(M, Y, A, w, K) {
  n <- length(Y); res <- Y - M[cbind(seq_len(n), A)]
  corr <- numeric(K); ICc <- matrix(0, n, K)
  for (k in 1:K) { i <- A == k; wi <- w[i]; corr[k] <- sum(wi * res[i]) / sum(wi); ICc[i, k] <- n * wi * (res[i] - corr[k]) / sum(wi) }
  mu <- colMeans(M) + corr
  est <- mu[-1] - mu[1]
  IF <- sapply(2:K, function(k) { dA <- M[, k] - M[, 1]; (dA - mean(dA)) + ICc[, k] - ICc[, 1] })
  list(est = est, IF = matrix(IF, n))
}

one_rep_a <- function(n, K, gamma, lambda) {
  d  <- gen_data(n, K, gamma, lambda)
  if (any(tabulate(d$A, K) < 5)) return(NULL)
  ps <- fit_gps(d, K, flag = lambda > 0)
  ea <- ps[cbind(seq_len(n), d$A)]
  Xl <- cbind(1, as.matrix(d[, paste0("X", 1:5)]))
  XD <- list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2))
  out <- list()
  for (y in c("hom", "het")) {
    Y <- d[[paste0("Y_", y)]]
    out[[length(out) + 1]] <- data.frame(method = "OW", om = "none", effect = y, arm = 2:K,
      est = estimate(Y, d$A, (1 / ea) / rowSums(1 / ps), K)[, "est"], se = NA, lam = NA)
    for (om in names(XD)) {
      M <- fit_M(Y, d$A, XD[[om]], K)
      fits <- lapply(CG, function(cc) aug_if(M, Y, d$A, if (is.infinite(cc)) 1 / ea else 1 / (ea + 1 / (cc * K^2)), K))
      Ti <- fits[[length(CG)]]
      for (j in seq_len(K - 1)) {
        t  <- sapply(fits, function(f) f$est[j])
        v  <- sapply(fits, function(f) sum(f$IF[, j]^2) / n^2)
        vd <- sapply(fits, function(f) sum((f$IF[, j] - Ti$IF[, j])^2) / n^2)
        b2 <- pmax(0, (t - t[length(CG)])^2 - vd)
        jm <- which.min(b2 + v)
        ok <- abs(t - t[length(CG)]) <= 2 * sqrt(vd)
        ok[length(CG)] <- TRUE
        jl <- min(which(rev(cumprod(rev(ok))) == 1))      # smallest lambda with all larger lambdas agreeing
        out[[length(out) + 1]] <- data.frame(method = c("ROW-mse", "ROW-lep"), om = om, effect = y, arm = j + 1,
          est = t[c(jm, jl)], se = sqrt(v[c(jm, jl)]), lam = CG[c(jm, jl)])
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
    o <- tryCatch(one_rep_a(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o
  }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res)
  res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL)
  saveRDS(res, f)
  cat(sprintf("Scenario %2d done %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
