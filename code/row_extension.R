############################################################################
# ROW: Regularised Overlap Weights -- simulation extension
#
#   h_lambda(X) = 1 / ( sum_k 1/e_k(X) + lambda ),   w = h_lambda(X) / e_A(X)
#   lambda = 0 -> overlap weights (ATO);  lambda -> Inf -> IPTW (ATE)
#   lambda is reported on the scale c = lambda / K^2 (S(X) = K^2 at equipoise)
#
# Reuses the data-generating mechanism and estimators of ow_vs_iptw_sim.R
# Usage: Rscript row_extension.R [n_rep] [n_cores] [outdir]
############################################################################
Sys.setenv(SIM_SOURCE_ONLY = "1")
source("ow_vs_iptw_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_row"
dir.create(OUTDIR, showWarnings = FALSE)

C_GRID <- c(0.5, 2, 8, 32)                       # lambda = c * K^2

row_weights <- function(ps, ea, K, c) {
  h <- 1 / (rowSums(1 / ps) + c * K^2)
  h / ea
}

# extend the weight set (one_rep() looks make_weights up in the global env)
make_weights_base <- make_weights
make_weights <- function(d, ps, K) {
  W  <- make_weights_base(d, ps, K)
  ea <- ps[cbind(seq_len(nrow(d)), d$A)]
  for (c in C_GRID) W[[sprintf("ROW-c%g", c)]] <- row_weights(ps, ea, K, c)
  W
}

true_row <- function(K, gamma, lambda, N = 2e6) {
  X <- gen_X(N); e <- true_ps(X, K, gamma, lambda); S <- rowSums(1 / e)
  out <- list()
  for (c in C_GRID) {
    h <- 1 / (S + c * K^2)
    for (k in 2:K) {
      tau <- DELTA[k] + ETA[k] * (X[, 1]^2 - 1)
      out[[length(out) + 1]] <- data.frame(K = K, arm = k, method = sprintf("ROW-c%g", c),
        T_hom = DELTA[k], T_het = sum(h * tau) / sum(h))
    }
  }
  do.call(rbind, out)
}

set.seed(20261115)
truths <- do.call(rbind, lapply(seq_len(nrow(grid)), function(s) {
  g <- grid[s, ]; if (g$n != 500) return(NULL)
  cbind(g[rep(1, 1), c("overlap", "positivity")], true_row(g$K, g$gamma, g$lambda), row.names = NULL)
}))
write.csv(truths, file.path(OUTDIR, "true_row_estimands.csv"), row.names = FALSE)

RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
for (s in seq_len(nrow(grid))) {
  g <- grid[s, ]; f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s))
  if (file.exists(f)) next
  set.seed(20261115 + s)                         # same seeds as the main study -> same datasets
  res <- mclapply(seq_len(N_REP), function(r) {
    o <- tryCatch(one_rep(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o
  }, mc.cores = CORES)
  res <- do.call(rbind, res)
  res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL)
  saveRDS(res, f)
  cat(sprintf("Scenario %2d/%d done %.1f min\n", s, nrow(grid), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
