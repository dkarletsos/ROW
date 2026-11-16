############################################################################
# Stress-test design for ROW (cap fixed in advance at 8K^2). Every assumption under pressure:
#  * outcome models SEVERELY misspecified: sin(), step function, X2*X4 and X4*X6 interactions
#    (neither the linear nor the flexible working model contains them)
#  * propensity model misspecified (true treatment model has X1^2 and X2*X5; fitted: main effects) - factor
#  * overlap: poor (gamma 1.6) vs very poor (gamma 2.2) + structural positivity violation
#  * heavy-tailed errors (t3), binary outcome; n = 500 and 2,000; K = 3, 4
# Reuses the estimator implementations of heldout_sim.R unchanged.
############################################################################
Sys.setenv(HELDOUT_SOURCE_ONLY = "1"); source("heldout_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 300L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_stress"
dir.create(OUTDIR, showWarnings = FALSE)

sgrid <- expand.grid(K = c(3, 4), level = c("poor", "verypoor_pv"), psmis = c(0, 1), n = c(500, 2000), stringsAsFactors = FALSE)
sgrid$gamma <- ifelse(sgrid$level == "poor", 1.6, 2.2); sgrid$pv <- as.numeric(sgrid$level != "poor"); sgrid$scenario <- seq_len(nrow(sgrid))

true_ps_s <- function(X, K, gamma, pv, psmis) {
  Xc <- sweep(X, 2, colMeans(X))
  lp <- cbind(0, sweep(gamma * Xc %*% t(B2[seq_len(K - 1), , drop = FALSE]), 2, AL2[seq_len(K - 1)], "+"))
  if (psmis) { nl <- cbind(X[, 1]^2 - 1, X[, 2] * X[, 5]); G <- rbind(c(0.9, -0.8), c(-0.7, 0.9), c(0.8, 0.7))
    lp[, -1] <- lp[, -1] + gamma * nl %*% t(G[seq_len(K - 1), , drop = FALSE]) }
  lp[, 1] <- lp[, 1] - 5 * pv * (X[, 2] > 1.2)
  e <- exp(lp - apply(lp, 1, max)); e / rowSums(e)
}
ms   <- function(X) 1.0 * (X[, 1] > 0.5) + X[, 2] * X[, 4] + 0.5 * X[, 3] - 0.31
mu0s <- function(X) 5 + 1.2 * sin(1.5 * X[, 1]) + 0.8 * X[, 2]^2 - 0.6 * X[, 3] + 0.5 * X[, 4] * X[, 6] + 0.5 * X[, 5] + 0.6 * X[, 8] + 1.0 * (X[, 1] > 0.5)
lpbs <- function(X) -0.8 + 0.8 * sin(1.5 * X[, 1]) + 0.5 * X[, 2]^2 - 0.5 * X[, 3] + 0.4 * X[, 8]
PSMIS <- 0
gen_data <- function(n, K, gamma, pv) {
  X <- gen_X(n); e <- true_ps_s(X, K, gamma, pv, PSMIS)
  A <- apply(e, 1, function(p) sample.int(K, 1, prob = p)); eps <- rnorm(n, 0, 2); m <- ms(X)
  list(X = X, A = A,
       Y = list(cont_hom    = mu0s(X) + D_C[A] + eps,
                cont_het    = mu0s(X) + D_C[A] + E_C[A] * m + eps,
                cont_het_t3 = mu0s(X) + D_C[A] + E_C[A] * m + 2 * rt(n, 3) / sqrt(3),
                bin_het     = rbinom(n, 1, plogis(lpbs(X) + D_B[A] + E_B[A] * m))))
}
truth_s <- function(K, gamma, pv, psmis, N = 1e6) {
  X <- gen_X(N); e <- true_ps_s(X, K, gamma, pv, psmis); h <- 1 / rowSums(1 / e); m <- ms(X); p1 <- plogis(lpbs(X) + E_B[1] * m)
  data.frame(arm = 2:K, cont_hom = D_C[2:K],
             cont_het = sapply(2:K, function(k) D_C[k] + E_C[k] * mean(m)),
             cont_het_t3 = sapply(2:K, function(k) D_C[k] + E_C[k] * mean(m)),
             bin_het = sapply(2:K, function(k) mean(plogis(lpbs(X) + D_B[k] + E_B[k] * m) - p1)),
             ATO_cont_het = sapply(2:K, function(k) D_C[k] + E_C[k] * sum(h * m) / sum(h)),
             share_minps_lt05 = mean(apply(e, 1, min) < 0.05), share_minps_lt01 = mean(apply(e, 1, min) < 0.01))
}

RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time(); set.seed(888)
if (!file.exists(file.path(OUTDIR, "truth.csv"))) {
  TR <- do.call(rbind, lapply(seq_len(nrow(sgrid)), function(s) { g <- sgrid[s, ]; cbind(scenario = s, truth_s(g$K, g$gamma, g$pv, g$psmis)) }))
  write.csv(TR, file.path(OUTDIR, "truth.csv"), row.names = FALSE); write.csv(sgrid, file.path(OUTDIR, "grid.csv"), row.names = FALSE)
  cat(sprintf("Truth done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
for (s in seq_len(nrow(sgrid))) {
  g <- sgrid[s, ]; f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s)); if (file.exists(f)) next
  PSMIS <- g$psmis; set.seed(50000 + s)
  res <- mclapply(seq_len(N_REP), function(r) { o <- tryCatch(one_rep(g$n, g$K, g$gamma, g$pv), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r; o }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res); res$scenario <- s; saveRDS(res, f)
  cat(sprintf("Scenario %d done %.1f min (%d reps)\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins")), length(unique(res$rep))))
}
