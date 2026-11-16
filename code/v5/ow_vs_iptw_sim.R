############################################################################
# Overlap Weighting vs Trimmed IPTW for Multi-Group Comparisons
# Simulation study -- PHUSE EU Connect 2026
#
# Extends the PHUSE EU Connect 2025 multigroup_iptw code (stabilised IPTW
# with a generalised propensity score (GPS) and weight trimming) with a
# head-to-head comparison against multinomial overlap weighting (Li & Li 2019).
#
# Dependencies: base R + nnet only (sandwich variance written out by hand)
# Usage:   Rscript ow_vs_iptw_sim.R [n_rep] [n_cores]
############################################################################

suppressPackageStartupMessages({ library(nnet); library(parallel) })

args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else max(1L, detectCores() - 1L)
OUTDIR <- if (length(args) >= 3) args[3] else "results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------
# 1. Data-generating mechanism
# --------------------------------------------------------------------------
# Five baseline confounders (3 continuous, 2 binary).
# Treatment: multinomial logit; gamma scales all slopes (overlap strength).
# Positivity violation: in the region X1 > 1 (~16% of patients) the
# reference arm (arm 1) is heavily disfavoured (logit penalty lambda).
# Outcome: continuous, with arm-specific effects that may vary with X1.

B <- rbind(                         # PS slopes for arms 2..4 vs arm 1
  c( 0.8, -0.5,  0.3,  0.4, -0.3),
  c(-0.6,  0.7, -0.4, -0.3,  0.5),
  c( 1.0,  0.6, -0.5,  0.2,  0.4))
ALPHA <- c(-0.6, -1.0, -0.8)       # intercepts: reference arm is the largest
DELTA <- c(0, -2.0, -3.0, -1.5)     # average effect vs arm 1 at X1 = 0
ETA   <- c(0, -1.2,  1.0, -0.8)     # effect modification by severity (X1^2 - 1)

gen_X <- function(n) {
  X1 <- rnorm(n); X2 <- 0.3 * X1 + sqrt(1 - 0.09) * rnorm(n); X3 <- rnorm(n)
  X4 <- rbinom(n, 1, 0.5); X5 <- rbinom(n, 1, plogis(-0.8 + 0.5 * X3))
  cbind(X1, X2, X3, X4, X5)
}

true_ps <- function(X, K, gamma, lambda) {
  Xc <- sweep(X, 2, c(0, 0, 0, 0.5, 0.33))      # centre binaries
  lp <- cbind(0, sweep(gamma * Xc %*% t(B[seq_len(K - 1), , drop = FALSE]), 2,
                       ALPHA[seq_len(K - 1)], "+"))
  lp[, 1] <- lp[, 1] - lambda * (X[, 1] > 1)     # positivity violation
  e <- exp(lp - apply(lp, 1, max)); e / rowSums(e)
}

gen_data <- function(n, K, gamma, lambda) {
  X <- gen_X(n); e <- true_ps(X, K, gamma, lambda)
  A <- apply(e, 1, function(p) sample.int(K, 1, prob = p))
  mu0 <- 10 + 1.0 * X[, 1] + 0.8 * X[, 2] - 0.5 * X[, 3] + 0.7 * X[, 4] + 0.5 * X[, 5]
  eps <- rnorm(n, 0, 2)
  Y_hom <- mu0 + DELTA[A] + eps                        # homogeneous effects
  Y_het <- mu0 + DELTA[A] + ETA[A] * (X[, 1]^2 - 1) + eps      # heterogeneous effects
  data.frame(X, A = A, Y_hom = Y_hom, Y_het = Y_het)
}

# True estimands (pairwise contrasts arm k vs arm 1) from a 2e6 super-population
true_estimands <- function(K, gamma, lambda, N = 2e6) {
  X <- gen_X(N); e <- true_ps(X, K, gamma, lambda)
  h <- 1 / rowSums(1 / e)                        # overlap tilting function
  tau <- sapply(2:K, function(k) (DELTA[k] + ETA[k] * (X[, 1]^2 - 1)) - 0)
  data.frame(arm = 2:K,
             ATE_hom = DELTA[2:K], ATO_hom = DELTA[2:K],
             ATE_het = colMeans(tau), ATO_het = colSums(tau * h) / sum(h),
             min_ps_p01 = unname(quantile(apply(e, 1, min), 0.01)),
             share_minps_lt05 = mean(apply(e, 1, min) < 0.05))
}

# --------------------------------------------------------------------------
# 2. Weighting methods
# --------------------------------------------------------------------------
fit_gps <- function(d, K, flag = FALSE) {
  # flag = TRUE: the structural positivity violation arises from a known
  # clinical rule (X1 > 1), so the analyst includes that indicator in the GPS model
  fml <- if (flag) factor(A, levels = 1:K) ~ X1 + X2 + X3 + X4 + X5 + I(X1 > 1) else
                   factor(A, levels = 1:K) ~ X1 + X2 + X3 + X4 + X5
  f <- suppressWarnings(multinom(fml, data = d,
                trace = FALSE, maxit = 500))
  p <- fitted(f); if (K == 2) p <- cbind(1 - p, p); p
}

trunc_w <- function(w, lo, hi) { q <- quantile(w, c(lo, hi)); pmin(pmax(w, q[1]), q[2]) }

make_weights <- function(d, ps, K) {
  n  <- nrow(d); idx <- cbind(seq_len(n), d$A)
  ea <- ps[idx]                                  # GPS of received arm
  marg <- tabulate(d$A, K)[d$A] / n
  sw <- marg / ea                                # stabilised IPTW (2025 approach)
  ow <- (1 / ea) / rowSums(1 / ps)               # multinomial overlap weights
  list(IPTW = sw,
       `IPTW-T1/99` = trunc_w(sw, 0.01, 0.99),   # 2025 primary analysis
       `IPTW-T5/95` = trunc_w(sw, 0.05, 0.95),
       OW = ow)
}

# Hajek arm means + robust (HC0-type) variance treating weights as known
estimate <- function(Y, A, w, K) {
  mu <- v <- numeric(K)
  for (k in 1:K) {
    i <- A == k; sw <- sum(w[i])
    mu[k] <- sum(w[i] * Y[i]) / sw
    v[k]  <- sum(w[i]^2 * (Y[i] - mu[k])^2) / sw^2
  }
  cbind(est = mu[-1] - mu[1], se = sqrt(v[-1] + v[1]))
}

ess_pct <- function(w, A, K) {
  100 * sum(sapply(1:K, function(k) { x <- w[A == k]; sum(x)^2 / sum(x^2) })) / length(w)
}

max_smd <- function(X, A, w, K) {
  # max |weighted SMD| over covariates and all pairwise arm contrasts;
  # denominator: unweighted pooled SD across arms
  m <- sapply(1:K, function(k) colSums(X[A == k, , drop = FALSE] * w[A == k]) / sum(w[A == k]))
  s <- sqrt(rowMeans(sapply(1:K, function(k) apply(X[A == k, , drop = FALSE], 2, var))))
  mx <- 0
  for (a in 1:(K - 1)) for (b in (a + 1):K) mx <- max(mx, abs(m[, a] - m[, b]) / s)
  mx
}

# --------------------------------------------------------------------------
# 3. One replicate
# --------------------------------------------------------------------------
one_rep <- function(n, K, gamma, lambda) {
  d  <- gen_data(n, K, gamma, lambda)
  if (any(tabulate(d$A, K) < 5)) return(NULL)
  ps <- fit_gps(d, K, flag = lambda > 0)
  W  <- make_weights(d, ps, K)
  X  <- as.matrix(d[, paste0("X", 1:5)])
  out <- list()
  for (m in names(W)) {
    w <- W[[m]]
    for (y in c("hom", "het")) {
      r <- estimate(d[[paste0("Y_", y)]], d$A, w, K)
      out[[length(out) + 1]] <- data.frame(method = m, effect = y, arm = 2:K,
        est = r[, "est"], se = r[, "se"], ess = ess_pct(w, d$A, K),
        smd = max_smd(X, d$A, w, K), wmax_ratio = max(w) / mean(w))
    }
  }
  do.call(rbind, out)
}

# --------------------------------------------------------------------------
# 4. Scenario grid and driver
# --------------------------------------------------------------------------
grid <- expand.grid(n = c(500, 2000), K = c(3, 4),
                    overlap = c("good", "moderate", "poor"),
                    positivity = c("none", "severe"), stringsAsFactors = FALSE)
grid$gamma  <- c(good = 0.5, moderate = 1.0, poor = 1.5)[grid$overlap]
grid$lambda <- c(none = 0, severe = 4)[grid$positivity]
grid$scenario <- seq_len(nrow(grid))

if (!identical(Sys.getenv("SIM_SOURCE_ONLY"), "1")) {
  set.seed(20261115)
  truths <- do.call(rbind, lapply(seq_len(nrow(grid)), function(s) {
    g <- grid[s, ]
    if (g$n != 500) return(NULL)                 # truth does not depend on n
    cbind(g[, c("K", "overlap", "positivity")], true_estimands(g$K, g$gamma, g$lambda), row.names = NULL)
  }))
  write.csv(truths, file.path(OUTDIR, "true_estimands.csv"), row.names = FALSE)

  RNGkind("L'Ecuyer-CMRG")
  t0 <- Sys.time()
  for (s in seq_len(nrow(grid))) {
    g <- grid[s, ]
    f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s))
    if (file.exists(f)) next
    set.seed(20261115 + s)                       # independent, reproducible stream per scenario
    res <- mclapply(seq_len(N_REP), function(r) {
      o <- tryCatch(one_rep(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
      if (!is.null(o)) o$rep <- r
      o
    }, mc.cores = CORES)
    res <- do.call(rbind, res)
    res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res)
    saveRDS(res, f)
    cat(sprintf("Scenario %2d/%d done (%s) %.1f min\n", s, nrow(grid),
                paste(g[1:4], collapse = "/"), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
