############################################################################
# Held-out validation of ROW (balance rule, cap fixed IN ADVANCE at 8K^2).
# A NEW data-generating design, not used for any tuning decision:
#   * 8 confounders (4 continuous incl. one strongly skewed, 4 binary), correlated
#   * new treatment model coefficients; positivity violation via X2 > 1.2 (reference arm)
#   * heterogeneity through a skewed modifier and an interaction (X3 + X2*X5)
#   * three outcome types: continuous homogeneous, continuous heterogeneous, BINARY (risk difference)
#   * n = 1,000 (a sample size not used in the main study)
# Same 11 estimators as the main study, same implementations.
# Usage: Rscript heldout_sim.R [n_rep] [cores] [outdir]
############################################################################
suppressPackageStartupMessages({ library(nnet); library(parallel) })
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 500L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_heldout"
dir.create(OUTDIR, showWarnings = FALSE)

grid <- expand.grid(K = c(3, 4), gamma = c(1.0, 1.6), pv = c(0, 1))
grid$scenario <- seq_len(nrow(grid)); grid$n <- 1000
grid$overlap <- ifelse(grid$gamma == 1, "moderate", "poor"); grid$positivity <- ifelse(grid$pv == 1, "violated", "respected")

B2 <- rbind(c( 0.6, -0.7,  0.5,  0.3, -0.4,  0.5,  0.3, -0.2),
            c(-0.5,  0.8, -0.6,  0.4,  0.5, -0.3, -0.4,  0.4),
            c( 0.7,  0.5,  0.6, -0.5,  0.3,  0.4, -0.3,  0.5))
AL2 <- c(-0.4, -0.9, -0.7)
D_C <- c(0, -1.5, -2.5, -1.0); E_C <- c(0, 1.5, -1.2, 1.0)    # continuous outcome
D_B <- c(0, -0.8, -0.4, 0.5);  E_B <- c(0, 0.7, -0.6, 0.5)    # binary outcome (logit scale)

gen_X <- function(n) {
  S <- matrix(0.3, 4, 4); diag(S) <- 1; Z <- matrix(rnorm(n * 4), n) %*% chol(S)
  X3 <- (exp(Z[, 3]) - exp(0.5)) / sqrt((exp(1) - 1) * exp(1))        # skewed, standardised
  X5 <- rbinom(n, 1, 0.4); X6 <- rbinom(n, 1, plogis(Z[, 1])); X7 <- rbinom(n, 1, 0.3); X8 <- rbinom(n, 1, plogis(-0.5 + 0.6 * Z[, 2]))
  cbind(X1 = Z[, 1], X2 = Z[, 2], X3 = X3, X4 = Z[, 4], X5 = X5, X6 = X6, X7 = X7, X8 = X8)
}
true_ps <- function(X, K, gamma, pv) {
  Xc <- sweep(X, 2, colMeans(X))
  lp <- cbind(0, sweep(gamma * Xc %*% t(B2[seq_len(K - 1), , drop = FALSE]), 2, AL2[seq_len(K - 1)], "+"))
  lp[, 1] <- lp[, 1] - 5 * pv * (X[, 2] > 1.2)
  e <- exp(lp - apply(lp, 1, max)); e / rowSums(e)
}
mod <- function(X) X[, 3] + X[, 2] * X[, 5] - 0.4
mu0c <- function(X) 5 + 0.8 * X[, 1] + 0.6 * X[, 2] - 0.7 * X[, 3] + 0.4 * X[, 4] + 0.5 * X[, 5] - 0.4 * X[, 6] + 0.3 * X[, 7] + 0.6 * X[, 8]
lpb  <- function(X) -0.5 + 0.5 * X[, 1] + 0.4 * X[, 2] - 0.5 * X[, 3] + 0.3 * X[, 5] + 0.4 * X[, 8]
gen_data <- function(n, K, gamma, pv) {
  X <- gen_X(n); e <- true_ps(X, K, gamma, pv)
  A <- apply(e, 1, function(p) sample.int(K, 1, prob = p)); eps <- rnorm(n, 0, 2)
  list(X = X, A = A,
       Y = list(cont_hom = mu0c(X) + D_C[A] + eps,
                cont_het = mu0c(X) + D_C[A] + E_C[A] * mod(X) + eps,
                bin_het  = rbinom(n, 1, plogis(lpb(X) + D_B[A] + E_B[A] * X[, 2]))))
}
truth <- function(K, gamma, pv, N = 1e6) {
  X <- gen_X(N); e <- true_ps(X, K, gamma, pv); m <- mod(X); p1 <- plogis(lpb(X))
  data.frame(arm = 2:K, cont_hom = D_C[2:K],
             cont_het = sapply(2:K, function(k) D_C[k] + E_C[k] * mean(m)),
             bin_het  = sapply(2:K, function(k) mean(plogis(lpb(X) + D_B[k] + E_B[k] * X[, 2]) - p1)),
             share_minps_lt05 = mean(apply(e, 1, min) < 0.05))
}

# ---- estimator building blocks (identical to the main study) ----
fit_gps <- function(X, A, K, pv) {
  D <- data.frame(X, A = factor(A, levels = 1:K)); if (pv) D$F <- as.numeric(D$X2 > 1.2)
  f <- suppressWarnings(multinom(A ~ ., data = D, trace = FALSE, maxit = 500)); fitted(f)
}
trunc_w <- function(w, lo, hi) { q <- quantile(w, c(lo, hi)); pmin(pmax(w, q[1]), q[2]) }
hajek <- function(Y, A, w, K) { mu <- sapply(1:K, function(k) { i <- A == k; sum(w[i] * Y[i]) / sum(w[i]) }); mu[-1] - mu[1] }
fit_M <- function(Y, A, Xd, K) { M <- matrix(NA_real_, length(Y), K)
  for (k in 1:K) { i <- A == k; b <- lm.fit(Xd[i, , drop = FALSE], Y[i])$coefficients; b[is.na(b)] <- 0; M[, k] <- Xd %*% b }; M }
aug_if <- function(M, Y, A, w, K) {
  n <- length(Y); res <- Y - M[cbind(seq_len(n), A)]; corr <- numeric(K); ICc <- matrix(0, n, K)
  for (k in 1:K) { i <- A == k; corr[k] <- sum(w[i] * res[i]) / sum(w[i]); ICc[i, k] <- n * w[i] * (res[i] - corr[k]) / sum(w[i]) }
  mu <- colMeans(M) + corr
  list(est = mu[-1] - mu[1], IF = sapply(2:K, function(k) { dA <- M[, k] - M[, 1]; (dA - mean(dA)) + ICc[, k] - ICc[, 1] }))
}
eb_weights <- function(F, tbar, maxit = 100) {
  Fc <- sweep(F, 2, tbar); th <- rep(0, ncol(F))
  obj <- function(t) { z <- Fc %*% t; log(sum(exp(z - max(z)))) + max(z) }
  for (it in 1:maxit) {
    z <- as.vector(Fc %*% th); p <- exp(z - max(z)); p <- p / sum(p)
    g <- colSums(Fc * p); if (max(abs(g)) < 1e-8) break
    H <- crossprod(Fc * sqrt(p)) - tcrossprod(g)
    step <- tryCatch(solve(H + diag(1e-8, ncol(F)), g), error = function(e) g)
    s <- 1; f0 <- obj(th); while (obj(th - s * step) > f0 && s > 1e-6) s <- s / 2
    th <- th - s * step
  }
  z <- as.vector(Fc %*% th); p <- exp(z - max(z)); p / sum(p)
}
maxsmd <- function(Fm, A, w, K) { mu <- colMeans(Fm); sdv <- apply(Fm, 2, sd)
  max(sapply(1:K, function(k) { i <- A == k; max(abs(colSums(Fm[i, , drop = FALSE] * w[i]) / sum(w[i]) - mu) / sdv) })) }

one_rep <- function(n, K, gamma, pv) {
  d <- gen_data(n, K, gamma, pv); A <- d$A; X <- d$X
  if (any(tabulate(A, K) < 10)) return(NULL)
  ps <- fit_gps(X, A, K, pv); ea <- ps[cbind(seq_len(n), A)]; h <- 1 / rowSums(1 / ps)
  # ROW: balance rule, cap fixed at 8K^2 (features: covariates + squares of non-binary covariates)
  Fm <- cbind(X, X[, 1:4]^2); cg <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)
  smd <- sapply(cg, function(cc) maxsmd(Fm, A, 1 / (ea + 1 / (cc * K^2)), K))
  ok <- which(smd <= 0.1); csel <- if (length(ok)) cg[min(ok)] else cg[which.min(smd)]
  wrow <- 1 / (ea + 1 / (csel * K^2))
  webm <- numeric(n); tb <- colMeans(X); for (k in 1:K) { i <- A == k; webm[i] <- eb_weights(X[i, , drop = FALSE], tb) }
  marg <- tabulate(A, K)[A] / n; wipw <- 1 / ea
  XD <- list(lin = cbind(1, X), flex = cbind(1, X, X[, 1:4]^2, X[, 2] * X[, 5]))
  out <- list()
  for (y in names(d$Y)) {
    Y <- d$Y[[y]]
    add <- function(m, om, est) out[[length(out) + 1]] <<- data.frame(method = m, om = om, effect = y, arm = 2:K, est = est, lam = csel)
    add("OW", "none", hajek(Y, A, h / ea, K)); add("EB-main", "none", hajek(Y, A, webm, K))
    for (om in names(XD)) {
      M <- fit_M(Y, A, XD[[om]], K); res <- Y - M[cbind(seq_len(n), A)]
      add("GCOMP", om, colMeans(M)[-1] - colMeans(M)[1])
      add("AIPW", om, aug_if(M, Y, A, wipw, K)$est)
      add("AIPW-T1/99", om, aug_if(M, Y, A, trunc_w(marg / ea, 0.01, 0.99), K)$est)
      add("AIPW-T5/95", om, aug_if(M, Y, A, trunc_w(marg / ea, 0.05, 0.95), K)$est)
      add("ROW", om, aug_if(M, Y, A, wrow, K)$est)
      add("AEB-main", om, aug_if(M, Y, A, webm, K)$est)
      wow <- h / ea; corr <- sapply(1:K, function(k) { i <- A == k; sum(wow[i] * res[i]) / sum(wow[i]) })
      muO <- colSums(M * (h / sum(h))) + corr; add("AOW", om, muO[-1] - muO[1])
      q <- quantile(wipw^2, c(0.05, 0.95)); LG <- exp(seq(log(q[1]), log(q[2]), length.out = 30))
      fits <- lapply(LG, function(l) aug_if(M, Y, A, wipw * l / (wipw^2 + l), K)); Ai <- aug_if(M, Y, A, wipw, K)
      su <- sapply(seq_len(K - 1), function(j) {
        t <- sapply(fits, function(f) f$est[j]); v <- sapply(fits, function(f) sum(as.matrix(f$IF)[, j]^2) / n^2)
        vd <- sapply(fits, function(f) sum((as.matrix(f$IF)[, j] - as.matrix(Ai$IF)[, j])^2) / n^2)
        t[which.min((abs(t - Ai$est[j]) + 2 * sqrt(vd))^2 + v)] })
      add("SU", om, su)
      pb <- pmax(ps, 0.01); Qs <- M
      for (k in 1:K) { i <- A == k; Hk <- 1 / pb[i, k]; eps <- sum(Hk * res[i]) / sum(Hk^2); Qs[, k] <- M[, k] + eps / pb[, k] }
      add("TMLE", om, colMeans(Qs)[-1] - colMeans(Qs)[1])
    }
  }
  do.call(rbind, out)
}

if (Sys.getenv("HELDOUT_SOURCE_ONLY") != "1") {
  RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
  set.seed(777)
  if (!file.exists(file.path(OUTDIR, "truth.csv"))) {
    TR <- do.call(rbind, lapply(seq_len(nrow(grid)), function(s) { g <- grid[s, ]; cbind(scenario = s, truth(g$K, g$gamma, g$pv)) }))
    write.csv(TR, file.path(OUTDIR, "truth.csv"), row.names = FALSE)
    cat(sprintf("Truth done %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  for (s in seq_len(nrow(grid))) {
    g <- grid[s, ]; f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s)); if (file.exists(f)) next
    set.seed(40000 + s)
    res <- mclapply(seq_len(N_REP), function(r) { o <- tryCatch(one_rep(g$n, g$K, g$gamma, g$pv), error = function(e) NULL)
      if (!is.null(o)) o$rep <- r; o }, mc.cores = CORES, mc.set.seed = TRUE)
    res <- do.call(rbind, res); res$scenario <- s; saveRDS(res, f)
    cat(sprintf("Scenario %d done %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
