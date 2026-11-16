############################################################################
# Additional comparators on identical datasets (same seeds as paired_sim.R):
#   AOW      : augmented overlap weights (target ATO; Mao, Li & Greene 2019) - scored against the ATE
#   EB       : entropy balancing of each arm to the full-population moments (X, X1^2..X3^2), weighting only
#   AEB      : augmented entropy balancing (outcome model + EB-weighted residual correction)
#   SU       : AIPW with optimistic weight shrinkage w*lam/(w^2+lam) (Su et al. 2020), lam chosen per contrast by
#              minimising (|bias estimate| + 2 SE)^2 + variance ("DRs-upper"), 30 geometric values
#   TMLE     : TMLE with linear fluctuation, GPS bounded below at 0.01
############################################################################
Sys.setenv(SIM_SOURCE_ONLY = "1")
source("ow_vs_iptw_sim.R")
args   <- commandArgs(trailingOnly = TRUE)
N_REP  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
CORES  <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_comp"
SCEN   <- if (length(args) >= 4) as.integer(strsplit(args[4], ",")[[1]]) else seq_len(nrow(grid))
dir.create(OUTDIR, showWarnings = FALSE)

fit_M <- function(Y, A, X, K) {
  M <- matrix(NA_real_, length(Y), K)
  for (k in 1:K) { i <- A == k; b <- lm.fit(X[i, , drop = FALSE], Y[i])$coefficients; b[is.na(b)] <- 0; M[, k] <- X %*% b }
  M
}
# augmented estimator of the ATE with arm-specific weights w (Hajek residual correction); returns est, IF
aug_if <- function(M, Y, A, w, K) {
  n <- length(Y); res <- Y - M[cbind(seq_len(n), A)]; corr <- numeric(K); ICc <- matrix(0, n, K)
  for (k in 1:K) { i <- A == k; corr[k] <- sum(w[i] * res[i]) / sum(w[i]); ICc[i, k] <- n * w[i] * (res[i] - corr[k]) / sum(w[i]) }
  mu <- colMeans(M) + corr
  list(est = mu[-1] - mu[1], IF = sapply(2:K, function(k) { dA <- M[, k] - M[, 1]; (dA - mean(dA)) + ICc[, k] - ICc[, 1] }))
}
se_of <- function(IF) sqrt(colSums(as.matrix(IF)^2)) / nrow(as.matrix(IF))

# entropy balancing of arm k to target moments tbar (dual Newton with step halving)
eb_weights <- function(F, tbar, maxit = 100) {
  Fc <- sweep(F, 2, tbar); th <- rep(0, ncol(F))
  obj <- function(t) log(sum(exp(Fc %*% t - max(Fc %*% t)))) + max(Fc %*% t)
  for (it in 1:maxit) {
    e <- as.vector(exp(Fc %*% th - max(Fc %*% th))); p <- e / sum(e)
    g <- colSums(Fc * p); if (max(abs(g)) < 1e-8) break
    H <- crossprod(Fc * sqrt(p)) - tcrossprod(g)
    step <- tryCatch(solve(H + diag(1e-8, ncol(F)), g), error = function(e) g)
    s <- 1; f0 <- obj(th); while (obj(th - s * step) > f0 && s > 1e-6) s <- s / 2
    th <- th - s * step
  }
  e <- as.vector(exp(Fc %*% th - max(Fc %*% th))); e / sum(e)
}

one_rep_c <- function(n, K, gamma, lambda) {
  d  <- gen_data(n, K, gamma, lambda)
  if (any(tabulate(d$A, K) < 5)) return(NULL)
  ps <- fit_gps(d, K, flag = lambda > 0)
  ea <- ps[cbind(seq_len(n), d$A)]; S <- rowSums(1 / ps); h <- 1 / S
  X  <- as.matrix(d[, paste0("X", 1:5)]); Fm <- cbind(X, X[, 1]^2, X[, 2]^2, X[, 3]^2)
  tbar <- colMeans(Fm)
  web <- webm <- numeric(n); tbm <- colMeans(X)
  for (k in 1:K) { i <- d$A == k; web[i] <- eb_weights(Fm[i, , drop = FALSE], tbar); webm[i] <- eb_weights(X[i, , drop = FALSE], tbm) }
  Xl <- cbind(1, X); XD <- list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2))
  wipw <- 1 / ea
  out <- list()
  for (y in c("hom", "het")) {
    Y <- d[[paste0("Y_", y)]]
    r <- estimate(Y, d$A, web, K)
    out[[length(out) + 1]] <- data.frame(method = "EB", om = "none", effect = y, arm = 2:K, est = r[, "est"], se = r[, "se"])
    r2 <- estimate(Y, d$A, webm, K)
    out[[length(out) + 1]] <- data.frame(method = "EB-main", om = "none", effect = y, arm = 2:K, est = r2[, "est"], se = r2[, "se"])
    out[[length(out) + 1]] <- data.frame(method = "OW", om = "none", effect = y, arm = 2:K,
      est = estimate(Y, d$A, (1 / ea) * h, K)[, "est"], se = NA)
    for (om in names(XD)) {
      M <- fit_M(Y, d$A, XD[[om]], K)
      # augmented OW (ATO target): plug-in over the h-tilted population + OW-weighted residual correction
      res <- Y - M[cbind(seq_len(n), d$A)]; hn <- h / sum(h); wow <- h / ea
      corr <- sapply(1:K, function(k) { i <- d$A == k; sum(wow[i] * res[i]) / sum(wow[i]) })
      muO <- colSums(M * hn) + corr
      out[[length(out) + 1]] <- data.frame(method = "AOW", om = om, effect = y, arm = 2:K, est = muO[-1] - muO[1], se = NA)
      a <- aug_if(M, Y, d$A, web, K)
      out[[length(out) + 1]] <- data.frame(method = "AEB", om = om, effect = y, arm = 2:K, est = a$est, se = se_of(a$IF))
      a2 <- aug_if(M, Y, d$A, webm, K)
      out[[length(out) + 1]] <- data.frame(method = "AEB-main", om = om, effect = y, arm = 2:K, est = a2$est, se = se_of(a2$IF))
      # Su et al. optimistic shrinkage, DRs-upper selection per contrast
      q <- quantile(wipw^2, c(0.05, 0.95)); LG <- exp(seq(log(q[1]), log(q[2]), length.out = 30))
      fits <- lapply(LG, function(l) aug_if(M, Y, d$A, wipw * l / (wipw^2 + l), K))
      Ai <- aug_if(M, Y, d$A, wipw, K)
      for (j in seq_len(K - 1)) {
        t  <- sapply(fits, function(f) f$est[j]); v <- sapply(fits, function(f) sum(as.matrix(f$IF)[, j]^2) / n^2)
        vd <- sapply(fits, function(f) sum((as.matrix(f$IF)[, j] - as.matrix(Ai$IF)[, j])^2) / n^2)
        crit <- (abs(t - Ai$est[j]) + 2 * sqrt(vd))^2 + v; jm <- which.min(crit)
        out[[length(out) + 1]] <- data.frame(method = "SU", om = om, effect = y, arm = j + 1, est = t[jm], se = sqrt(v[jm]))
      }
      # TMLE, linear fluctuation, GPS bounded at 0.01
      pb <- pmax(ps, 0.01); Qs <- M; IFt <- matrix(0, n, K)
      for (k in 1:K) { i <- d$A == k; Hk <- i / pb[, k]
        eps <- sum(Hk[i] * res[i]) / sum(Hk[i]^2); Qs[, k] <- M[, k] + eps / pb[, k] }
      mu <- colMeans(Qs)
      for (k in 1:K) { i <- d$A == k; IFt[, k] <- (Qs[, k] - mu[k]) + i / pb[, k] * (Y - Qs[, k]) }
      out[[length(out) + 1]] <- data.frame(method = "TMLE", om = om, effect = y, arm = 2:K, est = mu[-1] - mu[1],
        se = sqrt(colSums((IFt[, -1, drop = FALSE] - IFt[, 1])^2)) / n)
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
    o <- tryCatch(one_rep_c(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r
    o
  }, mc.cores = CORES, mc.set.seed = TRUE)
  res <- do.call(rbind, res)
  res <- cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL)
  saveRDS(res, f)
  cat(sprintf("Scenario %2d done %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
