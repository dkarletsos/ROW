############################################################################
# PET and AIPW-PET (Orihara, Komukai & Li 2026), extended to K arms, on the SAME datasets as the
# three ROW designs (same seeds, cores and per-replicate RNG use, so replicate r is the identical dataset).
# Multi-arm tilting family: w_beta(x) = h_O(x)^beta, h_O = 1 / sum_j 1/e_j(x)  (binary case: {e(1-e)}^beta up to a constant).
# Per beta: Hajek (augmented) estimator of the w_beta-weighted effect; polynomial (degree q) OLS fit across a grid
# of K_B = 25 equally spaced beta in [beta1, 0.99]; extrapolation to beta = 0 (intercept).
# Hyperparameters by the authors' Algorithm 1, per contrast: target variance = kappa * var(IPW) (kappa = 5/6, PET)
# or var(AIPW) (kappa = 1, AIPW-PET); for each q in 1..4 the smallest beta1 in {0.44,...,0.69} meeting the target
# (else the largest); then the largest q meeting the target (else the q with the smallest variance).
# Usage: Rscript pet_sim.R <design> [n_rep] [cores]
############################################################################
args <- commandArgs(trailingOnly = TRUE)
DESIGN <- args[1]; NREP_P <- if (length(args) >= 2) as.integer(args[2]) else 500L
CORES_P <- if (length(args) >= 3) as.integer(args[3]) else 2L
suppressPackageStartupMessages({ library(nnet); library(parallel) })
B1C <- seq(0.44, 0.69, by = 0.05); QC <- 1:4; KB <- 25
GRIDS <- lapply(B1C, function(b1) seq(b1, 0.99, length.out = KB))
BALL <- sort(unique(c(0, unlist(GRIDS))))
PALPHA <- lapply(GRIDS, function(g) lapply(QC, function(q) { V <- outer(g, 0:q, `^`); solve(crossprod(V), t(V))[1, ] }))

fit_M <- function(Y, A, Xd, K) { M <- matrix(NA_real_, length(Y), K)
  for (k in 1:K) { i <- A == k; b <- lm.fit(Xd[i, , drop = FALSE], Y[i])$coefficients; b[is.na(b)] <- 0; M[, k] <- Xd %*% b }; M }

# tilted (augmented) estimates and influence functions for all beta in BALL; M = NULL -> unaugmented IPW
tilt <- function(Y, A, ps, K, M = NULL) {
  n <- length(Y); h <- 1 / rowSums(1 / ps); H <- exp(outer(log(h), BALL)); hb <- colMeans(H)
  if (is.null(M)) M <- matrix(0, n, K)
  r <- Y - M[cbind(seq_len(n), A)]; mu <- matrix(0, K, length(BALL)); IF <- vector("list", K)
  for (k in 1:K) {
    nu <- colSums(H * M[, k]) / colSums(H)
    ik <- A == k; Wk <- H[ik, , drop = FALSE] / ps[ik, k]; rho <- colSums(Wk * r[ik]) / colSums(Wk); mu[k, ] <- nu + rho
    I1 <- sweep(H * outer(M[, k], nu, "-"), 2, hb, "/")
    I2 <- matrix(0, n, length(BALL)); I2[ik, ] <- sweep(Wk * outer(r[ik], rho, "-"), 2, colSums(Wk) / n, "/")
    IF[[k]] <- I1 + I2
  }
  list(tau = sweep(mu[-1, , drop = FALSE], 2, mu[1, ]), IF = IF)
}
pet_select <- function(tl, j, kappa) {
  n <- nrow(tl$IF[[1]]); D <- tl$IF[[j + 1]] - tl$IF[[1]]
  target <- kappa * sum(D[, 1]^2) / n^2                              # BALL[1] = 0: IPW / AIPW variance
  res <- matrix(NA, length(QC), 3)
  for (qi in seq_along(QC)) {
    v <- e <- numeric(length(B1C))
    for (bi in seq_along(B1C)) { idx <- match(GRIDS[[bi]], BALL); a <- PALPHA[[bi]][[qi]]
      e[bi] <- sum(tl$tau[j, idx] * a); v[bi] <- sum((D[, idx] %*% a)^2) / n^2 }
    ok <- which(v <= target); b <- if (length(ok)) min(ok) else length(B1C); res[qi, ] <- c(e[b], v[b], v[b] <= target)
  }
  qs <- if (any(res[, 3] == 1)) max(which(res[, 3] == 1)) else which.min(res[, 2])
  c(est = res[qs, 1], se = sqrt(res[qs, 2]), q = QC[qs])
}
pet_all <- function(Y, A, ps, K, M = NULL) {
  tl <- tilt(Y, A, ps, K, M); kap <- if (is.null(M)) 5 / 6 else 1
  t(sapply(seq_len(K - 1), function(j) pet_select(tl, j, kap)))
}

if (DESIGN == "main") {
  Sys.setenv(SIM_SOURCE_ONLY = "1"); source("ow_vs_iptw_sim.R"); G0 <- grid; SEED0 <- 20261115; MINARM <- 5
  one <- function(g) {
    d <- gen_data(g$n, g$K, g$gamma, g$lambda); K <- g$K; if (any(tabulate(d$A, K) < MINARM)) return(NULL)
    ps <- fit_gps(d, K, flag = g$lambda > 0); X <- as.matrix(d[, paste0("X", 1:5)]); Xl <- cbind(1, X)
    list(A = d$A, K = K, ps = ps, XD = list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2)),
         Y = list(hom = d$Y_hom, het = d$Y_het))
  }
} else {
  if (DESIGN == "heldout") { Sys.setenv(HELDOUT_SOURCE_ONLY = "1"); source("heldout_sim.R"); G0 <- grid; G0$psmis <- 0; SEED0 <- 40000 } else {
    src <- readLines("stress_sim.R"); i <- grep("^RNGkind", src); eval(parse(text = src[1:(i - 1)])); G0 <- sgrid; SEED0 <- 50000 }
  one <- function(g) {
    d <- gen_data(g$n, g$K, g$gamma, g$pv); K <- g$K; if (any(tabulate(d$A, K) < 10)) return(NULL)
    ps <- fit_gps(d$X, d$A, K, g$pv); X <- d$X
    list(A = d$A, K = K, ps = ps, XD = list(lin = cbind(1, X), flex = cbind(1, X, X[, 1:4]^2, X[, 2] * X[, 5])), Y = d$Y)
  }
}
run_rep <- function(g) {
  d <- one(g); if (is.null(d)) return(NULL); out <- list()
  for (y in names(d$Y)) {
    r <- pet_all(d$Y[[y]], d$A, d$ps, d$K)
    out[[length(out) + 1]] <- data.frame(method = "PET", om = "none", effect = y, arm = 2:d$K, est = r[, "est"], se = r[, "se"], q = r[, "q"])
    for (om in names(d$XD)) {
      M <- fit_M(d$Y[[y]], d$A, d$XD[[om]], d$K); r <- pet_all(d$Y[[y]], d$A, d$ps, d$K, M)
      out[[length(out) + 1]] <- data.frame(method = "AIPW-PET", om = om, effect = y, arm = 2:d$K, est = r[, "est"], se = r[, "se"], q = r[, "q"])
      out[[length(out) + 1]] <- data.frame(method = "GCOMP-check", om = om, effect = y, arm = 2:d$K,
                                           est = colMeans(M)[-1] - colMeans(M)[1], se = NA, q = NA)
    }
  }
  do.call(rbind, out)
}
if (Sys.getenv("PET_TIME") == "1") { set.seed(2); g <- G0[nrow(G0), ]; if (!is.null(G0$psmis)) PSMIS <- g$psmis
  t0 <- Sys.time(); o <- run_rep(g); print(Sys.time() - t0); print(o); quit() }
OUT <- file.path("results_pet", DESIGN); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
for (s in seq_len(nrow(G0))) {
  g <- G0[s, ]; f <- file.path(OUT, sprintf("raw_s%02d.rds", s)); if (file.exists(f)) next
  if (!is.null(g$psmis)) PSMIS <- g$psmis
  set.seed(SEED0 + s)
  res <- mclapply(seq_len(NREP_P), function(r) { o <- tryCatch(run_rep(g), error = function(e) NULL); if (!is.null(o)) o$rep <- r; o },
                  mc.cores = CORES_P, mc.set.seed = TRUE)
  res <- do.call(rbind, res); res$scenario <- s; saveRDS(res, f)
  cat(sprintf("%s scenario %d done %.1f min\n", DESIGN, s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
