############################################################################
# ROW (arm-specific regularised weights) with M-estimation sandwich SEs
#   w_k(x) = 1 / (e_k(x) + 1/lambda),  lambda = c * K^2  (c = Inf -> AIPW)
# Stacked estimating equations: multinomial-logit score (gamma), per-arm OLS (beta_k),
# plug-in mean nu_k, weighted residual mean rho_k; mu_k = nu_k + rho_k.
# Block-triangular system -> sandwich = variance of the influence function below.
# d w_k / d gamma_j = -w_k^2 * e_k * (1{k=j} - e_j) * z   (analytic)
############################################################################
gps_from_gamma <- function(G, Z) { lp <- cbind(0, Z %*% t(G)); e <- exp(lp - apply(lp, 1, max)); e / rowSums(e) }

fit_gamma <- function(d, K, flag) {
  Z <- if (flag) cbind(1, as.matrix(d[, paste0("X", 1:5)]), as.numeric(d$X1 > 1)) else cbind(1, as.matrix(d[, paste0("X", 1:5)]))
  f <- suppressWarnings(nnet::multinom(factor(d$A, levels = 1:K) ~ Z - 1, trace = FALSE, maxit = 1000, reltol = 1e-12))
  G <- coef(f); if (is.null(dim(G))) G <- matrix(G, 1)
  list(G = unname(as.matrix(G)), Z = Z)
}

row_sandwich <- function(Y, A, Z, Xo, G, K, c) {
  n <- length(Y); p <- ncol(Z); q <- ncol(Xo)
  ps <- gps_from_gamma(G, Z)
  S <- matrix(0, n, (K - 1) * p); Info <- matrix(0, (K - 1) * p, (K - 1) * p)
  for (j in 2:K) {
    cj <- ((j - 2) * p + 1):((j - 1) * p); S[, cj] <- Z * ((A == j) - ps[, j])
    for (l in 2:K) { cl <- ((l - 2) * p + 1):((l - 1) * p)
      Info[cj, cl] <- crossprod(Z * (ps[, j] * ((j == l) - ps[, l])), Z) / n }
  }
  IFg <- S %*% solve(Info)
  IFmu <- IFmu_naive <- matrix(0, n, K); mu <- numeric(K)
  for (k in 1:K) {
    i <- A == k; Xi <- Xo[i, , drop = FALSE]
    XtX <- crossprod(Xi) / n
    b <- solve(crossprod(Xi), crossprod(Xi, Y[i])); m <- as.vector(Xo %*% b); r <- Y - m
    IFb <- matrix(0, n, q); IFb[i, ] <- (Xi * r[i]) %*% solve(XtX)
    wk <- if (is.infinite(c)) 1 / ps[, k] else 1 / (ps[, k] + 1 / (c * K^2))
    w <- wk * i
    nu <- mean(m); rho <- sum(w * r) / sum(w); mu[k] <- nu + rho
    # derivative of mean(I w (r - rho)) w.r.t. gamma (arms 2..K blocks)
    base <- -wk^2 * ps[, k] * i * (r - rho)
    Gd <- unlist(lapply(2:K, function(j) colMeans(Z * (base * ((k == j) - ps[, j])))))
    Wm <- mean(w)
    IFnu  <- (m - nu) + as.vector(IFb %*% colMeans(Xo))
    IFrho <- (w * (r - rho) + as.vector(IFg %*% Gd) - as.vector(IFb %*% colMeans(Xo * w))) / Wm
    IFmu[, k] <- IFnu + IFrho
    IFmu_naive[, k] <- (m - nu) + w * (r - rho) / Wm
  }
  data.frame(arm = 2:K, est = mu[-1] - mu[1],
    se_naive = sapply(2:K, function(k) sqrt(sum((IFmu_naive[, k] - IFmu_naive[, 1])^2)) / n),
    se_sand  = sapply(2:K, function(k) sqrt(sum((IFmu[, k] - IFmu[, 1])^2)) / n))
}
