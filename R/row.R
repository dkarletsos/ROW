## ROW: regularised outcome-weighted estimation of the average treatment effect (ATE)
## for comparisons of two or more treatments.

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------
.gps_from_gamma <- function(G, Z) {
  lp <- cbind(0, Z %*% t(G))
  e <- exp(lp - apply(lp, 1, max))
  e / rowSums(e)
}

.fit_gps <- function(A, Z, K) {
  f <- suppressWarnings(nnet::multinom(factor(A, levels = seq_len(K)) ~ Z - 1,
                                       trace = FALSE, maxit = 2000, reltol = 1e-12))
  G <- stats::coef(f)
  if (is.null(dim(G))) G <- matrix(G, 1)
  unname(as.matrix(G))
}

.balance_features <- function(Z) {
  X <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  nonbin <- apply(X, 2, function(v) length(unique(v)) > 2)
  if (any(nonbin)) {
    sq <- X[, nonbin, drop = FALSE]^2
    colnames(sq) <- paste0(colnames(X)[nonbin], "^2")
    X <- cbind(X, sq)
  }
  X
}

.max_smd <- function(Fm, A, w, K) {
  mu <- colMeans(Fm); s <- apply(Fm, 2, stats::sd); s[s == 0] <- 1
  max(vapply(seq_len(K), function(k) {
    i <- A == k
    max(abs(colSums(Fm[i, , drop = FALSE] * w[i]) / sum(w[i]) - mu) / s)
  }, numeric(1)))
}

# ROW estimate and influence functions (naive and sandwich) for lambda = c * K^2
.row_core <- function(Y, A, Z, Xo, G, K, c) {
  n <- length(Y); p <- ncol(Z); q <- ncol(Xo)
  ps <- .gps_from_gamma(G, Z)
  S <- matrix(0, n, (K - 1) * p); Info <- matrix(0, (K - 1) * p, (K - 1) * p)
  for (j in 2:K) {
    cj <- ((j - 2) * p + 1):((j - 1) * p)
    S[, cj] <- Z * ((A == j) - ps[, j])
    for (l in 2:K) {
      cl <- ((l - 2) * p + 1):((l - 1) * p)
      Info[cj, cl] <- crossprod(Z * (ps[, j] * ((j == l) - ps[, l])), Z) / n
    }
  }
  IFg <- S %*% solve(Info)
  IFmu <- IFnv <- matrix(0, n, K); mu <- numeric(K)
  for (k in seq_len(K)) {
    i <- A == k; Xi <- Xo[i, , drop = FALSE]
    XtX <- crossprod(Xi) / n
    b <- qr.solve(crossprod(Xi), crossprod(Xi, Y[i]))
    m <- as.vector(Xo %*% b); r <- Y - m
    IFb <- matrix(0, n, q); IFb[i, ] <- (Xi * r[i]) %*% qr.solve(XtX)
    wk <- if (is.infinite(c)) 1 / ps[, k] else 1 / (ps[, k] + 1 / (c * K^2))
    w <- wk * i
    nu <- mean(m); rho <- sum(w * r) / sum(w); mu[k] <- nu + rho
    base <- -wk^2 * ps[, k] * i * (r - rho)
    Gd <- unlist(lapply(2:K, function(j) colMeans(Z * (base * ((k == j) - ps[, j])))))
    Wm <- mean(w)
    IFnu  <- (m - nu) + as.vector(IFb %*% colMeans(Xo))
    IFrho <- (w * (r - rho) + as.vector(IFg %*% Gd) - as.vector(IFb %*% colMeans(Xo * w))) / Wm
    IFmu[, k] <- IFnu + IFrho
    IFnv[, k] <- (m - nu) + w * (r - rho) / Wm
  }
  list(mu = mu,
       est = mu[-1] - mu[1],
       se = vapply(2:K, function(k) sqrt(sum((IFmu[, k] - IFmu[, 1])^2)) / n, numeric(1)),
       se_naive = vapply(2:K, function(k) sqrt(sum((IFnv[, k] - IFnv[, 1])^2)) / n, numeric(1)))
}

.prepare <- function(data, treatment, outcome, ps_formula, outcome_formula, reference) {
  stopifnot(is.data.frame(data), treatment %in% names(data), outcome %in% names(data))
  vars <- unique(c(treatment, outcome, all.vars(ps_formula), all.vars(outcome_formula)))
  d <- stats::na.omit(data[, vars, drop = FALSE])
  trt <- factor(d[[treatment]])
  if (!is.null(reference)) trt <- stats::relevel(trt, ref = as.character(reference))
  K <- nlevels(trt)
  if (K < 2) stop("the treatment must have at least two levels")
  A <- as.integer(trt)
  if (any(tabulate(A, K) < 5)) stop("each treatment level needs at least 5 observations")
  Z  <- stats::model.matrix(ps_formula, d)
  Xo <- stats::model.matrix(outcome_formula, d)
  list(d = d, Y = d[[outcome]], A = A, K = K, levels = levels(trt), Z = Z, Xo = Xo)
}

# ---------------------------------------------------------------------------
# exported functions
# ---------------------------------------------------------------------------

row_weights <- function(ps, A, lambda) {
  ps <- as.matrix(ps)
  ea <- ps[cbind(seq_along(A), A)]
  if (is.infinite(lambda)) 1 / ea else 1 / (ea + 1 / lambda)
}

row_ate <- function(data, treatment, outcome, ps_formula, outcome_formula = ps_formula,
                    lambda = "balance", cap = 8, smd_threshold = 0.1,
                    grid = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128),
                    reference = NULL, level = 0.95) {
  pr <- .prepare(data, treatment, outcome, ps_formula, outcome_formula, reference)
  K <- pr$K
  G <- .fit_gps(pr$A, pr$Z, K)
  ps <- .gps_from_gamma(G, pr$Z)
  ea <- ps[cbind(seq_along(pr$A), pr$A)]
  Fm <- .balance_features(pr$Z)
  smd_path <- NULL
  if (identical(lambda, "balance")) {
    cg <- sort(grid[grid <= cap])
    if (!length(cg)) stop("no grid value at or below 'cap'")
    smd_path <- vapply(cg, function(cc) .max_smd(Fm, pr$A, 1 / (ea + 1 / (cc * K^2)), K), numeric(1))
    ok <- which(smd_path <= smd_threshold)
    c_sel <- if (length(ok)) cg[min(ok)] else cg[which.min(smd_path)]
    smd_path <- data.frame(lambda_K2 = cg, max_smd = smd_path)
  } else {
    if (!is.numeric(lambda) || length(lambda) != 1 || lambda <= 0) stop("'lambda' must be \"balance\" or a positive number (on the lambda / K^2 scale)")
    c_sel <- lambda
  }
  r <- .row_core(pr$Y, pr$A, pr$Z, pr$Xo, G, K, c_sel)
  z <- stats::qnorm(1 - (1 - level) / 2)
  est <- data.frame(contrast = paste(pr$levels[-1], "vs", pr$levels[1]),
                    estimate = r$est, se = r$se,
                    lower = r$est - z * r$se, upper = r$est + z * r$se,
                    se_naive = r$se_naive, stringsAsFactors = FALSE)
  w <- if (is.infinite(c_sel)) 1 / ea else 1 / (ea + 1 / (c_sel * K^2))
  structure(list(estimates = est, means = stats::setNames(r$mu, pr$levels),
                 lambda_K2 = c_sel, lambda = c_sel * K^2, K = K, n = nrow(pr$d),
                 levels = pr$levels, max_smd = .max_smd(Fm, pr$A, w, K),
                 balance_path = smd_path, weights = w, ps = ps, level = level,
                 call = match.call()),
            class = "row")
}

row_path <- function(data, treatment, outcome, ps_formula, outcome_formula = ps_formula,
                     grid = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, Inf), reference = NULL) {
  pr <- .prepare(data, treatment, outcome, ps_formula, outcome_formula, reference)
  G <- .fit_gps(pr$A, pr$Z, pr$K)
  do.call(rbind, lapply(grid, function(cc) {
    r <- .row_core(pr$Y, pr$A, pr$Z, pr$Xo, G, pr$K, cc)
    data.frame(lambda_K2 = cc, contrast = paste(pr$levels[-1], "vs", pr$levels[1]),
               estimate = r$est, se = r$se, stringsAsFactors = FALSE)
  }))
}

print.row <- function(x, digits = 3, ...) {
  cat("ROW: regularised outcome-weighted estimate of the average treatment effect\n")
  cat(sprintf("n = %d, %d treatment levels (reference: %s)\n", x$n, x$K, x$levels[1]))
  rule <- if (is.null(x$balance_path)) "fixed" else "balance rule"
  cat(sprintf("lambda = %g K^2 = %g (%s); max SMD of weighted arms vs full population = %.3f\n",
              x$lambda_K2, x$lambda, rule, x$max_smd))
  cat(sprintf("Standard errors: M-estimation sandwich; %g%% confidence intervals\n\n", 100 * x$level))
  e <- x$estimates; e[, -1] <- round(e[, -1], digits)
  print(e[, c("contrast", "estimate", "se", "lower", "upper")], row.names = FALSE)
  invisible(x)
}

coef.row <- function(object, ...) stats::setNames(object$estimates$estimate, object$estimates$contrast)

confint.row <- function(object, parm, level = object$level, ...) {
  z <- stats::qnorm(1 - (1 - level) / 2); e <- object$estimates
  out <- cbind(e$estimate - z * e$se, e$estimate + z * e$se)
  dimnames(out) <- list(e$contrast, paste0(format(100 * c((1 - level) / 2, 1 - (1 - level) / 2)), " %"))
  out
}
