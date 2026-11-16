# Real-data illustration: NCDS (PSweight), 3 education levels, outcome log hourly wage
source("row_sandwich.R")
load("psw/data/NCDS.rda"); d <- NCDS
lev <- c("None", "O/eq", ">=A/eq"); d$A <- match(d$Dmult, lev); K <- 3; n <- nrow(d); Y <- d$wage
Z <- model.matrix(~ white + maemp + factor(scht) + qmab + qmab2 + qvab + qvab2 + paed_u + maed_u + agepa + agema + sib_u, d)
f <- suppressWarnings(nnet::multinom(factor(d$A) ~ Z - 1, trace = FALSE, maxit = 2000, reltol = 1e-12))
G <- unname(as.matrix(coef(f))); ps <- gps_from_gamma(G, Z); ea <- ps[cbind(1:n, d$A)]
minps <- apply(ps, 1, min)
ess <- function(w) 100 * sum(sapply(1:K, function(k) { x <- w[d$A == k]; sum(x)^2 / sum(x^2) })) / n
diag <- data.frame(n = n, min_gps_q01 = quantile(minps, .01), share_minps_lt05 = mean(minps < .05), share_minps_lt10 = mean(minps < .10),
  ess_iptw = ess(1 / ea), ess_ow = ess((1 / ea) / rowSums(1 / ps)), max_w_ratio = max(1 / ea) / mean(1 / ea), max_w = max(1 / ea))
print(round(diag, 3))
hajek <- function(w) { mu <- sapply(1:K, function(k) { i <- d$A == k; sum(w[i] * Y[i]) / sum(w[i]) })
  v <- sapply(1:K, function(k) { i <- d$A == k; sum(w[i]^2 * (Y[i] - mu[k])^2) / sum(w[i])^2 })
  data.frame(arm = 2:K, est = mu[-1] - mu[1], se = sqrt(v[-1] + v[1])) }
M <- matrix(NA, n, K); for (k in 1:K) { i <- d$A == k; b <- lm.fit(Z[i, ], Y[i])$coefficients; b[is.na(b)] <- 0; M[, k] <- Z %*% b }
aug <- function(w) { r <- Y - M[cbind(1:n, d$A)]; corr <- numeric(K); IC <- matrix(0, n, K)
  for (k in 1:K) { i <- d$A == k; corr[k] <- sum(w[i] * r[i]) / sum(w[i]); IC[i, k] <- n * w[i] * (r[i] - corr[k]) / sum(w[i]) }
  mu <- colMeans(M) + corr
  data.frame(arm = 2:K, est = mu[-1] - mu[1], se = sapply(2:K, function(k) { dA <- M[, k] - M[, 1]; sqrt(sum(((dA - mean(dA)) + IC[, k] - IC[, 1])^2)) / n })) }
res <- list()
res$IPTW <- hajek(1 / ea)
res$OW <- hajek((1 / ea) / rowSums(1 / ps))
res$GCOMP <- aug(rep(1, n))
res$AIPW <- aug(1 / ea)
tw <- (tabulate(d$A, K)[d$A] / n) / ea
q1 <- quantile(tw, c(.01, .99)); res$AIPW_T1_99 <- aug(pmin(pmax(tw, q1[1]), q1[2]))
q5 <- quantile(tw, c(.05, .95)); res$AIPW_T5_95 <- aug(pmin(pmax(tw, q5[1]), q5[2]))
# balance rule (outcome-blind): smallest lambda <= 8K^2 with max SMD <= 0.1, features = covariates + squares of numeric ones
num <- c("qmab", "qmab2", "qvab", "qvab2", "paed_u", "maed_u", "agepa", "agema", "sib_u")
FMb <- cbind(Z[, -1], as.matrix(d[, num])^2)
msmd <- function(w) { mu <- colMeans(FMb); s <- apply(FMb, 2, sd); s[s == 0] <- 1
  max(sapply(1:K, function(k) { i <- d$A == k; max(abs(colSums(FMb[i, ] * w[i]) / sum(w[i]) - mu) / s) })) }
cg <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)
sm <- sapply(cg, function(c) msmd(1 / (ea + 1 / (c * K^2))))
lam_sel <- if (any(sm <= 0.1)) cg[min(which(sm <= 0.1))] else cg[which.min(sm)]
cat("max SMD by lambda/K^2:", paste(cg, round(sm, 3), sep = ":", collapse = "  "), "\nselected lambda/K^2 =", lam_sel, "\n")
cat("max SMD: unweighted", round(msmd(rep(1, n)), 3), " IPTW", round(msmd(1 / ea), 3), "\n")
rb <- row_sandwich(Y, d$A, Z, Z, G, K, lam_sel); res$ROW <- data.frame(arm = 2:K, est = rb$est, se = rb$se_sand)
out <- do.call(rbind, lapply(names(res), function(m) data.frame(method = m, contrast = paste(lev[res[[m]]$arm], "vs", lev[1]), res[[m]][, c("est", "se")])))
out$lo <- out$est - 1.96 * out$se; out$hi <- out$est + 1.96 * out$se; out$lam_sel <- lam_sel
print(transform(out, est = round(est, 3), se = round(se, 3), lo = round(lo, 3), hi = round(hi, 3)), row.names = FALSE)
CG <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, Inf)
path <- do.call(rbind, lapply(CG, function(c) data.frame(lamK2 = c, row_sandwich(Y, d$A, Z, Z, G, K, c))))
path$contrast <- paste(lev[path$arm], "vs", lev[1]); print(path[, c("lamK2", "contrast", "est", "se_sand")], row.names = FALSE, digits = 3)
write.csv(out, "ncds_estimates.csv", row.names = FALSE); write.csv(path, "ncds_path.csv", row.names = FALSE); write.csv(diag, "ncds_diag.csv", row.names = FALSE)
