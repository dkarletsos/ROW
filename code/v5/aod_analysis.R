# Real-data illustration: AOD (twang), 3 treatment arms, outcome suf12 (substance use frequency, 12 months)
source("row_sandwich.R")
load("tw/data/AOD.RData"); d <- AOD
d$A <- as.integer(d$treat); K <- 3; n <- nrow(d)          # 1 = community (reference), 2 = metcbt5, 3 = scy
labs <- levels(d$treat)
cov <- c("illact", "crimjust", "subprob", "subdep", "white")
Z <- cbind(1, as.matrix(d[, cov]))
f <- suppressWarnings(nnet::multinom(factor(d$A) ~ Z - 1, trace = FALSE, maxit = 1000, reltol = 1e-12))
G <- unname(as.matrix(coef(f))); ps <- gps_from_gamma(G, Z)
ea <- ps[cbind(1:n, d$A)]; Y <- d$suf12
minps <- apply(ps, 1, min)
ess <- function(w) 100 * sum(sapply(1:K, function(k) { x <- w[d$A == k]; sum(x)^2 / sum(x^2) })) / n
diag <- data.frame(min_gps_q01 = quantile(minps, .01), share_minps_lt05 = mean(minps < .05), share_minps_lt10 = mean(minps < .10),
                   ess_iptw = ess(1 / ea), ess_ow = ess((1 / ea) / rowSums(1 / ps)), max_w_iptw_ratio = max(1/ea) / mean(1/ea))
print(round(diag, 3))
hajek <- function(w) { mu <- sapply(1:K, function(k) { i <- d$A == k; sum(w[i] * Y[i]) / sum(w[i]) })
  v <- sapply(1:K, function(k) { i <- d$A == k; sum(w[i]^2 * (Y[i] - mu[k])^2) / sum(w[i])^2 })
  data.frame(arm = 2:K, est = mu[-1] - mu[1], se = sqrt(v[-1] + v[1])) }
Xo <- Z
res <- list()
res$IPTW <- hajek(1 / ea)
res$OW <- hajek((1 / ea) / rowSums(1 / ps))
M <- matrix(NA, n, K); for (k in 1:K) { i <- d$A == k; b <- lm.fit(Xo[i, ], Y[i])$coefficients; M[, k] <- Xo %*% b }
aug <- function(w) { r <- Y - M[cbind(1:n, d$A)]; corr <- numeric(K); IC <- matrix(0, n, K)
  for (k in 1:K) { i <- d$A == k; corr[k] <- sum(w[i] * r[i]) / sum(w[i]); IC[i, k] <- n * w[i] * (r[i] - corr[k]) / sum(w[i]) }
  mu <- colMeans(M) + corr
  data.frame(arm = 2:K, est = mu[-1] - mu[1], se = sapply(2:K, function(k) { dA <- M[, k] - M[, 1]; sqrt(sum(((dA - mean(dA)) + IC[, k] - IC[, 1])^2)) / n })) }
res$GCOMP <- aug(rep(1, n))
res$AIPW <- aug(1 / ea)
tw <- (tabulate(d$A, K)[d$A] / n) / ea; q <- quantile(tw, c(.01, .99)); res$`AIPW-T1/99` <- aug(pmin(pmax(tw, q[1]), q[2]))
rs <- row_sandwich(Y, d$A, Z, Xo, G, K, 2); res$`ROW (2K², sandwich)` <- data.frame(arm = 2:K, est = rs$est, se = rs$se_sand)
out <- do.call(rbind, lapply(names(res), function(m) data.frame(method = m, contrast = paste(labs[res[[m]]$arm], "vs", labs[1]), res[[m]][, c("est", "se")])))
out$lo <- out$est - 1.96 * out$se; out$hi <- out$est + 1.96 * out$se
print(transform(out, est = round(est, 3), se = round(se, 3), lo = round(lo, 3), hi = round(hi, 3)), row.names = FALSE)
# lambda path (sandwich SE)
CG <- c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, Inf)
path <- do.call(rbind, lapply(CG, function(c) data.frame(lamK2 = c, row_sandwich(Y, d$A, Z, Xo, G, K, c))))
path$contrast <- paste(labs[path$arm], "vs", labs[1]); print(path[, c("lamK2", "contrast", "est", "se_sand")], row.names = FALSE, digits = 3)
write.csv(out, "aod_estimates.csv", row.names = FALSE); write.csv(path, "aod_path.csv", row.names = FALSE); write.csv(diag, "aod_diag.csv", row.names = FALSE)
# ---- balance rule (outcome-blind): smallest lambda <= 4K^2 with max SMD <= 0.1 (covariates + squares of continuous ones)
FMb <- cbind(as.matrix(d[, cov]), d$illact^2, d$crimjust^2, d$subprob^2, d$subdep^2)
msmd <- function(w) { mu <- colMeans(FMb); s <- apply(FMb, 2, sd)
  max(sapply(1:K, function(k) { i <- d$A == k; max(abs(colSums(FMb[i, ] * w[i]) / sum(w[i]) - mu) / s) })) }
cg <- c(0.125, 0.25, 0.5, 1, 2, 4)
sm <- sapply(cg, function(c) msmd(1 / (ea + 1 / (c * K^2))))
lam_sel <- if (any(sm <= 0.1)) cg[min(which(sm <= 0.1))] else cg[which.min(sm)]
cat("balance-rule lambda/K^2 =", lam_sel, " (max SMD by lambda:", round(sm, 3), ")\n")
rb <- row_sandwich(Y, d$A, Z, Xo, G, K, lam_sel)
rbo <- data.frame(method = "ROW-bal", contrast = paste(labs[rb$arm], "vs", labs[1]), est = rb$est, se = rb$se_sand)
rbo$lo <- rbo$est - 1.96 * rbo$se; rbo$hi <- rbo$est + 1.96 * rbo$se; print(rbo)
out2 <- rbind(out[!grepl("^ROW", out$method), ], rbo); out2$lam_sel <- lam_sel
write.csv(out2, "aod_estimates.csv", row.names = FALSE)
