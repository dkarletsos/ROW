# Clinical illustration: Rotterdam breast cancer cohort (survival::rotterdam), 3 adjuvant-therapy groups
suppressPackageStartupMessages(library(survival)); source("../../R/row.R")
d <- rotterdam; d <- d[!(d$hormon == 1 & d$chemo == 1) & d$nodes > 0, ]   # node-positive; "both" (n = 28) excluded
d$trt <- factor(ifelse(d$chemo == 1, "chemo", ifelse(d$hormon == 1, "hormonal", "none")), levels = c("none", "hormonal", "chemo"))
d$y5 <- ifelse(d$dtime <= 1826 & d$death == 1, 1, ifelse(d$dtime > 1826, 0, NA)); d <- d[!is.na(d$y5), ]
d$lpgr <- log1p(d$pgr); d$ler <- log1p(d$er); d$size2 <- as.numeric(d$size == "20-50"); d$size3 <- as.numeric(d$size == ">50")
d$lnodes <- log1p(d$nodes); d$year <- d$year - 1985
f  <- ~ age + meno + size2 + size3 + grade + lnodes + lpgr + ler + year
fo <- ~ age + meno + size2 + size3 + grade + lnodes + lpgr + ler + year + I(age^2) + I(lnodes^2) + I(lpgr^2) + I(ler^2)
cat("n =", nrow(d), "\n"); print(table(d$trt)); print(round(tapply(d$y5, d$trt, mean), 3))
fit <- row_ate(d, "trt", "y5", f, fo, reference = "none"); print(fit); print(fit$balance_path)
ps <- fit$ps; A <- as.integer(d$trt); K <- 3; n <- nrow(d); ea <- ps[cbind(1:n, A)]; minps <- apply(ps, 1, min)
diag <- data.frame(n = n, share_minps_lt05 = mean(minps < .05), share_minps_lt01 = mean(minps < .01), max_iptw = max(1 / ea),
                   max_ratio = max(1 / ea) / mean(1 / ea), lam_sel = fit$lambda_K2, max_smd = fit$max_smd)
print(round(diag, 3))
Y <- d$y5; Xo <- model.matrix(fo, d)
M <- matrix(NA, n, K); for (k in 1:K) { i <- A == k; b <- lm.fit(Xo[i, ], Y[i])$coefficients; b[is.na(b)] <- 0; M[, k] <- Xo %*% b }
hajek <- function(w) { mu <- sapply(1:K, function(k) { i <- A == k; sum(w[i] * Y[i]) / sum(w[i]) })
  v <- sapply(1:K, function(k) { i <- A == k; sum(w[i]^2 * (Y[i] - mu[k])^2) / sum(w[i])^2 }); data.frame(arm = 2:K, est = mu[-1] - mu[1], se = sqrt(v[-1] + v[1])) }
aug <- function(w, h = rep(1, n)) { r <- Y - M[cbind(1:n, A)]; corr <- numeric(K); IC <- matrix(0, n, K)
  for (k in 1:K) { i <- A == k; corr[k] <- sum(w[i] * r[i]) / sum(w[i]); IC[i, k] <- n * w[i] * (r[i] - corr[k]) / sum(w[i]) }
  mu <- colSums(M * h) / sum(h) + corr
  data.frame(arm = 2:K, est = mu[-1] - mu[1], se = sapply(2:K, function(k) { dA <- (M[, k] - M[, 1]) * h / mean(h); sqrt(sum(((dA - mean(dA)) + IC[, k] - IC[, 1])^2)) / n })) }
h <- 1 / rowSums(1 / ps); tw <- (tabulate(A, K)[A] / n) / ea; tr <- function(lo, hi) { q <- quantile(tw, c(lo, hi)); pmin(pmax(tw, q[1]), q[2]) }
res <- list(IPTW = hajek(1 / ea), OW = hajek(h / ea), AOW = aug(h / ea, h), `G-computation` = aug(rep(0, n) + 1e-12),
            AIPW = aug(1 / ea), `AIPW-T1/99` = aug(tr(.01, .99)), `AIPW-T5/95` = aug(tr(.05, .95)),
            ROW = data.frame(arm = 2:K, est = fit$estimates$estimate, se = fit$estimates$se))
res$`G-computation`$se <- NA
out <- do.call(rbind, lapply(names(res), function(m) data.frame(method = m, contrast = paste(levels(d$trt)[res[[m]]$arm], "vs none"), res[[m]][, c("est", "se")])))
print(transform(out, est = round(est, 3), se = round(se, 3)), row.names = FALSE)
p <- row_path(d, "trt", "y5", f, fo, reference = "none"); p$estimate <- round(p$estimate, 3); p$se <- round(p$se, 3); print(p)
write.csv(out, "rotterdam_estimates.csv", row.names = FALSE); write.csv(diag, "rotterdam_diag.csv", row.names = FALSE); write.csv(p, "rotterdam_path.csv", row.names = FALSE)
