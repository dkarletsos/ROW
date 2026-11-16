# Coverage of naive vs sandwich SEs for ROW (arm-specific), poor overlap, n = 2000
Sys.setenv(SIM_SOURCE_ONLY = "1"); source("ow_vs_iptw_sim.R"); source("row_sandwich.R")
args <- commandArgs(trailingOnly = TRUE)
N_REP <- if (length(args) >= 1) as.integer(args[1]) else 1000L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 2L
OUTDIR <- if (length(args) >= 3) args[3] else "results_rowsand"
dir.create(OUTDIR, showWarnings = FALSE)
CS <- c(2, 8, Inf); SCEN <- which(grid$n == 2000 & grid$overlap == "poor")
one <- function(n, K, gamma, lambda) {
  d <- gen_data(n, K, gamma, lambda); if (any(tabulate(d$A, K) < 5)) return(NULL)
  g <- fit_gamma(d, K, lambda > 0); Xl <- cbind(1, as.matrix(d[, paste0("X", 1:5)]))
  XD <- list(lin = Xl, quad = cbind(Xl, d$X1^2, d$X2^2, d$X3^2)); out <- list()
  for (om in names(XD)) for (c in CS) for (y in c("hom", "het"))
    out[[length(out) + 1]] <- data.frame(c = c, om = om, effect = y, row_sandwich(d[[paste0("Y_", y)]], d$A, g$Z, XD[[om]], g$G, K, c))
  do.call(rbind, out)
}
RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
for (s in SCEN) {
  g <- grid[s, ]; f <- file.path(OUTDIR, sprintf("raw_s%02d.rds", s)); if (file.exists(f)) next
  set.seed(20261115 + s)
  res <- do.call(rbind, mclapply(seq_len(N_REP), function(r) { o <- tryCatch(one(g$n, g$K, g$gamma, g$lambda), error = function(e) NULL)
    if (!is.null(o)) o$rep <- r; o }, mc.cores = CORES, mc.set.seed = TRUE))
  saveRDS(cbind(g[rep(1, nrow(res)), c("scenario", "n", "K", "overlap", "positivity")], res, row.names = NULL), f)
  cat(sprintf("Scenario %2d done %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
