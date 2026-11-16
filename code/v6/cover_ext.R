############################################################################
# Extended sandwich-coverage study for ROW (balance rule, cap 8K^2) and AIPW (lambda = Inf),
# using the R package implementation (ROWpkg/R/row.R: M-estimation sandwich with analytical derivatives).
# Designs: main (poor overlap, n = 500 / 2000, all outcome types and both outcome models),
#          held-out (all 8 scenarios), stress test (all 16 scenarios).
# Usage: Rscript cover_ext.R <design> [n_rep] [cores]
############################################################################
args <- commandArgs(trailingOnly = TRUE)
DESIGN <- args[1]; NREP_C <- if (length(args) >= 2) as.integer(args[2]) else 300L
CORES_C <- if (length(args) >= 3) as.integer(args[3]) else 2L
suppressPackageStartupMessages({ library(nnet); library(parallel) })
source("../../R/row.R")
OUT <- file.path("results_cover", DESIGN); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (DESIGN == "main") {
  Sys.setenv(SIM_SOURCE_ONLY = "1"); source("ow_vs_iptw_sim.R")
  G0 <- grid[grid$overlap == "poor", ]; tr <- read.csv("true_estimands.csv")
  SC <- lapply(seq_len(nrow(G0)), function(i) { g <- G0[i, ]
    t <- tr[tr$K == g$K & tr$overlap == g$overlap & tr$positivity == g$positivity, ]
    list(id = g$scenario, gen = function() { d <- gen_data(g$n, g$K, g$gamma, g$lambda)
      X <- as.matrix(d[, paste0("X", 1:5)]); Z <- cbind(1, X); if (g$lambda > 0) Z <- cbind(Z, F = as.numeric(X[, 1] > 1))
      list(A = d$A, K = g$K, Z = Z, XO = list(lin = cbind(1, X), flex = cbind(1, X, X[, 1:3]^2)),
           Y = list(cont_hom = d$Y_hom, cont_het = d$Y_het)) },
      truth = list(cont_hom = t$ATE_hom, cont_het = t$ATE_het)) })
} else {
  if (DESIGN == "heldout") { Sys.setenv(HELDOUT_SOURCE_ONLY = "1"); source("heldout_sim.R"); G0 <- grid; tr <- read.csv("results_heldout/truth.csv")
    G0$psmis <- 0 } else {
    src <- readLines("stress_sim.R"); i <- grep("^RNGkind", src); eval(parse(text = src[1:(i - 1)])); G0 <- sgrid; tr <- read.csv("results_stress/truth.csv") }
  SC <- lapply(seq_len(nrow(G0)), function(i) { g <- G0[i, ]; t <- tr[tr$scenario == g$scenario, ]
    list(id = g$scenario, psmis = g$psmis, gen = function() { d <- gen_data(g$n, g$K, g$gamma, g$pv); X <- d$X
      Z <- cbind(1, X); if (g$pv) Z <- cbind(Z, F = as.numeric(X[, 2] > 1.2))
      list(A = d$A, K = g$K, Z = Z, XO = list(lin = cbind(1, X), flex = cbind(1, X, X[, 1:4]^2, X[, 2] * X[, 5])), Y = d$Y) },
      truth = lapply(setNames(nm = setdiff(names(t), c("scenario", "arm", "ATO_cont_het", "share_minps_lt05", "share_minps_lt01"))), function(v) t[[v]])) })
}

one <- function(sc) {
  if (!is.null(sc$psmis)) PSMIS <<- sc$psmis
  d <- sc$gen(); K <- d$K; if (any(tabulate(d$A, K) < 10)) return(NULL)
  G <- .fit_gps(d$A, d$Z, K); ps <- .gps_from_gamma(G, d$Z); ea <- ps[cbind(seq_along(d$A), d$A)]
  Fm <- .balance_features(d$Z); cg <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)
  smd <- vapply(cg, function(cc) .max_smd(Fm, d$A, 1 / (ea + 1 / (cc * K^2)), K), numeric(1))
  ok <- which(smd <= 0.1); csel <- if (length(ok)) cg[min(ok)] else cg[which.min(smd)]
  out <- list()
  for (y in names(d$Y)) for (om in names(d$XO)) for (m in c("ROW", "AIPW")) {
    r <- .row_core(d$Y[[y]], d$A, d$Z, d$XO[[om]], G, K, if (m == "ROW") csel else Inf)
    out[[length(out) + 1]] <- data.frame(method = m, effect = y, om = om, arm = 2:K, est = r$est, se = r$se, se_naive = r$se_naive,
                                         truth = sc$truth[[y]], lam = csel)
  }
  do.call(rbind, out)
}

if (Sys.getenv("COVER_TIME") == "1") { set.seed(1); t0 <- Sys.time(); o <- one(SC[[length(SC)]]); print(Sys.time() - t0); print(head(o)); quit() }
RNGkind("L'Ecuyer-CMRG"); t0 <- Sys.time()
for (sc in SC) {
  f <- file.path(OUT, sprintf("raw_s%02d.rds", sc$id)); if (file.exists(f)) next
  set.seed(60000 + sc$id + 100 * match(DESIGN, c("main", "heldout", "stress")))
  res <- do.call(rbind, mclapply(seq_len(NREP_C), function(r) { o <- tryCatch(one(sc), error = function(e) NULL); if (!is.null(o)) o$rep <- r; o },
                                 mc.cores = CORES_C, mc.set.seed = TRUE))
  res$scenario <- sc$id; saveRDS(res, f)
  cat(sprintf("%s scenario %d done %.1f min\n", DESIGN, sc$id, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
