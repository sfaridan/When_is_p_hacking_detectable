################################################################################
#  reproduce_kudrin_fig2_minimum.R
#
#  Reproduces the bottom row ("Minimum") of Kudrin (2025, "Testing for and Evaluating the
#  Extent of Selective Reporting", Fig. 2) and adds the projection test.
#
#  What is simulated (Kudrin Sec. 5; EKW 2025 Secs. 5.1.1 and 5.3.1):
#    covariate selection with K = 3 controls, N = 200, minimum approach;
#    panels h = 1, h = 2 and h ~ Pi_hat;  S = 1,000 t-statistics per meta-sample;
#    tau = 0, 0.1, ..., 1;  1,000 Monte Carlo repetitions per (panel, tau).
#    As in EKW/Kudrin, one million honest (all-controls) and one million minimum-p-value
#    t-statistics are simulated once per panel; each meta-sample draws S of them with
#    replacement, each from the p-hacked pool with probability tau.
#
#  Uses your code only:
#    reproduce_kudrin_simulations_functions.R  kudrin_draw()            the DGP
#    kudrin_tests.R                            kudrin_tests()           T_inf, T_H, T_H^hybrid
#    functions_all_weighted.R                  run_test_weight_fast()   projection test
#                                              create_basis_stable(), compute_epsilons()
#    Tests.R (EKW)                             only if CFG$run_ekw = TRUE
#
#  Outputs in OUT_DIR:
#    fig2_minimum_reproduction.pdf/.png  his formatting, with the projection test added
#    fig2_minimum_vs_published.pdf/.png  your curves (solid) over his published ones (faint)
#    fig2_minimum_published_redrawn.pdf  his published curves in this format (drawn first)
#    results_fig2_minimum.csv/.rds       rejection rates by panel, tau and test
#    comparison_with_published.csv       your curves minus his, with simulation-error flags
#  Finished cells are saved in OUT_DIR/cells, so an interrupted run resumes where it stopped.
#
#  Run once with QUICK <- TRUE (minutes) to check that everything works, then with FALSE.
################################################################################

rm(list = ls())

## ---- 0. Settings -------------------------------------------------------------------------

ROOT     <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable" #"C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"

CODE_DIR <- file.path(ROOT, "code")
OUT_DIR  <- file.path(ROOT, "simulations", "kudrin_fig2_minimum_J0_20_check")

QUICK <- FALSE   # TRUE: small run that checks the pipeline; FALSE: Kudrin's settings

CFG <- list(
  S              = 1000,                  # t-statistics per meta-sample (his Fig. 2)
  taus           = seq(0, 1, by = 0.1),   # extent of selective reporting
  reps           = 1000,                  # Monte Carlo repetitions per (panel, tau), as in his Fig. 2
  panels         = c("h1", "h2", "hPi"),  # h = 1, h = 2, h ~ Pi_hat
  alpha          = 0.05,
  
  # DGP, passed to kudrin_draw()
  n_pool         = 1e6,                   # studies per pool (EKW and Kudrin use one million)
  K              = 3,
  N              = 200,
  gamma_max      = 0.8,
  dgp            = "regression",
  # Pi_hat = Gamma(shape 1.4134, 1.2980). His footnote 13, and KUDRIN_GAMMA in your DGP file,
  # treat 1.2980 as a RATE (mean 1.09). His Fig. 2 is reproduced only when 1.2980 is a SCALE
  # (mean 1.83), presumably MATLAB's gamrnd(a, b) convention: in a Python mirror of this design,
  # T_inf at tau = 0.5 is 0.31 with the scale and 0.67 with the rate, against his 0.34.
  gamma_shape    = 1.4134,
  gamma_second   = 1.2980,
  gamma_param    = "scale",               # "scale" (matches his figure) or "rate" (his footnote)
  
  # Kudrin tests: your kudrin_tests()
  kudrin_boots   = 199,                   # T_inf bootstrap draws (his Monte Carlo used warp-speed)
  eig_floor      = 1e-8,                  # hybrid metric floor (your default). Extra values, e.g.
  # c(1e-8, 1e-4), add hybrid-only sensitivity columns.
  
  # Projection test: your run_test_weight_fast(), with the paper's settings
  run_projection = TRUE,
  num_coeffs     = 41,                    # J = 40
  J0             = 20,                    # weights w_j = min(1, (1/p)^(j - J0))
  p_weight       = 2,
  shift          = 1.96,
  L              = 6.5,
  numgrid        = 3000,
  proj_boots     = 100,
  proj_mode      = "full",                # "full": proj_boots bootstrap draws per meta-sample, as in
  # the paper; "warp": one draw per meta-sample and a pooled
  # critical value (Giacomini, Politis and White 2013)
  epsilon_U_pre  = NULL,                  # NULL: compute_epsilons(L, numgrid) once and cache it
  # (slow); or 0.00063, the value in the paper
  
  # EKW tests for the blue line (Tests.R)
  run_ekw        = FALSE,                 # FALSE: plot his published "Max Power of EKW Tests" line
  ekw_J          = 15,                    # bins on (0, 0.15], as in EKW (2025) Table 1
  
  # Computing
  n_cores        = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  seed           = 20261022,
  resume         = TRUE,
  run_selfchecks = TRUE
)

if (QUICK) {
  CFG$reps          <- 100
  CFG$n_pool        <- 2e5
  CFG$kudrin_boots  <- 99
  CFG$proj_mode     <- "warp"
  CFG$epsilon_U_pre <- 0.00063
  OUT_DIR           <- paste0(OUT_DIR, "_quick")
}


## ---- 1. Your code ------------------------------------------------------------------------

CODE_FILES <- file.path(CODE_DIR, c(if (CFG$run_ekw) "Tests.R",
                                    "reproduce_kudrin_simulations_functions.R",
                                    "functions_all_weighted.R",
                                    "kudrin_tests.R"))
missing_files <- CODE_FILES[!file.exists(CODE_FILES)]
if (length(missing_files) > 0) stop("Not found: ", paste(missing_files, collapse = ", "))
for (f in CODE_FILES) source(f)
stopifnot(exists("kudrin_draw"), exists("kudrin_tests"), exists("run_test_weight_fast"),
          exists("create_basis_stable"), exists("compute_epsilons"))
dir.create(file.path(OUT_DIR, "cells"), recursive = TRUE, showWarnings = FALSE)

if (CFG$run_selfchecks) {
  kudrin_selfcheck()          # DGP: vectorised t-statistics equal lm() for every specification
  kudrin_tests_selfcheck()    # Kudrin tests: Si, I(a), E[G_tilde], QP, full runs
}


## ---- 2. Kudrin's published curves -------------------------------------------------------
## Extracted from the vector paths of Fig. 2 in the November 14, 2025 version (accurate to the
## plotting precision, about 0.002). His tests are plotted at tau = 0, 0.1, ..., 1 and the EKW
## maximum at tau = 0, 0.05, ..., 1.

PUB <- list(
  tau     = seq(0, 1, by = 0.1),
  tau_ekw = seq(0, 1, by = 0.05),
  h1  = list(T_inf      = c(0.043, 0.055, 0.135, 0.219, 0.318, 0.527, 0.630, 0.850, 0.927, 0.972, 0.996),
             T_H        = c(0.052, 0.054, 0.093, 0.111, 0.157, 0.227, 0.306, 0.417, 0.532, 0.630, 0.718),
             T_H_hybrid = c(0.055, 0.055, 0.090, 0.100, 0.135, 0.380, 0.829, 0.993, 1.000, 1.000, 1.000),
             EKW_max    = c(0.084, 0.077, 0.068, 0.065, 0.075, 0.072, 0.069, 0.069, 0.067, 0.067, 0.059,
                            0.064, 0.067, 0.067, 0.068, 0.068, 0.069, 0.089, 0.086, 0.106, 0.122)),
  h2  = list(T_inf      = c(0.053, 0.038, 0.061, 0.042, 0.037, 0.043, 0.031, 0.037, 0.047, 0.046, 0.054),
             T_H        = c(0.065, 0.063, 0.066, 0.050, 0.071, 0.069, 0.066, 0.076, 0.085, 0.081, 0.096),
             T_H_hybrid = c(0.101, 0.080, 0.078, 0.055, 0.076, 0.072, 0.069, 0.082, 0.120, 0.195, 0.394),
             EKW_max    = c(0.068, 0.068, 0.063, 0.059, 0.067, 0.061, 0.063, 0.059, 0.069, 0.067, 0.070,
                            0.069, 0.069, 0.060, 0.067, 0.071, 0.066, 0.073, 0.067, 0.079, 0.096)),
  hPi = list(T_inf      = c(0.054, 0.057, 0.092, 0.143, 0.219, 0.342, 0.425, 0.572, 0.725, 0.863, 0.949),
             T_H        = c(0.056, 0.062, 0.077, 0.103, 0.128, 0.166, 0.200, 0.245, 0.318, 0.385, 0.452),
             T_H_hybrid = c(0.057, 0.063, 0.078, 0.102, 0.123, 0.164, 0.201, 0.248, 0.319, 0.388, 0.449),
             EKW_max    = c(0.050, 0.067, 0.069, 0.059, 0.064, 0.058, 0.053, 0.056, 0.056, 0.069, 0.057,
                            0.057, 0.057, 0.057, 0.057, 0.044, 0.051, 0.050, 0.047, 0.047, 0.044))
)

## For orientation only: a Python mirror of this design (paper settings for the projection
## test, warp-speed bootstrap, 1,000 repetitions) gave these projection-test rejection rates:
##   h = 1: .073 .040 .130 .188 .399 .622 .895 .983 1.00 1.00 1.00
##   h = 2: .064 .053 .052 .039 .047 .099 .062 .096 .098 .133 .201


## ---- 3. Pools of honest and p-hacked t-statistics (your kudrin_draw) --------------------

draw_effects <- function(panel, n) {
  if (panel == "h1") return(rep(1, n))
  if (panel == "h2") return(rep(2, n))
  if (panel == "hPi") {
    if (CFG$gamma_param == "scale") return(rgamma(n, shape = CFG$gamma_shape, scale = CFG$gamma_second))
    return(rgamma(n, shape = CFG$gamma_shape, rate = CFG$gamma_second))
  }
  stop("Unknown panel: ", panel)
}

pool_key <- function() {
  paste(CFG$n_pool, CFG$K, CFG$N, CFG$gamma_max, CFG$dgp, CFG$gamma_shape, CFG$gamma_second,
        CFG$gamma_param, CFG$alpha, CFG$seed, paste(CFG$panels, collapse = ","))
}

build_pools <- function(file) {
  if (file.exists(file)) {
    pools <- readRDS(file)
    if (identical(attr(pools, "key"), pool_key())) return(pools)
  }
  set.seed(CFG$seed)
  pools <- list()
  for (pn in CFG$panels) {
    message(sprintf("Simulating the %s pool (%s studies) ...", pn,
                    format(CFG$n_pool, big.mark = ",", scientific = FALSE)))
    secs <- system.time(
      d <- kudrin_draw(draw_effects(pn, CFG$n_pool), K = CFG$K, N = CFG$N,
                       gamma_max = CFG$gamma_max, alpha = CFG$alpha, dgp = CFG$dgp)
    )[["elapsed"]]
    pools[[pn]] <- abs(d[, c("honest", "minimum")])
    message(sprintf("  done in %.0f s", secs))
  }
  attr(pools, "key") <- pool_key()
  saveRDS(pools, file)
  pools
}

# Reference shares from an independent replica of EKW's Step 2/Step 3 design (2 million studies,
# Pi_hat with scale 1.2980). Differences above 0.005 are flagged.
check_pools <- function(pools) {
  ref <- rbind(h1  = c(0.1136, 0.2532, 0.2642, 0.4847),
               h2  = c(0.3040, 0.4970, 0.6080, 0.7958),
               hPi = c(0.2957, 0.4361, 0.4843, 0.6468))
  colnames(ref) <- c("honest p<=.05", "honest p<=.15", "minimum p<=.05", "minimum p<=.15")
  for (pn in names(pools)) {
    ph  <- 2 * pnorm(-pools[[pn]][, "honest"])
    pm  <- 2 * pnorm(-pools[[pn]][, "minimum"])
    got <- c(mean(ph <= 0.05), mean(ph <= 0.15), mean(pm <= 0.05), mean(pm <= 0.15))
    msg <- sprintf("Pool %-3s: %s", pn, paste(sprintf("%s %.4f", colnames(ref), got), collapse = ", "))
    comparable <- pn %in% rownames(ref) && !(pn == "hPi" && CFG$gamma_param != "scale")
    if (comparable) {
      dev <- max(abs(got - ref[pn, ]))
      msg <- paste0(msg, sprintf("  [max deviation from reference %.4f]", dev))
      if (dev > 0.005) warning("Pool ", pn, " differs from the reference by ", round(dev, 4))
    }
    message(msg)
  }
}


## ---- 4. Projection-test constants (your create_basis_stable / compute_epsilons) ----------

proj_constants <- function(file) {
  key <- paste(CFG$num_coeffs, CFG$L, CFG$numgrid,
               if (is.null(CFG$epsilon_U_pre)) "computed" else CFG$epsilon_U_pre)
  if (file.exists(file)) {
    pc <- readRDS(file)
    if (identical(pc$key, key)) return(pc)
  }
  eps <- CFG$epsilon_U_pre
  if (is.null(eps)) {
    message("Computing epsilon_U(L = ", CFG$L, ", numgrid = ", CFG$numgrid, ") once; this is slow ...")
    eps <- compute_epsilons(CFG$L, nx = CFG$numgrid)$epsilon_U
  }
  pc <- list(key = key, epsilon_U_pre = eps,
             U = create_basis_stable(CFG$num_coeffs, L = CFG$L, numgrid = CFG$numgrid, sigma_Y = 1))
  saveRDS(pc, file)
  pc
}


## ---- 5. One meta-sample ------------------------------------------------------------------

draw_meta <- function(pool, tau, S) {
  idx  <- sample.int(nrow(pool), S, replace = TRUE)
  hack <- runif(S) < tau
  ifelse(hack, pool[idx, "minimum"], pool[idx, "honest"])
}

hyb_name <- function(f) paste0("T_H_hybrid_floor", formatC(f, format = "e", digits = 0))

# EKW (2022) tests on (0, 0.15], called as in your run_sims(); solver failures count as
# non-rejections (EKW 2025, footnote 19). The robust caliper test is not in Tests.R.
ekw_rejections <- function(tstats) {
  ps   <- 2 * pnorm(-abs(tstats))
  pmn  <- min(ps[ps > 0])
  pmx  <- 0.15
  J    <- CFG$ekw_J
  safe <- function(expr) {
    v <- tryCatch(suppressWarnings(expr), error = function(e) NA_real_)
    if (length(v) != 1 || !is.numeric(v)) NA_real_ else as.numeric(v)
  }
  pv <- c(EKW_Binomial      = safe(Binomial(ps, 0.04, 0.05, "c")),
          EKW_Discontinuity = safe(Discontinuity_test(ps, 0.05)),
          EKW_Fisher        = safe(Fisher(ps, pmn, pmx)),
          EKW_LCM           = safe(LCM(ps, pmn, pmx)),
          EKW_CS1           = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 1, 0)),
          EKW_CSUB          = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 0, 1)),
          EKW_CS2B          = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 2, 1)))
  setNames(as.numeric(!is.na(pv) & pv <= CFG$alpha), names(pv))
}

# All tests on one meta-sample. Returns a named numeric vector (1 = reject).
one_rep <- function(r, panel, tau) {
  tstats <- draw_meta(POOLS[[panel]], tau, CFG$S)
  out    <- c(T_inf = NA_real_, T_H = NA_real_, T_H_hybrid = NA_real_, hyb_restricted = NA_real_)
  
  kt <- tryCatch(kudrin_tests(tstats, boots = CFG$kudrin_boots, alpha = CFG$alpha,
                              eig_floor = CFG$eig_floor[1]),
                 error = function(e) NULL)
  if (!is.null(kt)) {
    out[c("T_inf", "T_H", "T_H_hybrid")] <- as.numeric(kt$pval[c("T_inf", "T_H", "T_H_hybrid")] <= CFG$alpha)
    out[["hyb_restricted"]] <- as.numeric(kt$restricted)
  }
  for (f in CFG$eig_floor[-1]) {             # hybrid only; T_H and the hybrid use no bootstrap
    k2 <- tryCatch(kudrin_tests(tstats, boots = 1, alpha = CFG$alpha, eig_floor = f),
                   error = function(e) NULL)
    out[hyb_name(f)] <- if (is.null(k2)) NA_real_ else as.numeric(k2$pval[["T_H_hybrid"]] <= CFG$alpha)
  }
  
  if (CFG$run_projection) {
    pr <- tryCatch(run_test_weight_fast(data.frame(t = tstats, title = as.character(seq_along(tstats))),
                                        numcoeffs = CFG$num_coeffs, sigma_Y = 1, shift = CFG$shift,
                                        L = CFG$L, numgrid = CFG$numgrid,
                                        boots = if (CFG$proj_mode == "warp") 1 else CFG$proj_boots,
                                        U = PROJ$U, epsilon_U_pre = PROJ$epsilon_U_pre,
                                        p = CFG$p_weight, J0 = CFG$J0),
                   error = function(e) NULL)
    proj <- c(proj = NA_real_, proj_resid = NA_real_, proj_boot1 = NA_real_, proj_epsU = NA_real_)
    if (!is.null(pr)) {
      proj <- c(proj       = as.numeric(pr$pval <= CFG$alpha),
                proj_resid = as.numeric(pr$resid),
                proj_boot1 = as.numeric(pr$resids_boot[1]),
                proj_epsU  = as.numeric(pr$epsilon_U))
    }
    out <- c(out, proj)
  }
  
  if (CFG$run_ekw) out <- c(out, ekw_rejections(tstats))
  out
}


## ---- 6. Monte Carlo ----------------------------------------------------------------------

cell_key <- function() {
  paste(pool_key(), CFG$S, CFG$reps, CFG$kudrin_boots, paste(CFG$eig_floor, collapse = ","),
        CFG$run_projection, CFG$proj_mode, CFG$proj_boots, CFG$num_coeffs, CFG$J0, CFG$p_weight,
        CFG$shift, CFG$L, CFG$numgrid, CFG$run_ekw, CFG$ekw_J)
}

# PSOCK cluster, so it also runs on Windows; each worker sources your files.
start_cluster <- function() {
  if (CFG$n_cores <= 1) return(NULL)
  cl <- parallel::makePSOCKcluster(CFG$n_cores)
  parallel::clusterExport(cl, "CODE_FILES", envir = globalenv())
  invisible(parallel::clusterEvalQ(cl, { for (f in CODE_FILES) source(f); NULL }))
  parallel::clusterExport(cl, c("CFG", "POOLS", "PROJ", "draw_meta", "hyb_name", "ekw_rejections",
                                "one_rep"), envir = globalenv())
  parallel::clusterSetRNGStream(cl, CFG$seed)
  cl
}

run_cell <- function(panel, tau, cl) {
  file <- file.path(OUT_DIR, "cells", sprintf("cell_%s_tau%03d.rds", panel, round(100 * tau)))
  if (CFG$resume && file.exists(file)) {
    cell <- readRDS(file)
    if (identical(cell$key, cell_key())) return(cell)
  }
  t0 <- Sys.time()
  if (is.null(cl)) {
    L <- lapply(seq_len(CFG$reps), one_rep, panel = panel, tau = tau)
  } else {
    L <- parallel::parLapply(cl, seq_len(CFG$reps), one_rep, panel = panel, tau = tau)
  }
  nm   <- names(L[[1]])
  M    <- do.call(rbind, lapply(L, function(v) v[nm]))
  cell <- list(key = cell_key(), panel = panel, tau = tau, M = M,
               secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(cell, file)
  cell
}


## ---- 7. Rejection rates and comparison with the published curves ------------------------

summarise_cell <- function(cell) {
  M    <- cell$M
  rate <- function(col) mean(M[, col], na.rm = TRUE)
  s <- data.frame(panel = cell$panel, tau = cell$tau, reps = nrow(M),
                  T_inf = rate("T_inf"), T_H = rate("T_H"), T_H_hybrid = rate("T_H_hybrid"),
                  hybrid_restricted = rate("hyb_restricted"),
                  kudrin_failed = mean(is.na(M[, "T_inf"])), secs = cell$secs)
  for (f in CFG$eig_floor[-1]) s[[hyb_name(f)]] <- rate(hyb_name(f))
  if (CFG$run_projection) {
    if (CFG$proj_mode == "warp") {
      ok     <- is.finite(M[, "proj_resid"]) & is.finite(M[, "proj_boot1"])
      crit   <- quantile(M[ok, "proj_boot1"], 1 - CFG$alpha, names = FALSE)
      s$proj <- mean(M[ok, "proj_resid"] > crit + M[ok, "proj_epsU"])
    } else {
      s$proj <- rate("proj")
    }
    s$proj_failed <- mean(is.na(M[, "proj_resid"]))
  }
  if (CFG$run_ekw) {
    ek <- grep("^EKW_", colnames(M), value = TRUE)
    for (e in ek) s[[e]] <- rate(e)
    s$EKW_max <- max(unlist(s[ek]))
  }
  s
}

# Your rejection rates minus his, on the tau grid of his figure. The flag counts taus where the
# gap exceeds 3 standard errors (both simulations plus 0.002 of digitising error).
compare_with_published <- function(res) {
  out <- list()
  for (pn in intersect(unique(res$panel), names(PUB))) {
    r    <- res[res$panel == pn, ]
    keep <- abs(r$tau * 10 - round(r$tau * 10)) < 1e-8
    r    <- r[keep, ]
    if (nrow(r) == 0) next
    i <- round(r$tau * 10) + 1
    for (tn in c("T_inf", "T_H", "T_H_hybrid")) {
      ours <- r[[tn]]
      pub  <- PUB[[pn]][[tn]][i]
      se   <- sqrt(ours * (1 - ours) / r$reps + pub * (1 - pub) / 1000 + 0.002^2)
      out[[length(out) + 1]] <- data.frame(
        panel = pn, test = tn, n_tau = nrow(r),
        mean_diff = mean(ours - pub), max_abs_diff = max(abs(ours - pub)),
        tau_of_max = r$tau[which.max(abs(ours - pub))],
        n_beyond_3se = sum(abs(ours - pub) > 3 * se))
    }
  }
  do.call(rbind, out)
}


## ---- 8. Figure in Kudrin's format --------------------------------------------------------
## Colours, line styles and markers are taken from the PDF of his Fig. 2: T_inf dark red,
## solid, crosses; T_H yellow, dash-dot, asterisks; T_H^hybrid green, dash-dot, thick circles;
## EKW maximum blue, dashed; 5% level black, dotted; axes 0.15 grey with inward ticks on all
## four sides; panel boxes about 1.29:1. The projection test is added in purple.

STYLE <- list(
  sig        = list(col = "#000000", lty = "13", lwd = 0.8),
  EKW_max    = list(col = "#0000FF", lty = "44", lwd = 1.0),
  T_inf      = list(col = rgb(0.634766, 0.0778809, 0.184082), lty = "solid", lwd = 1.0,
                    pch = 4, cex = 0.75, plwd = 1.0),
  T_H        = list(col = rgb(0.929688, 0.693359, 0.125), lty = "C363", lwd = 1.0,
                    pch = 8, cex = 0.75, plwd = 1.0),
  T_H_hybrid = list(col = rgb(0, 0.501953, 0), lty = "C363", lwd = 1.0,
                    pch = 1, cex = 0.95, plwd = 1.6),
  proj       = list(col = rgb(0.494, 0.184, 0.556), lty = "solid", lwd = 1.4,
                    pch = 15, cex = 0.7, plwd = 1.0)
)
AXIS_COL <- grey(0.15)
PANEL_TITLE <- list(h1  = expression(italic(h) == 1*","~italic(K) == 3),
                    h2  = expression(italic(h) == 2*","~italic(K) == 3),
                    hPi = expression(italic(h) %~% hat(Pi)*","~italic(K) == 3))

lighten <- function(col, f = 0.65) {     # opaque tint (works on every device)
  m <- col2rgb(col) / 255
  rgb(t(m + (1 - m) * f))
}

draw_series <- function(x, y, st, faint = FALSE) {
  ok   <- is.finite(y)
  col  <- if (faint) lighten(st$col) else st$col
  lwd  <- if (faint) 0.8 * st$lwd else st$lwd
  lines(x[ok], y[ok], col = col, lty = st$lty, lwd = lwd)
  if (!is.null(st$pch)) {
    plwd <- if (faint) 0.8 * st$plwd else st$plwd
    points(x[ok], y[ok], col = col, pch = st$pch, cex = st$cex, lwd = plwd, xpd = TRUE)
  }
}

draw_panel <- function(res, pn, overlay_published, add_projection, legend_here) {
  r <- res[res$panel == pn, ]
  r <- r[order(r$tau), ]
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  at   <- c(0, 0.25, 0.5, 0.75, 1)
  labs <- c("0", "0.25", "0.50", "0.75", "1")
  for (side in 1:4) {
    axis(side, at = at, labels = if (side <= 2) labs else FALSE, lwd = 0, lwd.ticks = 0.5,
         col.ticks = AXIS_COL, col.axis = AXIS_COL, cex.axis = 0.82)
  }
  box(lwd = 0.5, col = AXIS_COL)
  mtext(expression("Extent of Selective Reporting ("*tau*")"), side = 1, line = 1.25,
        cex = 0.82, family = "serif")
  mtext("Power", side = 2, line = 2.0, cex = 0.82, family = "serif", las = 0)
  mtext(PANEL_TITLE[[pn]], side = 3, line = 0.2, cex = 0.82, family = "serif")
  
  abline(h = CFG$alpha, col = STYLE$sig$col, lty = STYLE$sig$lty, lwd = STYLE$sig$lwd)
  if ("EKW_max" %in% names(r)) {
    draw_series(r$tau, r$EKW_max, STYLE$EKW_max)
  } else {
    draw_series(PUB$tau_ekw, PUB[[pn]]$EKW_max, STYLE$EKW_max)
  }
  tests <- c("T_inf", "T_H", "T_H_hybrid")
  if (overlay_published) {
    for (tn in tests) draw_series(PUB$tau, PUB[[pn]][[tn]], STYLE[[tn]], faint = TRUE)
  }
  for (tn in tests) draw_series(r$tau, r[[tn]], STYLE[[tn]])
  show_proj <- add_projection && ("proj" %in% names(r))
  if (show_proj) draw_series(r$tau, r$proj, STYLE$proj)
  
  if (legend_here) {
    keys <- c("sig", "EKW_max", tests, if (show_proj) "proj")
    labs_leg <- c(expression("Significance level, 5%"), expression("Max Power of EKW Tests"),
                  expression(italic(T)[infinity]), expression(italic(T)[italic(H)]),
                  expression(italic(T)[italic(H)]^hybrid),
                  expression("Projection test"))[seq_along(keys)]
    sty <- STYLE[keys]
    legend("topleft", inset = c(0.02, 0.025), legend = labs_leg,
           col    = sapply(sty, function(s) s$col),
           lty    = sapply(sty, function(s) s$lty),
           lwd    = sapply(sty, function(s) s$lwd),
           pch    = sapply(sty, function(s) if (is.null(s$pch)) NA else s$pch),
           pt.cex = sapply(sty, function(s) if (is.null(s$cex)) 1 else s$cex),
           pt.lwd = sapply(sty, function(s) if (is.null(s$plwd)) 1 else s$plwd),
           cex = 0.62, bg = "white", box.lwd = 0.4, box.col = AXIS_COL, seg.len = 3,
           x.intersp = 0.6, y.intersp = 0.9)
  }
}

plot_fig2_minimum <- function(res, stem, overlay_published = FALSE, add_projection = TRUE,
                              legend = TRUE) {
  height <- if (overlay_published) 2.2 else 2.0
  for (ext in c("pdf", "png")) {
    file <- paste0(stem, ".", ext)
    if (ext == "pdf") {
      pdf(file, width = 6.6, height = height, pointsize = 8)
    } else if (isTRUE(capabilities("cairo"))) {
      png(file, width = 6.6, height = height, units = "in", res = 300, pointsize = 8, type = "cairo")
    } else {
      png(file, width = 6.6, height = height, units = "in", res = 300, pointsize = 8)
    }
    par(mfrow = c(1, 3))
    par(cex = 1, mar = c(2.5, 2.9, 1.1, 0.7), oma = c(if (overlay_published) 1.2 else 0, 0, 1.5, 0),
        mgp = c(1.4, 0.3, 0), tcl = 0.25, las = 1, family = "sans")
    for (k in seq_along(CFG$panels)) {
      draw_panel(res, CFG$panels[k], overlay_published, add_projection, legend_here = legend && k == 2)
    }
    mtext("Minimum", side = 3, outer = TRUE, line = 0.25, font = 2, cex = 1.3, family = "serif")
    if (overlay_published) {
      mtext("Faint: Kudrin (2025, Fig. 2) as published.   Solid: this replication.",
            side = 1, outer = TRUE, line = 0.2, cex = 0.8, family = "serif")
    }
    dev.off()
  }
  invisible(paste0(stem, c(".pdf", ".png")))
}


# His published curves as a results table (to redraw his panels with this code)
published_as_res <- function() {
  do.call(rbind, lapply(intersect(CFG$panels, names(PUB)), function(pn)
    data.frame(panel = pn, tau = PUB$tau, reps = 1000, T_inf = PUB[[pn]]$T_inf,
               T_H = PUB[[pn]]$T_H, T_H_hybrid = PUB[[pn]]$T_H_hybrid)))
}


## ---- 9. Run ------------------------------------------------------------------------------

# Redraw his published panels first, so the formatting can be compared with his PDF right away
plot_fig2_minimum(published_as_res(), file.path(OUT_DIR, "fig2_minimum_published_redrawn"),
                  add_projection = FALSE, legend = FALSE)

POOLS <- build_pools(file.path(OUT_DIR, "pools.rds"))
check_pools(POOLS)
PROJ <- if (CFG$run_projection) proj_constants(file.path(OUT_DIR, "proj_constants.rds")) else NULL
if (!is.null(PROJ)) {
  message(sprintf("Projection test: epsilon_U = %.5f per sqrt(n) (paper: 0.00063).", PROJ$epsilon_U_pre))
}

# Pilot: run time and the sup-norm bias bound
set.seed(CFG$seed + 1)
invisible(one_rep(1, CFG$panels[1], 1))       # warm-up: builds kudrin_tests' cached bounds
pilot <- system.time(one_rep(1, CFG$panels[1], 1))[["elapsed"]]
bbar  <- kudrin_tests(draw_meta(POOLS[[1]], 0, CFG$S), boots = 1)$bias_bound_inf
n_cells <- length(CFG$panels) * length(CFG$taus)
message(sprintf("Sup-norm bias bound %.3e (expected about 6.55e-04). One meta-sample takes %.1f s,", bbar, pilot))
message(sprintf("so %d cells x %d repetitions should take about %.1f hours on %d core(s).",
                n_cells, CFG$reps, pilot * CFG$reps * n_cells / 3600 / max(1, CFG$n_cores), CFG$n_cores))

cl  <- start_cluster()
res <- tryCatch({
  acc <- NULL
  for (pn in CFG$panels) {
    for (tau in CFG$taus) {
      s   <- summarise_cell(run_cell(pn, tau, cl))
      acc <- rbind(acc, s)
      message(sprintf("%-3s tau = %.1f:  T_inf %.3f   T_H %.3f   hybrid %.3f%s   (%.0f s)",
                      pn, tau, s$T_inf, s$T_H, s$T_H_hybrid,
                      if (CFG$run_projection) sprintf("   projection %.3f", s$proj) else "", s$secs))
    }
  }
  acc
}, finally = {
  if (!is.null(cl)) parallel::stopCluster(cl)
})

saveRDS(res, file.path(OUT_DIR, "results_fig2_minimum.rds"))
write.csv(res, file.path(OUT_DIR, "results_fig2_minimum.csv"), row.names = FALSE)

cmp <- compare_with_published(res)
write.csv(cmp, file.path(OUT_DIR, "comparison_with_published.csv"), row.names = FALSE)
message("Your Kudrin tests minus his published curves:")
print(cmp, digits = 3)
if (length(CFG$eig_floor) > 1) {
  message("Hybrid sensitivity to eig_floor (rejection rates):")
  print(res[, c("panel", "tau", "T_H_hybrid", sapply(CFG$eig_floor[-1], hyb_name))], digits = 3)
}

plot_fig2_minimum(res, file.path(OUT_DIR, "fig2_minimum_reproduction"))
plot_fig2_minimum(res, file.path(OUT_DIR, "fig2_minimum_vs_published"), overlay_published = TRUE)
message("Figures and tables written to ", OUT_DIR)