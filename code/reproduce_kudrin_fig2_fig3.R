################################################################################
#  reproduce_kudrin_fig2_fig3.R
#
#  Reproduces all of Kudrin (2025, "Testing for and Evaluating the Extent of Selective
#  Reporting") Fig. 2 (K = 3, S = 1,000) and Fig. 3 (K = 3, S = 5,000): both rows
#  (Thresholding and Minimum), all three panels, and adds the projection test. It replaces
#    reproduce_Kudrin_figure_plus_proj.R            (Fig. 2, Minimum row)
#    reproduce_Kudrin_figure_threshold_plus_proj.R  (Fig. 2, Thresholding row)
#  with the same DGP, tests and settings in one run, and draws the figures in a new style.
#
#  What is simulated (Kudrin Sec. 5; EKW 2025 Secs. 5.1.1 and 5.3.1):
#    covariate selection with K = 3 controls, N = 200;
#    rows: thresholding (all K controls; if p > 0.05, the best of all (K-1)-subsets; ...; stop
#          at the first subset size whose best p-value is <= 0.05; if nothing is significant,
#          the smallest p-value overall) and minimum (smallest p-value over all 2^K
#          specifications);
#    panels h = 1, h = 2 and h ~ Pi_hat;  S = 1,000 (Fig. 2) and S = 5,000 (Fig. 3);
#    tau = 0, 0.1, ..., 1;  1,000 Monte Carlo repetitions per (S, row, panel, tau).
#    As in EKW/Kudrin, one million studies are simulated once per panel. One call of your
#    kudrin_draw() gives each study's honest (all controls), thresholded and minimum
#    t-statistic; with the seed of the two earlier scripts these are the same studies as theirs.
#    Each meta-sample draws S studies with replacement and reports, for each, the p-hacked
#    t-statistic with probability tau. Both rows and both figures draw from the same pools.
#
#  Common random numbers: every (S, panel, tau) cell has its own seed, shared by the two rows,
#  so at every tau both rows draw the same studies and the same hacking indicators; at tau = 0
#  they test identical honest meta-samples and agree exactly. His Fig. 2 has the same feature
#  for T_H, T_H^hybrid and the EKW maximum (section 2).
#
#  Uses your code only:
#    reproduce_kudrin_simulations_functions.R  kudrin_draw()            the DGP
#    kudrin_tests.R                            kudrin_tests()           T_inf, T_H, T_H^hybrid
#    functions_all_weighted.R                  run_test_weight_fast()   projection test
#                                              create_basis_stable(), compute_epsilons()
#    kudrin_figures_engine.R                   meta-samples, test calls, cells, figures (put it
#                                              in CODE_DIR with the files above)
#    Tests.R (EKW)                             only if CFG$run_ekw = TRUE
#
#  Outputs in OUT_DIR (N = 2 for Fig. 2 and 3 for Fig. 3):
#    figN_published_redrawn.pdf/.png  his published curves in the new style (drawn first)
#    figN_reproduction.pdf/.png       this replication, with the projection test
#    figN_vs_published.pdf/.png       this replication (solid) over his published curves (faint)
#    results_fig2_fig3.csv/.rds       rejection rates by S, row, panel, tau and test
#    comparison_with_published.csv    your curves minus his, with simulation-error flags
#  Fig. 2's cells run first, and its figures are written as soon as they finish. Finished cells
#  are saved in OUT_DIR/cells, so an interrupted run resumes where it stopped.
#
#  Run once with QUICK <- TRUE (minutes) to check that everything works, then with FALSE.
################################################################################

rm(list = ls())

## ---- 0. Settings -------------------------------------------------------------------------

ROOT     <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"

CODE_DIR <- file.path(ROOT, "code")
OUT_DIR  <- file.path(ROOT, "simulations", "kudrin_fig2_fig3")

QUICK <- FALSE   # TRUE: small run that checks the pipeline; FALSE: Kudrin's settings

CFG <- list(
  S_values       = c(1000, 5000),         # t-statistics per meta-sample: his Fig. 2 and Fig. 3
  taus           = seq(0, 1, by = 0.1),   # extent of selective reporting
  reps           = 1000,                  # Monte Carlo repetitions per (S, row, panel, tau)
  rows           = c("threshold", "minimum"),   # his top and bottom rows
  panels         = c("h1", "h2", "hPi"),  # h = 1, h = 2, h ~ Pi_hat
  alpha          = 0.05,                  # test level, and the significance threshold of the search

  # DGP, passed to kudrin_draw()
  n_pool         = 1e6,                   # studies per pool (EKW and Kudrin use one million)
  K              = 3,
  N              = 200,
  gamma_max      = 0.8,
  dgp            = "regression",
  # Pi_hat = Gamma(shape 1.4134, 1.2980). His footnote 13, and KUDRIN_GAMMA in your DGP file,
  # treat 1.2980 as a RATE (mean 1.09). His Fig. 2 is reproduced only when 1.2980 is a SCALE
  # (mean 1.83), presumably MATLAB's gamrnd(a, b) convention: in a Python mirror of the Minimum
  # row, T_inf at tau = 0.5 (h ~ Pi_hat) is 0.31 with the scale and 0.67 with the rate, against his 0.34.
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

  # EKW tests for the grey line (Tests.R)
  run_ekw        = FALSE,                 # FALSE: plot his published "Max Power of EKW Tests" line
  ekw_J          = 15,                    # bins on (0, 0.15], as in EKW (2025) Table 1
  # If TRUE, note that your DGP stops the search when the Student-t p-value (df = 196) is
  # <= 0.05, i.e. |t| >= 1.972. The reported p-values 2 Phi(-|t|) therefore bunch below 0.0486
  # and their density is flat across 0.05, which weakens the EKW tests anchored exactly at 0.05
  # (Discontinuity; Binomial on [0.04, 0.05]). The other EKW tests, Kudrin's tests (bins end at
  # 0.04 and 0.05, so the bunching stays inside one bin) and the projection test are
  # essentially unaffected.

  # Computing
  n_cores        = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  seed           = 20261022,              # the earlier scripts' seed, so the pools are the same
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

CODE_FILES  <- file.path(CODE_DIR, c(if (CFG$run_ekw) "Tests.R",
                                     "reproduce_kudrin_simulations_functions.R",
                                     "functions_all_weighted.R",
                                     "kudrin_tests.R"))
ENGINE_FILE <- file.path(CODE_DIR, "kudrin_figures_engine.R")
missing_files <- c(CODE_FILES, ENGINE_FILE)[!file.exists(c(CODE_FILES, ENGINE_FILE))]
if (length(missing_files) > 0) stop("Not found: ", paste(missing_files, collapse = ", "))
for (f in CODE_FILES) source(f)
source(ENGINE_FILE)
stopifnot(exists("kudrin_draw"), exists("kudrin_tests"), exists("run_test_weight_fast"),
          exists("create_basis_stable"), exists("compute_epsilons"), exists("plot_figure"))
dir.create(file.path(OUT_DIR, "cells"), recursive = TRUE, showWarnings = FALSE)

if (CFG$run_selfchecks) {
  kudrin_selfcheck()          # DGP: vectorised t-statistics equal lm() for every specification
  kudrin_tests_selfcheck()    # Kudrin tests: Si, I(a), E[G_tilde], QP, full runs
}


## ---- 2. Kudrin's published curves -------------------------------------------------------
## Figs. 2 and 3 of the November 14, 2025 version, read from the vector paths of the panel
## objects. Fig. 2 is identical to PUB in the two earlier scripts. His three tests are plotted
## at tau = 0, 0.1, ..., 1, and every value lies within 2e-6 of a multiple of 0.001: rejection
## frequencies out of 1,000 repetitions. The EKW maximum is plotted at tau = 0, 0.05, ..., 1 in
## Fig. 2 and at tau = 0, 0.1, ..., 1 in Fig. 3, where its values are multiples of 0.0002
## (5,000 repetitions). In the h ~ Pi_hat panel of Fig. 2 (Thresholding), T_H is hidden behind
## T_H^hybrid.
##
## At tau = 0 both of his rows test the same honest meta-samples. In Fig. 2, T_H, T_H^hybrid and
## the EKW maximum agree exactly across rows, but T_inf does not (0.036 vs 0.043, 0.026 vs
## 0.053, 0.045 vs 0.054), presumably because of his warp-speed bootstrap. In Fig. 3 only T_H
## agrees exactly; T_inf, T_H^hybrid (0.031 vs 0.020, 0.017 vs 0.021, 0.024 vs 0.025) and the
## EKW maximum do not. So his T_inf curves, and in Fig. 3 his hybrid curves, carry noise of their
## own that the binomial standard errors of compare_with_published() do not include: isolated
## flags for those tests at small tau are not by themselves a discrepancy.

PUB <- list(
  "1000" = list(
    figure  = 2,
    tau     = seq(0, 1, by = 0.1),
    tau_ekw = seq(0, 1, by = 0.05),
    reps    = 1000,
    rows    = list(
      threshold = list(
        h1  = list(T_inf      = c(0.036, 0.045, 0.127, 0.202, 0.348, 0.586, 0.760, 0.906, 0.953, 0.995, 1.000),
                   T_H        = c(0.052, 0.066, 0.121, 0.172, 0.292, 0.451, 0.593, 0.734, 0.857, 0.922, 0.961),
                   T_H_hybrid = c(0.055, 0.069, 0.127, 0.171, 0.311, 0.679, 0.975, 1.000, 1.000, 1.000, 1.000),
                   EKW_max    = c(0.084, 0.086, 0.080, 0.074, 0.081, 0.091, 0.110, 0.119, 0.132, 0.147, 0.169,
                              0.192, 0.227, 0.255, 0.279, 0.327, 0.384, 0.468, 0.547, 0.660, 0.777)),
        h2  = list(T_inf      = c(0.026, 0.049, 0.145, 0.419, 0.771, 0.921, 0.986, 0.997, 1.000, 1.000, 1.000),
                   T_H        = c(0.065, 0.093, 0.214, 0.427, 0.721, 0.897, 0.984, 0.998, 1.000, 1.000, 1.000),
                   T_H_hybrid = c(0.101, 0.143, 0.291, 0.603, 0.939, 0.996, 1.000, 1.000, 1.000, 1.000, 1.000),
                   EKW_max    = c(0.068, 0.068, 0.075, 0.088, 0.122, 0.156, 0.245, 0.378, 0.544, 0.723, 0.861,
                              0.950, 0.986, 0.997, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)),
        hPi = list(T_inf      = c(0.045, 0.047, 0.086, 0.223, 0.471, 0.698, 0.875, 0.968, 0.990, 0.999, 1.000),
                   T_H        = c(0.056, 0.075, 0.147, 0.259, 0.437, 0.656, 0.800, 0.917, 0.962, 0.994, 0.999),
                   T_H_hybrid = c(0.057, 0.075, 0.148, 0.263, 0.434, 0.657, 0.795, 0.921, 0.962, 0.986, 0.996),
                   EKW_max    = c(0.050, 0.069, 0.071, 0.073, 0.080, 0.085, 0.100, 0.120, 0.151, 0.180, 0.215,
                              0.256, 0.293, 0.349, 0.408, 0.484, 0.548, 0.611, 0.674, 0.737, 0.796))),
      minimum = list(
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
                              0.057, 0.057, 0.057, 0.057, 0.044, 0.051, 0.050, 0.047, 0.047, 0.044))))),
  "5000" = list(
    figure  = 3,
    tau     = seq(0, 1, by = 0.1),
    tau_ekw = seq(0, 1, by = 0.1),
    reps    = 1000,
    rows    = list(
      threshold = list(
        h1  = list(T_inf      = c(0.036, 0.109, 0.462, 0.849, 0.980, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
                   T_H        = c(0.023, 0.054, 0.336, 0.739, 0.963, 0.998, 1.000, 1.000, 1.000, 1.000, 1.000),
                   T_H_hybrid = c(0.031, 0.076, 0.351, 0.754, 0.967, 1.000, 0.999, 1.000, 1.000, 1.000, 1.000),
                   EKW_max    = c(0.0516, 0.0824, 0.1612, 0.2646, 0.4100, 0.6246, 0.8162, 0.9358, 0.9846, 0.9996, 1.0000)),
        h2  = list(T_inf      = c(0.024, 0.182, 0.818, 0.997, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
                   T_H        = c(0.022, 0.201, 0.856, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
                   T_H_hybrid = c(0.017, 0.221, 0.859, 0.998, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
                   EKW_max    = c(0.0546, 0.0838, 0.2878, 0.8580, 0.9996, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000)),
        hPi = list(T_inf      = c(0.028, 0.089, 0.441, 0.890, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
                   T_H        = c(0.019, 0.090, 0.497, 0.920, 0.993, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
                   T_H_hybrid = c(0.024, 0.110, 0.505, 0.910, 0.999, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
                   EKW_max    = c(0.0452, 0.0496, 0.1050, 0.2502, 0.5168, 0.7950, 0.9504, 0.9954, 0.9998, 1.0000, 1.0000))),
      minimum = list(
        h1  = list(T_inf      = c(0.043, 0.158, 0.444, 0.791, 0.945, 0.998, 1.000, 1.000, 1.000, 1.000, 1.000),
                   T_H        = c(0.023, 0.038, 0.115, 0.300, 0.551, 0.801, 0.935, 0.986, 0.998, 0.999, 1.000),
                   T_H_hybrid = c(0.020, 0.055, 0.134, 0.331, 0.565, 0.799, 0.944, 0.990, 0.998, 1.000, 1.000),
                   EKW_max    = c(0.0514, 0.0512, 0.0536, 0.0512, 0.0542, 0.0550, 0.0530, 0.0602, 0.1162, 0.2662, 0.5484)),
        h2  = list(T_inf      = c(0.027, 0.024, 0.032, 0.036, 0.036, 0.036, 0.042, 0.047, 0.070, 0.073, 0.073),
                   T_H        = c(0.022, 0.026, 0.021, 0.022, 0.028, 0.031, 0.044, 0.053, 0.053, 0.058, 0.070),
                   T_H_hybrid = c(0.021, 0.033, 0.024, 0.025, 0.037, 0.029, 0.039, 0.039, 0.046, 0.051, 0.089),
                   EKW_max    = c(0.0518, 0.0534, 0.0538, 0.0514, 0.0512, 0.0466, 0.0470, 0.0444, 0.0490, 0.0496, 0.0854)),
        hPi = list(T_inf      = c(0.038, 0.081, 0.176, 0.413, 0.781, 0.932, 0.997, 1.000, 1.000, 1.000, 1.000),
                   T_H        = c(0.019, 0.027, 0.056, 0.097, 0.232, 0.391, 0.586, 0.796, 0.918, 0.985, 0.996),
                   T_H_hybrid = c(0.025, 0.025, 0.055, 0.112, 0.231, 0.405, 0.615, 0.800, 0.929, 0.980, 0.999),
                   EKW_max    = c(0.0474, 0.0476, 0.0474, 0.0486, 0.0510, 0.0468, 0.0468, 0.0422, 0.0466, 0.0482, 0.0496)))))
)


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
  paste("kudrin_fig2_fig3", CFG$n_pool, CFG$K, CFG$N, CFG$gamma_max, CFG$dgp, CFG$gamma_shape,
        CFG$gamma_second, CFG$gamma_param, CFG$alpha, CFG$seed, paste(CFG$panels, collapse = ","))
}

# Same seed and the same kudrin_draw() calls, in the same order, as the two earlier scripts, so
# the same studies; they kept two of these three columns each.
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
    pools[[pn]] <- abs(d[, c("honest", "threshold", "minimum")])
    message(sprintf("  done in %.0f s", secs))
  }
  attr(pools, "key") <- pool_key()
  saveRDS(pools, file)
  pools
}

# Reference shares of 2 Phi(-|t|) <= cutoff, from independent replicas of this design: the
# Thresholding script's (10 million studies per panel; Student-t p-values in the search, as in
# your DGP; Pi_hat with scale 1.2980) for the honest and thresholded columns, and the Minimum
# script's (2 million studies) for the minimum column. Thresholding and minimum searches have the
# SAME shares at every cutoff >= 0.05 (both report a significant specification whenever one
# exists, and the same overall minimum otherwise), so the thresholded column is also checked at
# p <= .01 and .04, where it has far fewer small p-values. Differences above 0.005 are flagged.
check_pools <- function(pools) {
  cuts <- c(honest = 0.05, honest = 0.15, threshold = 0.01, threshold = 0.04, threshold = 0.05,
            minimum = 0.05, minimum = 0.15)
  ref <- rbind(h1  = c(0.1138, 0.2533, 0.0498, 0.2061, 0.2643, 0.2642, 0.4847),
               h2  = c(0.3040, 0.4973, 0.1917, 0.5183, 0.6083, 0.6080, 0.7958),
               hPi = c(0.2962, 0.4370, 0.2197, 0.4268, 0.4850, 0.4843, 0.6468))
  colnames(ref) <- c("honest p<=.05", "honest p<=.15", "threshold p<=.01", "threshold p<=.04",
                     "threshold p<=.05", "minimum p<=.05", "minimum p<=.15")
  for (pn in names(pools)) {
    got <- sapply(seq_along(cuts), function(j)
      mean(2 * pnorm(-pools[[pn]][, names(cuts)[j]]) <= cuts[[j]]))
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


## ---- 4. Rows, panels and cells ------------------------------------------------------------

ROWS <- list(threshold = list(id = "threshold", title = "Thresholding", column = "threshold"),
             minimum   = list(id = "minimum",   title = "Minimum",      column = "minimum"))[CFG$rows]
PANELS <- list(
  h1  = list(id = "h1",  title = as.expression(bquote(italic(h) == 1*","~italic(K) == .(CFG$K)))),
  h2  = list(id = "h2",  title = as.expression(bquote(italic(h) == 2*","~italic(K) == .(CFG$K)))),
  hPi = list(id = "hPi", title = as.expression(bquote(italic(h) %~% hat(Pi)*","~italic(K) == .(CFG$K))))
)[CFG$panels]

# Both rows of a panel draw from that panel's pool
CELLS <- make_cells(ROWS, PANELS, pool_of = function(rw, pn) pn$id)


## ---- 5. Run ------------------------------------------------------------------------------

# His published curves in the new style first, so the style can be checked right away
for (S in CFG$S_values) {
  P <- PUB[[as.character(S)]]
  plot_figure(published_as_res(P, S), S, ROWS, PANELS,
              file.path(OUT_DIR, sprintf("fig%d_published_redrawn", P$figure)), pub = P,
              add_projection = FALSE)
}

POOLS <- build_pools(file.path(OUT_DIR, "pools.rds"))
check_pools(POOLS)
POOL_KEYS <- setNames(as.list(paste(pool_key(), names(POOLS))), names(POOLS))
PROJ <- if (CFG$run_projection) proj_constants(file.path(OUT_DIR, "proj_constants.rds")) else NULL
if (!is.null(PROJ)) {
  message(sprintf("Projection test: epsilon_U = %.5f per sqrt(n) (paper: 0.00063).", PROJ$epsilon_U_pre))
}

pilot_estimate(CELLS)

cl  <- start_cluster()
res <- tryCatch({
  acc <- NULL
  for (S in CFG$S_values) {
    P    <- PUB[[as.character(S)]]
    acc  <- rbind(acc, run_cells(Filter(function(cc) cc$S == S, CELLS), cl))
    saveRDS(acc, file.path(OUT_DIR, "results_fig2_fig3.rds"))
    stem <- file.path(OUT_DIR, sprintf("fig%d", P$figure))
    plot_figure(acc, S, ROWS, PANELS, paste0(stem, "_reproduction"), pub = P)
    plot_figure(acc, S, ROWS, PANELS, paste0(stem, "_vs_published"), pub = P, overlay_published = TRUE,
                note = sprintf("Faint: Kudrin (2025, Fig. %d) as published.   Solid: this replication.",
                               P$figure))
    message(sprintf("Fig. %d (S = %d) written to %s", P$figure, as.integer(S), OUT_DIR))
  }
  acc
}, finally = {
  if (!is.null(cl)) parallel::stopCluster(cl)
})

saveRDS(res, file.path(OUT_DIR, "results_fig2_fig3.rds"))
write.csv(res, file.path(OUT_DIR, "results_fig2_fig3.csv"), row.names = FALSE)

cmp <- do.call(rbind, lapply(CFG$S_values, function(S) compare_with_published(res, S, PUB[[as.character(S)]])))
write.csv(cmp, file.path(OUT_DIR, "comparison_with_published.csv"), row.names = FALSE)
message("Your Kudrin tests minus his published curves:")
print(cmp, digits = 3)
tau0_check(res)
if (length(CFG$eig_floor) > 1) {
  message("Hybrid sensitivity to eig_floor (rejection rates):")
  print(res[, c("S", "row", "panel", "tau", "T_H_hybrid", sapply(CFG$eig_floor[-1], hyb_name))], digits = 3)
}
message("Figures and tables written to ", OUT_DIR)
