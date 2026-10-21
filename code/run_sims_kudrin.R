rm(list=ls())
gc()

library(Matrix)
library(tictoc)

root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
this_script        <- rstudioapi::getSourceEditorContext()$path
simulation_results <- paste0(root,"/simulations")

source(paste0(root,"/code/Tests.R"))
source(paste0(root,"/code/functions_all_weighted.R"))   # patched run_sims() with the Kudrin hooks
source(paste0(root,"/code/reproduce_kudrin_simulations_functions.R")) 

### One-time checks of the DGP code
kudrin_selfcheck(K = 3); kudrin_selfcheck(K = 7)   # t-stats == lm(), and the run_sims() hook runs
# kudrin_fig7_check()                              # compare with Kudrin (2025), Fig. 7

### Kudrin (2025), Fig. 3, bottom-middle panel: minimum approach, h = 2, K = 3, S = 5000.
### Every test there (T_inf, T_H, T_H^hybrid, max power of EKW tests) has rejection rate <= 0.08 for all tau.
parms <- expand.grid(nsims       = 500,            # number of sim repetitions
                     nboots      = 100,            # number of bootstraps per sim
                     n           = 5000,           # meta-sample size (Kudrin's S)
                     cv          = 1.96,           # critical value to shift by
                     sigma_Y     = 1,              # set this to one
                     prob_hack   = seq(0, 1, by = 0.1),   # Kudrin's tau: fraction of studies that p-hack
                     num_coeffs  = 41,             # the larger the less regularized
                     nu          = 99999,          # overwritten below (only feeds the nu_resid bookkeeping)
                     theta       = 1,              # no publication bias (as in Kudrin)
                     numgrid     = 3000,           # number of hn grid to make U
                     L           = 6.5,            # width of grid of hs to make U
                     pi0_shape   = "point",        # h_s = h_center for every study
                     h_center    = 2,
                     sigma_h     = 1,
                     hack_type   = "kudrin_min",   # "kudrin_thresh" = general-to-specific thresholding
                     kudrin_K    = 3,              # number of candidate controls
                     kudrin_N    = 200,            # observations per study
                     kudrin_dgp  = "regression",   # "asymptotic" = limit experiment with exactly normal t-stats
                     omit_proj   = FALSE,
                     omit_EWK    = TRUE,           # FALSE adds the EKW tests (Kudrin's dashed line)
                     shift       = 1.96,
                     J0          = 30,
                     p           = 2,              # w_j = min{1, 0.5^(j - 30)}
                     seed        = 3,
                     stringsAsFactors = FALSE)
# Residual df of the honest (all-controls) regression, so nu_resid / mean_cvs_full report the
# Student-t Delta of the paper's Section 5.2 for this design:
parms$nu <- parms$kudrin_N - parms$kudrin_K - 1
setwd(simulation_results)
out_file  <- run_sims(parms, "sims_weighted_kudrin_fig3")
out_parms <- readRDS(out_file)
print(out_parms[, c("kudrin_K", "pi0_shape", "h_center", "n", "prob_hack", "reject_rate")])

### The other minimum-approach panels where Kudrin's tests have little or no power:
### h = 2 (all K, S in {1000, 5000}; Figs. 2, 3, 8, 9) and h ~ Gamma for K in {5, 7} (Figs. 8, 9).
### The full grid below is 60 parameterizations (~20 hours at ~2.5 s per sim); trim as needed.
parms <- expand.grid(nsims       = 50, #500,
                     nboots      = 100,
                     n           = 5000, #c(1000, 5000),
                     cv          = 1.96,
                     sigma_Y     = 1,
                     prob_hack   = c(1,0), #c(0, 0.25, 0.5, 0.75, 1),
                     num_coeffs  = 41,
                     nu          = 99999,
                     theta       = 1,
                     numgrid     = 3000,
                     L           = 6.5,
                     pi0_shape   = c("point", "kudrin_gamma"),   # h = 2, or h ~ Gamma(1.4134, rate 1.2980)
                     h_center    = 2,                            # ignored by "kudrin_gamma"
                     sigma_h     = 1,
                     hack_type   = "kudrin_min",
                     kudrin_K    = 3, #c(3, 5, 7),
                     kudrin_N    = 200,
                     kudrin_dgp  = "regression",
                     omit_proj   = FALSE,
                     omit_EWK    = TRUE,
                     shift       = 1.96,
                     J0          = 30,
                     p           = 2,
                     seed        = 3,
                     stringsAsFactors = FALSE)
parms$nu <- parms$kudrin_N - parms$kudrin_K - 1
setwd(simulation_results)
out_file  <- run_sims(parms, "sims_weighted_kudrin_lowpower")
out_parms <- readRDS(out_file)
print(out_parms[, c("kudrin_K", "pi0_shape", "h_center", "n", "prob_hack", "reject_rate")])