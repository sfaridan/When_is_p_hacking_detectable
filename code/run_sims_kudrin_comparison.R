rm(list=ls())
gc()

library(Matrix)
library(tictoc)

root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
this_script        <- rstudioapi::getSourceEditorContext()$path
simulation_results <- paste0(root,"/simulations")

#source(paste0(root,"/code/Tests.R"))
source(paste0(root,"/code/functions_all_weighted.R"))   # patched run_sims() with the Kudrin hooks
source(paste0(root,"/code/reproduce_kudrin_simulations_functions.R")) 
source(paste0(root,"/code/kudrin_tests.R")) 


#### Stefan's parameterization
parms <- expand.grid(nsims       = 40, #500,
                     nboots      = 100,
                     n           = 10000, #c(1000, 5000),
                     cv          = 1.96,
                     sigma_Y     = 1,
                     prob_hack   = c(1), #c(0, 0.25, 0.5, 0.75, 1),
                     num_coeffs  = 41,
                     nu          = 99999,
                     theta       = 1,
                     numgrid     = 3000,
                     L           = 6.5,
                     pi0_shape   = c("point"), #c("point", "kudrin_gamma"),   # h = 2, or h ~ Gamma(1.4134, rate 1.2980)
                     h_center    = c(2),                            # ignored by "kudrin_gamma"
                     sigma_h     = 0.7,
                     hack_type   = "kudrin_min",
                     kudrin_K    = 3, #c(3, 5, 7),
                     kudrin_N    = 200,
                     kudrin_dgp  = "regression",
                     omit_proj   = FALSE,
                     omit_EWK    = TRUE,
                     kudrin_tests= TRUE,
                     shift       = 1.96,
                     J0          = 30,
                     p           = 2,
                     seed        = 3,
                     stringsAsFactors = FALSE)
parms$nu <- parms$kudrin_N - parms$kudrin_K - 1
setwd(simulation_results)
out_file  <- run_sims(parms, "sims_weighted_kudrin_Fig3_bottom_withTmax")
out_parms <- readRDS(out_file)
print(out_parms[, c("kudrin_K", "pi0_shape", "h_center", "n", "prob_hack", "reject_rate","Kudrin_T_inf", "Kudrin_T_H", "Kudrin_T_H_hybrid")])