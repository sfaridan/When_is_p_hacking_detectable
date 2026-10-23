rm(list=ls())
gc()

library(Matrix)
library(tictoc)

root               <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
this_script        <- rstudioapi::getSourceEditorContext()$path
simulation_results <- paste0(root,"/simulations")

source(paste0(root,"/code/Tests.R")) #EKW tests - their code
source(paste0(root,"/code/reproduce_kudrin_simulations_functions.R")) #reverse-engineer Kudrin (2024)'s simulations
source(paste0(root,"/code/functions_all_weighted.R")) 
source(paste0(root,"/code/kudrin_tests.R")) 


### Choose parameters for the simulations
parms <- expand.grid(nsims       = 500,            # number of sim repetitions
                     nboots      = 100,           # number of bootstraps per sim
                     n           = c(10000),           # meta-sample size
                     cv          = 1.96,           # critical value to shift by
                     sigma_Y     = 1,              # set this to one
                     prob_hack   = c(0,1),              #probability to take the larger of two t-scores
                     num_coeffs  = c(41),       # the larger the less regularized
                     nu          = c(99999),    # dof for true dgp
                     theta       = c(1),       # probability of reporting t when |t|<cv
                     numgrid     = 3000,           # number of hn grid to make U
                     L           = 6.5,            # width of grid of hs to make U  
                     pi0_shape   = c("hard"),       # shape of distribution of h
                     h_center    = c(2),              # center of true effect distribution
                     sigma_h     = 0.7,              # increases sd of h
                     hack_type   = "max",         # maximization hacking (threshold is default)  
                     omit_proj   = FALSE,
                     omit_EWK    = FALSE,
                     kudrin_tests= TRUE,
                     shift       = 1.96,
                     J0          = 20,
                     p           = 2,
                     seed        = 3 
)
### Run the simulations 
setwd(simulation_results)
out_file<- run_sims(parms,"sims_weighted_hard_comparison")
out_parms <- readRDS(out_file)
print(out_parms)


