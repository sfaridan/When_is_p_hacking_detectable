rm(list=ls())
gc()

#install.packages("osqp")                 # CRAN package for OSQP interface :contentReference[oaicite:2]{index=2}
#install.packages("Matrix")               # Sparse‐matrix support :contentReference[oaicite:3]{index=3}
library(osqp)
library(Matrix)
library(tictoc)

root               <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path
simulation_results <- paste0(root,"/simulations")

#source("C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking/Tests.R")
source(paste0(root,"/code/Tests.R"))
source(paste0(root,"/code/functions_all.R"))


### Choose parameters for the simulations
parms <- expand.grid(nsims       = 200,            # number of sim repetitions
                     nboots      = 100,           # number of bootstraps per sim
                     n           = c(5000),           # meta-sample size
                     cv          = 1.96,           # critical value to shift by
                     sigma_Y     = 1,              # set this to one
                     prob_hack   = c(1),              #probability to take the larger of two t-scores
                     num_coeffs  = c(31),       # the larger the less regularized
                     nu          = c(99999),    # dof for true dgp
                     theta       = c(1),       # probability of reporting t when |t|<cv
                     numgrid     = 3000,           # number of hn grid to make U
                     L           = 6.5,            # width of grid of hs to make U  
                     pi0_shape   = c("null","double_normal","poisson","chi2"),       # shape of distribution of h
                     h_center    = c(2),              # center of true effect distribution
                     sigma_h     = 1,              # increases sd of h
                     bimodal     = c(FALSE),          # separate distribution of h by 2
                     hack_type   = "threshold",         # maximization hacking (threshold is default)  
                     omit_proj   = FALSE,
                     omit_EWK    = FALSE,
                     shift       =1.96,
                     seed        =2 
)
### Run the simulations 
setwd(simulation_results)
out_file<- run_sims(parms,"sims_power_newsolver_LAPTOP")
out_parms <- readRDS(out_file)
print(out_parms)



