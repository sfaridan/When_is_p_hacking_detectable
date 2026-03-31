#install.packages("osqp")                 # CRAN package for OSQP interface :contentReference[oaicite:2]{index=2}
#install.packages("Matrix")               # Sparse‐matrix support :contentReference[oaicite:3]{index=3}
library(osqp)
library(Matrix)
library(tictoc)

root               <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path
raw_simulations    <- paste0(root,"/simulations") #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking/simulations/raw simulations"
simulation_results <- paste0(root,"/simulations")

#source("C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking/Tests.R")
source(paste0(root,"/code/Tests.R"))
source(paste0(root,"/code/functions_sims.R"))


### Choose parameters for the simulations
parms <- expand.grid(nsims       =   500,            # number of sim repetitions
                     nboots      =   100,           # number of bootstraps per sim
                     n           = c(100000),           # meta-sample size
                     cv          = 1.96,           # critical value to shift by
                     sigma_Y     = 1,              # set this to one
                     prob_hack   = c(1),              #probability to take the larger of two t-scores
                     num_coeffs  = c(30),       # the larger the less regularized
                     nu          = c(99999),    # dof for true dgp
                     theta       = c(1),       # probability of reporting t when |t|<cv
                     numgrid     =  3000,           # number of hn grid to make U
                     L           = 6.5,            # width of grid of hs to make U  
                     h           = 1.96,              # expectation of true effect distribution
                     sigma_h     = 0.7,            # standard deviation of true effect distribution
                     expo        = FALSE,
                     smooth_hack = TRUE,
                     omit_proj   = FALSE,
                     omit_EWK    = FALSE,
                     shift       =1,
                     seed        =1,
                     fastboot    =c(FALSE)         # compute cvs once for each parameterization to speed up the simulation
)
### Run the simulations 
setwd(simulation_results)
out_file<- run_sims(parms,"sims_size_hard_to_detect")
out_parms <- readRDS(out_file)
print(out_parms)



