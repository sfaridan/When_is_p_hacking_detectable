###### Creates table of simulations for threshold p-hacking ######


rm(list = ls())
gc()

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
kept_sims <- paste0(root,"/simulations/to keep/two rounds size power n 5k")
setwd(kept_sims)

#average together the n=200 and n=300 size sims
sims_hard <- readRDS("sims_hard_laptop2_nsims500_2026-05-11_20-15-59.rds")

table_hard <- sims_hard[,c("prob_hack","pi0_shape","reject_rate", "rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK",    "CS1_EWK",   "CS2B_EWK")]

library(knitr)

table_out <- table_hard
num_cols <- sapply(table_out, is.numeric)
table_out[num_cols] <- lapply(table_out[num_cols], function(x) round(x, 3))

latex_tab <- kable(
  table_out,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  caption = "Simulation results"
)

writeLines(latex_tab, file.path(paste0(root,"/tables/"), "simulation_results.tex"))
latex_tab



######