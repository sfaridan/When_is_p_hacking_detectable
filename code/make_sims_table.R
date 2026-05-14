###### Creates table of simulations for threshold p-hacking ######


rm(list = ls())
gc()

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
kept_sims <- paste0(root,"/simulations/to keep/two rounds size power n 5k")
setwd(kept_sims)

#average together the n=200 and n=300 size sims
size_round1 <- readRDS("sims_size_newsolver_LAPTOP4_nsims200_2026-05-08_01-05-03.rds")
size_round2 <- readRDS("sims_size_secondround_LAPTOP4_nsims300_2026-05-09_06-56-17.rds")
vars_num <- c("prob_hack", "reject_rate","rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK", "CS1_EWK", "CS2B_EWK")
vars_all <- c("prob_hack","pi0_shape", "reject_rate","rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK", "CS1_EWK", "CS2B_EWK")
n1 <- size_round1$nsims[1]
n2 <- size_round2$nsims[2]
table_size1 <-  size_round1[,vars_num]*(n1/(n1+n2))+ size_round2[,vars_num]*(n2/(n1+n2))
table_size1$pi0_shape <-  size_round1$pi0_shape 
table_size1 <- table_size1[,vars_all]

#average together the n=200 and n=300 power sims
power_round1 <- readRDS("sims_power_newsolver_LAPTOP4_nsims200_2026-05-08_01-24-54.rds")
power_round2 <- readRDS("sims_power_secondround_LAPTOP4_nsims300_2026-05-09_08-41-46.rds")
n1 <- power_round1$nsims[1]
n2 <- power_round2$nsims[2]
table_power1 <-  power_round1[,vars_num]*(n1/(n1+n2))+ power_round2[,vars_num]*(n2/(n1+n2))
table_power1$pi0_shape <-  power_round1$pi0_shape 
table_power1 <- table_power1[,vars_all]

# add the h=2 size simulation
size_h2 <- readRDS("sims_size_moredgps_LAPTOP3_nsims500_2026-05-09_18-07-18.rds")
table_size_h2 <- size_h2[1,vars_all]
table_size_h2$pi0_shape <- "H=0.5"

# add the h=2 power simulation
power_h2 <- readRDS("sims_power_moredgps_LAPTOP3_nsims500_2026-05-09_23-29-09.rds")
table_power_h2 <- power_h2[1,vars_all]
table_power_h2$pi0_shape <- "H=0.5"

#add nearzero size simulation
size_unif <- readRDS("sims_size_nearzero_LAPTOP2_nsims500_2026-05-10_15-09-40.rds")
table_size_unif <- size_unif[,vars_all]

#add nearzero size simulation
power_unif <- readRDS("sims_power_nearzero_LAPTOP2_nsims500_2026-05-11_01-44-10.rds")
table_power_unif <- power_unif[,vars_all]

table_all <- rbind(table_size1,table_size_h2,table_size_unif, table_power1,table_power_h2,table_power_unif)

#adjust names
table_all$pi0_shape <- as.character(table_all$pi0_shape)
table_all$pi0_shape[table_all$pi0_shape == "null"] <- "H=0"
table_all$pi0_shape[table_all$pi0_shape == "uniform"] <- "Unif(-0.5,0.5)"
table_all$pi0_shape[table_all$pi0_shape == "double_normal"] <- "Mix Normals"
table_all$pi0_shape[table_all$pi0_shape == "normal"] <- "Norm(0,1)"
table_all$pi0_shape[table_all$pi0_shape == "chi2"] <- "chi^2(2)"
table_all$pi0_shape[table_all$pi0_shape == "poisson"] <- "Poisson(2)"
table_all


library(knitr)

table_out <- table_all
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