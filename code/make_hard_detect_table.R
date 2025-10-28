#Make table showing that the projection test can detect things that other tests can't
# Stefan Faridani
# October 14th, 2025

library(xtable)

root               <- "C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path
tables    <- paste0(root,"/tables") 
raw_simulations    <- paste0(root,"/simulations/raw simulations") #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking/simulations/raw simulations"
simulation_results <- paste0(root,"/simulations")


hard_detect_1 <- readRDS(paste0(simulation_results,"/to keep/main_sims1_nsims500_2025-10-12_04-55-40.rds"))
hard_detect_0 <- readRDS(paste0(simulation_results,"/to keep/aux_sims_1_nsims500_2025-10-12_00-14-40.rds"))

hard_detect <- rbind(hard_detect_1,hard_detect_0)

#The table
hard_detect <-round(hard_detect[,c("prob_hack","reject_rate_eps1", "CS2B_EWK","CS1_EWK", "rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK" )],2)

latex_table <- xtable(hard_detect,caption=this_script, label = "tab:hard_detect")
output_file <- paste0(tables,"/hard_detect.tex")
print(latex_table, file = output_file, include.rownames = FALSE)
cat("LaTeX table saved to", output_file, "\n")