#Make table Comparing power of the projection test across methods
# Stefan Faridani
# October 14th, 2025

library(xtable)

root               <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
this_script        <- rstudioapi::getSourceEditorContext()$path
tables    <- paste0(root,"/tables") 
simulation_results <- paste0(root,"/simulations")


comparison_1 <- readRDS(paste0(simulation_results,"/to keep/power_comparison_aux3_nsims500_2025-10-27_11-17-27.rds"))
comparison_0 <- readRDS(paste0(simulation_results,"/to keep/power_comparison_sims3_nsims500_2025-10-27_12-26-28.rds"))

comparison <- rbind(comparison_1,comparison_0)
comparison <- comparison[order(comparison$prob_hack), ]

#The table
comparison <-round(comparison[,c("prob_hack","reject_rate", "CS2B_EWK","CS1_EWK", "rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK" )],2)

latex_table <- xtable(comparison,caption=this_script, label = "tab:comparison")
output_file <- paste0(tables,"/comparison.tex")
print(latex_table, file = output_file, include.rownames = FALSE)
cat("LaTeX table saved to", output_file, "/n")