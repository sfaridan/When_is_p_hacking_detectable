#### =====================================================================================================
#### run_cochrane_and_mm_scripts.R: every script of the two applications, in order.
#### Each script starts with rm(list = ls()), which would delete this driver's own variables, so the driver keeps
#### them inside local(), where the scripts cannot reach them. The scripts themselves are sourced into the global
#### environment (source()'s default), which cochrane_power_experiment.R needs: its parallel workers get their
#### objects from there. The run stops at the first error; the scripts that ran before it keep their output.
#### =====================================================================================================
rm(list = ls())
gc(reset = TRUE)
options(warn = 1)                      #show each warning when it happens, not after the whole run

local({
  code_dir <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable/code"
  scripts  <- c("application_cochrane.R",       #Cochrane: the sample, its report, every test
                "application_mm_analysis.R",    #Brodeur: the sample, its report, every test
                "cochrane_exhibits.R",          #Cochrane: appendix table and T_inf gap figure
                "application_mm_output.R",      #Brodeur: appendix tables and figures
                "main_applications_table.R",    #the main-text table, from both applications
                "cochrane_power_experiment.R")  #last: the longest, and nothing above needs it
  missing <- scripts[!file.exists(file.path(code_dir, scripts))]
  if (length(missing)) stop("Not in ", code_dir, ": ", paste(missing, collapse = ", "))
  for (f in scripts) {
    t0 <- Sys.time()
    cat("\n================ ", f, ", started ", format(t0, "%Y-%m-%d %H:%M"), " ================\n", sep = "")
    source(file.path(code_dir, f))
    graphics.off()                     #close any graphics device a script left open
    cat(sprintf("================ %s done in %.1f minutes\n", f, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
})