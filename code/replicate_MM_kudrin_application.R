# 1. Clear all objects from the global environment
rm(list = ls())

# 2. Force R to release the freed memory to your OS
gc(reset = TRUE)

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_all_weighted.R"))

#Load in the data
data_location <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data"
setwd(data_location)
MM_raw <- haven::read_dta('MM data fromkudrin.dta')

source(paste0(root, "/code/kudrin_cleaning.R"))
MM_dropped <- prep_kudrin("MM data fromkudrin.dta")   # raw t on his sample: for the replication block
set.seed(1)
#MM_data    <- deround_kudrin(MM_dropped)              # one de-rounding draw: for your own analysis

boots <- 10000

## get as close as we can to Kudrin 2025 ##

source(paste0(root,"/code/kudrin_tests.R"))
#kudrin_tests_selfcheck() #optional one-time check of kudrin_tests.R

run_kudrin <- function(data){
  kt <- kudrin_tests(data$t, cluster = data$title, boots = boots) #cluster by article, as in the projection test
  data.frame(n = length(data$t), articles = length(unique(data$title)), as.list(kt$pval),
             min = min(kt$pval), restricted = kt$restricted) #restricted = FALSE means T_H_hybrid = T_H
}
run_kudrin_nocluster <- function(data){
  kt <- kudrin_tests(data$t, cluster = NULL, boots = boots) #cluster by article, as in the projection test
  data.frame(n = length(data$t), articles = length(unique(data$title)), as.list(kt$pval),
             min = min(kt$pval), restricted = kt$restricted) #restricted = FALSE means T_H_hybrid = T_H
}

set.seed(1) #T_inf is bootstrapped; reset the seed so these results don't depend on the code above
kudrin_rct_noderound <- run_kudrin(MM_dropped[MM_dropped$method=="RCT",])
kudrin_iv_noderound   <- run_kudrin(MM_dropped[MM_dropped$method=="IV",])
kudrin_did_noderound  <- run_kudrin(MM_dropped[MM_dropped$method=="DID",])
kudrin_rdd_noderound  <- run_kudrin(MM_dropped[MM_dropped$method=="RDD",])

set.seed(1) #T_inf is bootstrapped; reset the seed so these results don't depend on the code above
kudrin_rct_noderound_nocluster <- run_kudrin_nocluster(MM_dropped[MM_dropped$method=="RCT",])
kudrin_iv_noderound_nocluster   <- run_kudrin_nocluster(MM_dropped[MM_dropped$method=="IV",])
kudrin_did_noderound_nocluster  <- run_kudrin_nocluster(MM_dropped[MM_dropped$method=="DID",])
kudrin_rdd_noderound_nocluster  <- run_kudrin_nocluster(MM_dropped[MM_dropped$method=="RDD",])

empirical_kudrin_noderound_nocluster  <- cbind(methods = c("RCT","IV","DID","RDD"), round(rbind(kudrin_rct_noderound_nocluster,kudrin_iv_noderound_nocluster,kudrin_did_noderound_nocluster,kudrin_rdd_noderound_nocluster),4))
empirical_kudrin_noderound  <- cbind(methods = c("RCT","IV","DID","RDD"), round(rbind(kudrin_rct_noderound,kudrin_iv_noderound,kudrin_did_noderound,kudrin_rdd_noderound),4))

empirical_kudrin_replicate <- empirical_kudrin_noderound
empirical_kudrin_replicate$T_inf <- empirical_kudrin_noderound_nocluster$T_inf
empirical_kudrin_replicate$min <- with(empirical_kudrin_replicate, pmin(T_inf, T_H, T_H_hybrid))

setwd(paste0(root,"/results"))
saveRDS(empirical_kudrin_replicate, file = "empirical_kudrin_replicate")
print(readRDS("empirical_kudrin_replicate"))



source(paste0(root,"/code/Tests.R"))

run_ewk <- function(ps,ids=1:length(ps),pmax=0.15){
  pmin <- 0
  out <- list()
  
  out$lcms_EWK     <- LCM(ps,pmin,pmax)
  out$Fisher_EWK   <- Fisher(ps,pmin,pmax)
  out$disconts_EWK <- Discontinuity_test(ps,.05)
  out$binomial_EWK <- Binomial(ps,.04, .05, "c")
  out$CS1_EWK_s    <- CoxShi(ps, ids, pmin, pmax, 15, 1, 0)
  out$CS2_EWK_s    <- CoxShi(ps, ids, pmin, pmax, 15, 2, 1)
  
  out$min <- min(out$lcms_EWK,out$Fisher_EWK,out$disconts_EWK, out$binomial_EWK,  out$CS1_EWK_s, out$CS2_EWK_s)
  
  return(out)
}

rct_dropped     <- MM_dropped$method == "RCT"
iv_dropped     <- MM_dropped$method == "IV"
did_dropped     <- MM_dropped$method == "DID"
rdd_dropped     <- MM_dropped$method == "RDD"

ewk_rct_rep <- run_ewk(2 * pnorm(-abs(MM_dropped$t[rct_dropped])), ids = as.integer(factor(MM_dropped$title[rct_dropped])))
ewk_iv_rep  <- run_ewk(2 * pnorm(-abs(MM_dropped$t[iv_dropped])),  ids = as.integer(factor(MM_dropped$title[iv_dropped])))
ewk_did_rep <- run_ewk(2 * pnorm(-abs(MM_dropped$t[did_dropped])), ids = as.integer(factor(MM_dropped$title[did_dropped])))
ewk_rdd_rep <- run_ewk(2 * pnorm(-abs(MM_dropped$t[rdd_dropped])), ids = as.integer(factor(MM_dropped$title[rdd_dropped])))
empirical_ewk_rep <-cbind(ewk_rct_rep,ewk_iv_rep,ewk_did_rep,ewk_rdd_rep)
setwd(paste0(root,"/results"))
saveRDS(empirical_ewk_rep, file = paste0("EWK_replicate"))
print(readRDS("EWK_replicate"))
