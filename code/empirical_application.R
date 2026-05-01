# created 4 1 2026

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_all.R"))

#Load in the data
data_location <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data"
setwd(data_location)
MM_data <- haven::read_dta('MM data.dta')

#### Prep data

#Set the seed (for consistent de-rounding)
set.seed(1) 

#de-round via mus and sd
MM_data$mu    <- de_round(MM_data$mu)
MM_data$sd    <- de_round(MM_data$sd)
MM_data$t_der <- MM_data$mu/MM_data$sd
MM_data$t[!is.na(MM_data$t_der)] <- MM_data$t_der[!is.na(MM_data$t_der)] #replace t with derounded whenever it exists

#de-round the raw t-score itself when mu, sd are not available
MM_data$t[is.na(MM_data$t_der)] <- de_round(MM_data$t[is.na(MM_data$t_der)])

if (sum(is.infinite(MM_data$t))+sum(is.na(MM_data$t))+sum(is.nan(MM_data$t)) > 0){
  stop("Tscores not clean")
}




### Run the tests

cv         <- 1.96
sigma_Y    <- 1 
shift      <- cv
L          <- 6.5 #grid support [-L,L] \cup 9999
numgrid    <- 3000 #total number of elements in grid (density)
boots      <- 50 #bootstrap repetitions

#Run the projection test
J <- 20
coeffs <- J+1
rct_results <- run_test(MM_data[MM_data$method=="RCT",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
did_results <- run_test(MM_data[MM_data$method=="DID",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
rdd_results <- run_test(MM_data[MM_data$method=="RDD",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
iv_results  <- run_test(MM_data[MM_data$method=="IV",] ,numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)


#Save results
empirical_results <- as.data.frame(do.call(rbind, list(rct_results, did_results,rdd_results, iv_results)))
empirical_results$methods <- c("RCT","DID","RDD","IV")
setwd(paste0(root,"/results"))
name <- paste0("MM_tangcone_results_",J,"_J")
saveRDS(empirical_results, file = name)
print(readRDS(name))



#Run the projection test with more coefficients
J      <- 30
coeffs <- J+1
rct_results <- run_test(MM_data[MM_data$method=="RCT",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
did_results <- run_test(MM_data[MM_data$method=="DID",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
rdd_results <- run_test(MM_data[MM_data$method=="RDD",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
iv_results  <- run_test(MM_data[MM_data$method=="IV",] ,numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
#Save results
empirical_results_30 <- as.data.frame(do.call(rbind, list(rct_results, did_results,rdd_results, iv_results)))
empirical_results_30$methods <- c("RCT","DID","RDD","IV")
setwd(paste0(root,"/results"))
name <- paste0("MM_tangcone_results_",J,"_J")
saveRDS(empirical_results_30, file = name)
print(readRDS(name))



#Run EWK


source(paste0(root,"/code/Tests.R"))

run_ewk <- function(ps,pmax=0.15){
  pmin <- min(ps[ps>0])
  out <- list()
  
  out$lcms_EWK     <- LCM(ps,pmin,pmax)
  out$Fisher_EWK   <- Fisher(ps,pmin,pmax)
  out$disconts_EWK <- Discontinuity_test(ps,.05)
  out$binomial_EWK <- Binomial(ps,.04, .05, "c")
  out$CS1_EWK_s    <- CoxShi(ps, 1:length(ps), pmin, pmax, 30, 1, 0)
  out$CS2_EWK_s    <- CoxShi(ps, 1:length(ps), pmin, pmax, 30, 2, 1)
  
  out$min <- min(out$lcms_EWK,out$Fisher_EWK,out$disconts_EWK, out$binomial_EWK,  out$CS1_EWK_s, out$CS2_EWK_s)
  
  return(out)
}
ewk_rct <- run_ewk(2*pnorm(-abs(MM_data$t[MM_data$method=="RCT"])))
ewk_iv  <- run_ewk(2*pnorm(-abs(MM_data$t[MM_data$method=="IV"])))
ewk_did <- run_ewk(2*pnorm(-abs(MM_data$t[MM_data$method=="DID"])))
ewk_rdd <- run_ewk(2*pnorm(-abs(MM_data$t[MM_data$method=="RDD"])))

empirical_ewk <-cbind(ewk_rct,ewk_iv,ewk_did,ewk_rdd)
setwd(paste0(root,"/results"))
saveRDS(empirical_ewk, file = paste0("EWK_results"))