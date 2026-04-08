# created 4 1 2026

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_sims.R"))
source(paste0(root,"/code/functions_application.R"))

#Load in the data
data_location <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data"
setwd(data_location)
MM_data <- haven::read_dta('MM data.dta')

#### Prep data

#Set the seed (for consistent de-rounding)
set.seed(1) 

#drop observations where the point estimate or standard error are missing (can't de-round) 
to_drop <-  is.na(MM_data$mu) |  is.na(MM_data$sd) # | MM_data$sd == 0  
to_drop[is.na(to_drop)] <- 1
sum(to_drop)

#de-round
MM_data$mu <- de_round(MM_data$mu)
MM_data$sd <- de_round(MM_data$sd)
MM_data$t_der <- MM_data$mu/MM_data$sd
MM_data$t[!is.na(MM_data$t_der)] <- MM_data$t_der[!is.na(MM_data$t_der)] #replace t with derounded whenever it exists

#de-round everything else
MM_data$t[is.na(MM_data$t_der)] <- de_round(MM_data$t[is.na(MM_data$t_der)])


### Run the tests

cv         <- 1.96
sigma_Y    <- 1 
shift      <- cv
L          <- 6.5 #grid support [-L,L] \cup 9999
numgrid    <- 3000 #total number of elements in grid (density)
boots      <- 200 #bootstrap repetitions


#Run the projection test
rct_results <- run_test(MM_data[MM_data$method=="RCT",],numcoeffs=20,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
did_results <- run_test(MM_data[MM_data$method=="DID",],numcoeffs=20,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
rdd_results <- run_test(MM_data[MM_data$method=="RDD",],numcoeffs=20,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
iv_results  <- run_test(MM_data[MM_data$method=="IV",] ,numcoeffs=20,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)


#Save results
empirical_results <- as.data.frame(do.call(rbind, list(rct_results, did_results,rdd_results, iv_results)))
empirical_results$methods <- c("RCT","DID","RDD","IV")
setwd(paste0(root,"/results"))
saveRDS(empirical_results, file = paste0("MM_results_",20,"_coeffs"))
print(readRDS(paste0("MM_results_",20,"_coeffs")))



#Run the projection test with more coefficients
coeffs<-30
rct_results <- run_test(MM_data[MM_data$method=="RCT",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
did_results <- run_test(MM_data[MM_data$method=="DID",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
rdd_results <- run_test(MM_data[MM_data$method=="RDD",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
iv_results  <- run_test(MM_data[MM_data$method=="IV",] ,numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots)
#Save results
empirical_results_30 <- as.data.frame(do.call(rbind, list(rct_results, did_results,rdd_results, iv_results)))
empirical_results_30$methods <- c("RCT","DID","RDD","IV")
setwd(paste0(root,"/results"))
saveRDS(empirical_results_30, file = paste0("MM_results_",coeffs,"_coeffs"))
print(readRDS(paste0("MM_results_",coeffs,"_coeffs")))

print(readRDS(paste0("MM_results_",20,"_coeffs")))
print(readRDS(paste0("MM_results_",30,"_coeffs")))


#Run EWK
source(paste0(root,"/code/Tests.R"))
ps <- pnorm(-abs(MM_data$t[MM_data$method=="RCT"]))
pmax <- 0.15
pmin <- min(ps[ps>0])
lcms_EWK <- LCM(ps,pmin,pmax)
Fisher_EWK <- Fisher(ps,pmin,pmax)
disconts_EWK <- Discontinuity_test(ps,.05)
binomial_EWK <- Binomial(ps,.04, .05, "c")
CS1_EWK_s <- CoxShi(ps, 1:length(ps), pmin, pmax, 30, 1, 0)
CS2_EWK_s <-CoxShi(ps, 1:length(ps), pmin, pmax, 30, 2, 1)
c(lcms_EWK,Fisher_EWK,disconts_EWK,binomial_EWK,CS1_EWK_s,CS2_EWK_s)