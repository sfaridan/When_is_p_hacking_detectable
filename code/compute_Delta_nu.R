
Delta_nu_maxh <- function(nu){
  hgrid <- c((0:1000)/1000,9999)
  ints <- rep(0,length(hgrid))
  for (hh in 1:length(hgrid)){
    h <- hgrid[hh]
    integrand <- function(x) {
      dnorm(x) * (dt(x-h, df = nu) - dnorm(x-h))^2
    }
    ints[hh] <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)
  }
  return(max(ints))
}



# Define the integrand function
integrand <- function(x) {
  dnorm(x) * (dt(x, df = 14) - dnorm(x))^2
}
# Perform the numerical integration
result <- integrate(integrand, lower = -Inf, upper = Inf)
# Take the square root
sqrt(result$value)

#Read in effective sample sizes
root <- "C:/Users/sfaridani6/Documents/R/smoothness test for p-hacking"
rct_ess <- read.csv(paste0(root,"/Effective Sample Sizes for Brodeur et al (2020) - RCTs.csv"))
rct_ess$number[rct_ess$number=="40,000"] <- 40000
rct_ess_numbers <- as.numeric(rct_ess[rct_ess$number != "","number"])
did_ess <- read.csv(paste0(root,"/Effective Sample Sizes for Brodeur et al (2020) - DIDs.csv"))
iv_ess <- read.csv(paste0(root,"/Effective Sample Sizes for Brodeur et al (2020) - IVs.csv"))
did_ess$numerical[did_ess$numerical=="300 SEAs in the 2001 DHS and more than 300 SEAs in the 2007"]<- 300
did_ess_numbers <- as.numeric(gsub("?","",gsub(",","",did_ess[1:100,"numerical"])))
did_ess_numbers <- did_ess_numbers[!is.na(did_ess_numbers)]
iv_ess$numerical[iv_ess$numerical=="192 countries 1243 goods"]<- 192
iv_ess_numbers <- as.numeric(gsub("?","",gsub(",","",iv_ess[1:100,"numerical"])))
iv_ess_numbers <- iv_ess_numbers[!is.na(iv_ess_numbers)]

deltas_rcts <- 0*rct_ess_numbers
for (i in 1:length(rct_ess_numbers)){
  deltas_rcts[i] <-Delta_nu_maxh(rct_ess_numbers[i])
  print(c(i, deltas_rcts[i] ))
}
deltas_dids <- 0*did_ess_numbers
for (i in 1:length(did_ess_numbers)){
  deltas_dids[i] <-Delta_nu_maxh(did_ess_numbers[i])
}
deltas_ivs <- 0*iv_ess_numbers
for (i in 1:length(iv_ess_numbers)){
  deltas_ivs[i] <-Delta_nu_maxh(iv_ess_numbers[i])
}

c("RCTS: ", mean(deltas_rcts), ", DIDs: ", mean(deltas_dids), ", IVs: ",mean(deltas_ivs))
c("RCTS: ", median(rct_ess_numbers), ", DIDs: ", median(did_ess_numbers), ", IVs: ",median(iv_ess_numbers))

Delta_nu_maxh(120)
