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
  integrand <- function(x) {
    dnorm(x) * (dt(x, df = rct_ess_numbers[i]) - dnorm(x))^2
  }
  deltas_rcts[i] <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)
}
deltas_dids <- 0*did_ess_numbers
for (i in 1:length(did_ess_numbers)){
  integrand <- function(x) {
    dnorm(x) * (dt(x, df = did_ess_numbers[i]) - dnorm(x))^2
  }
  deltas_dids[i] <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)
  print(i)
}
deltas_ivs <- 0*iv_ess_numbers
for (i in 1:length(iv_ess_numbers)){
  integrand <- function(x) {
    dnorm(x) * (dt(x, df = iv_ess_numbers[i]) - dnorm(x))^2
  }
  deltas_ivs[i] <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)
}
c("RCTS: ", mean(deltas_rcts), ", DIDs: ", mean(deltas_dids), ", IVs: ",mean(deltas_ivs))
