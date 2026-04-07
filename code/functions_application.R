


#returns number of decimal places in x for de-rounding
decimalplaces <- function(x) {
  options(scipen=999)
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

de_round <- function(x){
  
  for (i in 1:length(x)){
    if(!is.na(x[i])){
      decimals <- decimalplaces(x[i])
      x[i] <- 10^(-decimals)*( 10^decimals*x[i] +(runif(1)-0.5)  ) 
    }
  }
  
  return( x )
}

run_test<- function(data,numcoeffs,sigma_Y=1,shift=1.96,L=6.5,numgrid=3000,boots=150){

  #Projection basis
  U          <- create_basis(numcoeffs,L=L,numgrid=numgrid,sigma_Y=sigma_Y)
  delta      <- (2*L)/numgrid
  
  eps1<- max_diff_over_delta(delta,L=L,nx=numgrid)
  integrand2 <- function(x) {
    dnorm(x) * (dnorm(x-L))^2
  }
  eps2<-  sqrt(integrate(integrand2, lower = -10, upper = 10)$value)
  epsilon_U <- max(eps1,eps2)
  Uboot <- U
  solver <-  setup_projection_solver(U)
  
  ts        <- c(abs(data$t),-abs(data$t))-shift #symmeterize and then shift
  coeffs_ts <- get_coeffs( ts,sigma_Y=sigma_Y,numcoeffs=numcoeffs)
  projection <- compute_residual_fast(coeffs_ts,solver) 
  rh        <- projection$residual 
  
  #The bootstrap
  orig_resid    <- rh # test statistic
  orig_coeffs   <- coeffs_ts #coeffs_ts
  projpoint     <- U%*%projection$alpha_opt
  resids_boot   <- rep(NA,boots)
  for (b in 1:boots){
    
    # 1. Split the data into a list of data.frames, one per article:
    article_list <- split(data, data$title)
    
    # 2. Number of clusters (articles):
    m <- length(article_list)
    
    # 3. Resample m clusters with replacement:
    sampled_list <- sample(article_list, size = m, replace = TRUE)
    
    # 4. Re‐combine into one bootstrapped dataset:
    boot_data <- do.call(rbind, sampled_list)
    
    ts_boot                <- c(abs(boot_data$t),-abs(boot_data$t))-1.96
    coeffs_boot            <- get_coeffs( ts_boot,sigma_Y=sigma_Y,numcoeffs=numcoeffs)
    
    #sharper critical value
    v_boot                 <- projpoint+(coeffs_boot - orig_coeffs)
    resids_boot[b]         <- compute_residual_fast(v_boot,solver,alpha_start=projection$alpha_opt)$residual # compute_residual(v_boot,U)$residual
    
    #resids_boot[b]         <-  compute_residual(coeffs_boot,U)$residual-orig_resid
    print(paste0("boot ", b, " of ", boots))
    print(c( orig_resid,quantile(resids_boot[1:b],0.95),quantile(resids_boot[1:b],0.90)))
  }
  pval       <- mean(resids_boot+epsilon_U>orig_resid)
  breakdown  <- orig_resid - quantile(resids_boot,0.95)-epsilon_U
  
  return(list(n=length(data$t),articles=length(unique(data$title)),resid=orig_resid,epsilon_U=epsilon_U,boot95=quantile(resids_boot,0.95),pval=pval,breakdown=breakdown))
}



