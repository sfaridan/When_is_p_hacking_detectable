

#deprecated as of 5/2
run_test_slow<- function(data,numcoeffs,sigma_Y=1,shift=1.96,L=6.5,numgrid=3000,boots=150,U=NULL,seed=1){
  
  set.seed(seed)
  start_time <- Sys.time()
  
  #Projection basis
  if(is.null(U)){ #option to pass in U 
    U          <- create_basis(numcoeffs,L=L,numgrid=numgrid,sigma_Y=sigma_Y)
  }
  
  
  n          <- length(data$t)
  epsilon_U  <- sqrt(n)*compute_epsilons(L,nx=numgrid)$epsilon_U
  solver     <-  setup_projection_solver(U)
  
  ts                <- c(abs(data$t),-abs(data$t))-shift #symmetrize and then shift
  coeffs_orig       <- get_coeffs( ts,sigma_Y=sigma_Y,numcoeffs=numcoeffs)
  projection        <- compute_residual_fast(coeffs_orig,solver) 
  orig_resid        <- sqrt(n)*projection$residual 
  
  #The bootstrap
  projpoint           <- U%*%projection$alpha_opt
  resids_boot         <- rep(NA,boots)
  for (b in 1:boots){
    tic()
    warmstart <- 0*projection$alpha_opt
    
    # 1. Split the data into a list of data.frames, one per article:
    article_list <- split(data, data$title)
    
    # 2. Number of clusters (articles):
    m <- length(article_list)
    
    # 3. Resample m clusters with replacement:
    sampled_list <- sample(article_list, size = m, replace = TRUE)
    
    # 4. Re‐combine into one bootstrapped dataset:
    boot_data <- do.call(rbind, sampled_list)
    
    ts_boot                <- c(abs(boot_data$t),-abs(boot_data$t))-shift # symmeterize and center
    coeffs_boot            <- get_coeffs( ts_boot,sigma_Y=sigma_Y,numcoeffs=numcoeffs)
    
    #project onto tangent cone
    estar                  <- sqrt(n)*(coeffs_boot - coeffs_orig)
    
    
    #numerical tangent cone: line (25) of Fang and Santos (2019)
    sn                     <- n^(-1/3)
    proj_pertrubed         <- compute_residual_fast(coeffs_orig+estar*sn,solver,alpha_start = projection$alpha_opt)$residual
    resids_boot[b]         <- (proj_pertrubed-orig_resid/sqrt(n))/sn #numerical estimator of tangent cone
    
    toc()
    #print(paste0("boot ", b, " of ", boots))
    #print("Resid, 95%, 90%, eps, p, Bhat ")
    #print(c( orig_resid,quantile(resids_boot[1:b],0.95),quantile(resids_boot[1:b],0.90),epsilon_U,mean(orig_resid < resids_boot[1:b]+epsilon_U),(orig_resid - quantile(resids_boot[1:b],0.95)-epsilon_U)/sqrt(n)))
  }
  pval             <- mean(orig_resid < resids_boot+epsilon_U)
  breakdown        <- max(c(0,(orig_resid - quantile(resids_boot,0.95)-epsilon_U)/sqrt(n)))
  
  return(list(n=length(data$t),articles=length(unique(data$title)),resid=orig_resid,epsilon_U=epsilon_U,boot95=quantile(resids_boot,0.95),pval=pval,breakdown=breakdown,projpoint =projpoint,
              runtime = Sys.time()-start_time))
}

