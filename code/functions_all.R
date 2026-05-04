
library(osqp)
library(Matrix)
library(tictoc)

#The normalized generalized Hermite polynomials (Carrasco 2011)
hermite_general<- function(x,j,sigma_Y){
  ans <- 0*x
  for (l in 0:floor(j/2)){
    ans <- ans+(-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l)*(x/sigma_Y)^(j-2*l)
  }
  return(ans/sqrt(factorial(j)))
}

#takes in population without publication bias and truncates it
# if absolute value of t-score is above cv, you stay
# if |T|<cv then stay with probability theta
trunc_population <- function(untrunc_population, cv, theta){
  randkeep <- runif(length(untrunc_population))
  truncated <- untrunc_population[ abs(untrunc_population)>= cv | randkeep<=theta  ]
  return(truncated)
}

# project t-curve onto orthonormal basis polynomials
get_coeffs <- function(data,sigma_Y=1,numcoeffs =100){
  coeff_vec <- rep(0,numcoeffs)
  sig <- sigma_Y
  for (j in 0:(numcoeffs-1)){
    coeff_vec[j+1] <- mean( hermite_general(data,j,sigma_Y)*dnorm(data/sig)/sig )  
  }
  return(coeff_vec)
}

# Create the basis for the convex projection
create_basis <- function(num_coeffs,L=6,numgrid=1500,sigma_Y=1){
  pis <- c(seq(-L, L, length.out = numgrid),0,9999)
  U <- matrix(rep(NA,num_coeffs*length(pis)) ,nrow = num_coeffs)
  for (j in 1:num_coeffs){
    for(i in 1:length(pis)) {
      hh <- pis[i]
      U[j,i]= sqrt(sigma_Y^2/(1+sigma_Y^2))^(j-1)*dnorm(hh/sqrt(sigma_Y^2+1))/sqrt(sigma_Y^2+1)*hermite_general(hh,j-1,sqrt(sigma_Y^2+1))
    }
  }
  return(U)
}


## For given x,y compute:
## sqrt( ∫ φ(t) (φ(t - x) - φ(t - y))^2 dt )
kernel_diff <- function(x, y, rel.tol = 1e-8) {
  integrand <- function(t) {
    dnorm(t) * (dnorm(t - x) - dnorm(t - y))^2
  }
  val <- integrate(integrand, lower = -Inf, upper = Inf, rel.tol = rel.tol)$value
  sqrt(val)
}

## Approximate max_{|x - y| <= delta/2} of that quantity
## over a search box x,y ∈ [-L, L] with a grid of size nx × ny
max_diff_over_delta <- function(L = 6, nx = 3000,ny=100,
                                rel.tol = 1e-8) {
  
  delta      <- (2*L)/(nx-1)
  xs <- seq(-L, L, length.out = nx)
  max_val <- 0
  
  for (x in xs) {
    ys <- seq(max(-L, x - delta/2), min(L, x + delta/2), length.out = ny) #seq(x+delta/2, x-delta/2,length.out = ny)
    for (y in ys) {
      
      val <- kernel_diff(x, y, rel.tol = rel.tol)
      if (val > max_val) {
        max_val <- val
      }
      
    }
  }
  
  max_val
}

compute_eps2 <- function(L){
  integrand2 <- function(x) {
    dnorm(x) * (dnorm(x-L))^2
  }
  eps2<-  sqrt(integrate(integrand2, lower = -Inf, upper = Inf)$value)
  return(eps2)
}

compute_epsilons<- function(L, nx = 3000){
 eps1 <-  max_diff_over_delta(L=L,nx=nx)
  eps2 <- compute_eps2(L) # error from finite L
  return(list(eps1=eps1,eps2=eps2,epsilon_U = eps1+eps2))
  
}

draw_exp_matrix <- function(n, k, rate = 1) {
  
  if (!is.numeric(n) || !is.numeric(k) || n <= 0 || k <= 0 ||
      n != as.integer(n) || k != as.integer(k)) {
    stop("n and k must be positive integers.")
  }
  
  if (!is.numeric(rate) || rate <= 0) {
    stop("rate must be a positive number.")
  }
  
  # Generate n * k samples from Exponential(rate), and convert to matrix
  matrix(rexp(n * k, rate = rate), nrow = n, ncol = k)
}

draw_unif_matrix <- function(n, k) {
  
  if (!is.numeric(n) || !is.numeric(k) || n <= 0 || k <= 0 || n != as.integer(n) || k != as.integer(k)) {
    stop("n and k must be positive integers.")
  }
  
  # Generate n * k samples from the uniform, and convert to matrix
  matrix(runif(n * k, min = -sqrt(3),max=sqrt(3)), nrow = n, ncol = k)
}


compute_residual_fast <- function(v, solver, alpha_start = NULL) {
  q <- as.numeric(-2 * crossprod(solver$U, v))
  
  solver$model$Update(q = q)   # model@Update(...) in osqp >= 1.0
  if (!is.null(alpha_start)) {
    solver$model$WarmStart(x = alpha_start)
  }
  
  result <- solver$model$Solve()
  alpha  <- result$x
  
  # residual from quadratic form (see below)
  cvec     <- as.numeric(crossprod(solver$U, v))
  resid_sq <- drop(crossprod(v)) - 2 * sum(cvec * alpha) +
    sum(alpha * as.numeric((crossprod(solver$U) %*% alpha)))
  list(residual = sqrt(max(resid_sq, 0)), alpha_opt = alpha)
}


setup_projection_solver <- function(U) {
  nvar <- ncol(U)
  
  P_sp <- as(2 * crossprod(U), "dgCMatrix")
  A    <- rbind(Matrix(1, 1, nvar), Diagonal(nvar))
  A_sp <- as(as(A, "TsparseMatrix"), "dgCMatrix")
  l    <- c(1, rep(0, nvar))
  u    <- c(1, rep(Inf, nvar))
  
  pars <- osqpSettings(
    verbose = FALSE,
    warm_start = TRUE   # or warm_starting = TRUE depending on version
  )
  
  model <- osqp(P = P_sp, q = rep(0, nvar), A = A_sp, l = l, u = u, pars = pars)
  
  list(
    U = U,
    P = P_sp,
    model = model
  )
}


t_to_p_normal_approx <- function(t_score) {
  # Input validation
  if (!is.numeric(t_score)) {
    stop("t_score must be numeric.")
  }
  
  # Use the standard normal distribution to approximate the two-tailed p-value
  p_value <- 2 * pnorm(-abs(t_score))
  return(p_value)
}

run_sims<- function(parms,sim_file_prefix="sim_parms_"){
  tic()
  num_parameterizations <- nrow(parms)
  
  #Loop over parameterizations
  for (parm in 1:num_parameterizations) {
    
    #Set the seed
    set.seed(parms$seed[parm])
    
    #Read in the parameterization
    n          <- parms$n[parm] #5000    # meta-sample size
    cv         <- parms$cv[parm] #1.96    # critical value to shift by
    theta      <- parms$theta[parm] #0.7     # probability of reporting t when |t|<cv
    sigma_Y    <- 1 
    num_coeffs <- parms$num_coeffs[parm] #20
    nsims      <- parms$nsims[parm] #200 #simulations
    nboots     <- parms$nboots[parm] #400 #bootstrap repetitions
    nu         <- parms$nu[parm] #9999999
    numgrid    <- parms$numgrid[parm] #1500
    L          <- parms$L[parm] #6.5
    h_center   <- parms$h_center[parm] # 2
    sigma_h    <- parms$sigma_h[parm]
    bimodal    <- parms$bimodal[parm]
    pi0_shape        <- parms$pi0_shape[parm]
    prob_hack  <- parms$prob_hack[parm] # 2
    hack_type  <- parms$hack_type[parm]
    
    start_time <- Sys.time()

    
   
    integrand <- function(x) { dnorm(x) * (dt(x, df = nu) - dnorm(x))^2 }
    nu_resid  <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value) #this is Delta if we assume t~student-t(nu)
    
    
    
    #Precompute to accelerate the projection test
    U             <- create_basis(num_coeffs,L=L,numgrid=numgrid,sigma_Y=sigma_Y)
    solver        <- setup_projection_solver(U)
    epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U
    
    #Set up Elliott et al pval vectors
    lcms_EWK        <- rep(NA,nsims)
    disconts_EWK    <- rep(NA,nsims)
    CS1_EWK         <- rep(NA,nsims)
    Fisher_EWK      <- rep(NA,nsims)
    CS2B_EWK        <- rep(NA,nsims)
    binomial_EWK    <- rep(NA,nsims)
    
    #Begin the simulation loop
    resids    <- rep(NA,nsims)
    boot95    <- rep(NA,nsims)
    rejects   <- rep(NA,nsims)
    for(sim in 1:nsims){
      
      sim_start  <- Sys.time()
        print("******************************")
        print(paste0("Parm: ", parm, " of ",num_parameterizations, ", Sim: ", sim, " of ", parms$nsims[parm] ))
        print(Sys.time())

        
      
      rands   <- runif(n) #random numbers to determine who p-hacks
      
      
      if(pi0_shape == "normal") { hnoise <- rnorm(n)}
      if(pi0_shape == "uniform"){ hnoise <- unif(n)}
      if(pi0_shape == "point")  { hnoise <- rep(0,n)}
      if(pi0_shape == "chi2")   { hnoise <- rnorm(n)^2}
      
      if(bimodal){  
        rr <- runif(n)
        hnoise[rr < 0.5]   <- hnoise[rr < 0.5]  - 1
        hnoise[rr >= 0.5]  <- hnoise[rr >= 0.5] + 1
      }
      
      hs      <- h_center + hnoise*sigma_h # latent true effects
      
      #Draw individual t-scores
        ts1     <- rt(n,nu)+hs
        ts2     <- rt(n,nu)+hs
      
      if(hack_type == "max"){ #maximization p-hacking with no threshold 
        ts3     <- rt(n,nu)+hs
        ts4     <- rt(n,nu)+hs
        ts5     <- rt(n,nu)+hs
        ts6     <- rt(n,nu)+hs
        ts7     <- rt(n,nu)+hs
        ts8     <- rt(n,nu)+hs
        ts_pre<- pmax(ts1,ts2,ts3,ts4,ts5,ts6,ts7,ts8 )*(1*(rands<= prob_hack))+ts1*(rands>prob_hack)  #maximization p-hacking that doesn't add any discontinuities
      }
      if(hack_type=="threshold"){
        ts_pre <- ts1*(rands>prob_hack | (ts1)>cv)+(rands<= prob_hack & (ts1)<=cv)*pmax(ts1,ts2 ) #Draw 1 t-score. Report it if you don't phack or if it is significant and positive. Otherwise draw a second and report the max of the two
      }
      
      if(theta != 1){  #putlication bias if theta < 1
        ts_pre <- trunc_population(ts_pre,cv,theta) 
      }
       
      ts <-  ts_pre # finalize observed t-scores
      
      print(paste0("H: ", pi0_shape, ", center=", h_center, ", sigma_h=",sigma_h,", bim=",bimodal))
      print(summary(hs))
      print(paste0("frac hacking: ", prob_hack,", hacktype: ", hack_type))
      print(summary(ts))
      
      #Run tests from Elliott et al.
      if(parms$omit_EWK[parm]==FALSE){
        ewk_start  <- Sys.time()
        
        ps <- t_to_p_normal_approx(ts_pre) #convert to p-values
        
        maxp <- 0.15
        pmin <- min(ps[ps>0])
        lcms_EWK[sim] <- LCM(ps,pmin,maxp)
        Fisher_EWK[sim] <- Fisher(ps,pmin,maxp)
        disconts_EWK[sim] <- Discontinuity_test(ps,.05)
        binomial_EWK[sim] <- Binomial(ps,.04, .05, "c")
        
        #Cox-shi tests occasionally fail to converge. We ignore these simulation iterations.
        #CoxShi(ps, 1:n, pmin, maxp, 30, 1, 0) #for debugging to see why CoxShi sometimes fails to converge
        CS1_EWK_s <- try({CoxShi(ps, 1:length(ps), pmin, maxp, 30, 1, 0)},silent=TRUE)
        if(inherits(CS1_EWK_s,"try-error")){  
          CS1_EWK[sim]<- NA
          print("cs1 NA")
        }else{CS1_EWK[sim] <- CS1_EWK_s}
        CS2_EWK_s <- try({ CoxShi(ps, 1:length(ps), pmin, maxp, 30, 2, 1)},silent=TRUE)
        if(inherits(CS2_EWK_s,"try-error")){  
          CS2B_EWK[sim]<- NA
          print("cs2 NA")
        }else{
          CS2B_EWK[sim] <- CS2_EWK_s
        }
        print(paste0("EWK time: ", Sys.time() -ewk_start ))
      }
      
      #Projection method
      if(parms$omit_proj[parm]==FALSE){
        proj_start  <- Sys.time()
        
        data        <- data.frame(t=ts,title = as.character(1:length(ts)))
        sim_results <- run_test(data,num_coeffs,sigma_Y=sigma_Y,shift=cv,L=L,numgrid=numgrid,boots=nboots,
                                U=U, solver=solver, epsilon_U_pre=epsilon_U_pre)
        
        epsilon_U   <- sim_results$epsilon_U
        resids[sim] <- sim_results$resid
        boot95[sim] <- sim_results$boot95
        
        rejects[sim] <- (sim_results$pval <= .05) 
        print(paste0("proj time: ", Sys.time() -proj_start ))
      }
      
      
      if(sim / 1 == floor(sim/1)){
        #print(paste0("parm: ", parm, " of ", num_parameterizations,", sim: ", sim, " of ", nsims))
        print(parms[parm,])
        print( paste0("Rejections: CS1: ", mean(CS1_EWK[1:sim]<.05,na.rm=TRUE), ", CS2B: ", mean(CS2B_EWK[1:sim]<.05,na.rm=TRUE), ", proj: ", mean(rejects[1:sim])))
        }
  
      print(paste0("sim time: ", Sys.time() -sim_start ))
      
    } # ends sim loop
    
    rejects_no_nu <- resids > boot95+epsilon_U
    
    #Record the results
    parms$epsilon_U[parm]          <- epsilon_U
    parms$nu_resid[parm]           <- nu_resid
    parms$reject_rate_no_adj[parm] <- mean(resids > boot95)
    parms$reject_rate_no_nu[parm]  <- mean(rejects_no_nu)
    parms$reject_rate[parm]        <- mean(rejects)
    parms$mean_resid[parm]         <- mean(resids)
    parms$breakdown[parm]          <- mean( resids[rejects_no_nu] - boot95[rejects_no_nu]-epsilon_U )
    parms$boot95[parm]             <- mean(boot95)
    parms$mean_cvs_no_nu[parm]     <- mean(boot95+epsilon_U)
    parms$mean_cvs_full[parm]      <- mean(boot95+epsilon_U+nu_resid)
    parms$finished_at[parm]        <-  format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    parms$time_taken[parm]         <-  Sys.time()-start_time
    parms$script_path[parm]        <-  this_script
    
    #Record Elliott et al. results
    parms$rej_lcm_EWK[parm]  <- mean(lcms_EWK < .05)
    parms$Fisher_EWK[parm]   <- mean(Fisher_EWK < .05)
    parms$disconts_EWK[parm] <- mean(disconts_EWK < .05)
    parms$binomial_EWK[parm] <- mean(binomial_EWK < .05)
    parms$CS1_EWK[parm]      <- mean(CS1_EWK < .05,na.rm=TRUE)
    parms$CS2B_EWK[parm]     <- mean(CS2B_EWK < .05,na.rm=TRUE)
    
    # Save full results of this parameterization
    setwd(paste0(simulation_results,"/raw"))
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    filename <- sprintf("sim_raw_%s_nu%s_coeffs%d_theta%.2f.rds", 
                        timestamp,nu, num_coeffs, theta)
    results <- list(
      resids = resids,
      boot95 = boot95,
      rejects = rejects
    )
    saveRDS(results, file = filename)
    results <- readRDS(filename)
    
    print(parms[1:parm,])
    print(paste0(parm, " out of ", num_parameterizations))
    
    # Save summary results every parameterization
    setwd(simulation_results)
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    filename <- paste0(sim_file_prefix,dim(parms)[1],"_nsims",parms$nsims[1],"_",timestamp,".rds")
    saveRDS(parms, file = filename)
    print("******************************")
    print("Results here:")
    print(filename)
    
  } # ends parm loop
  
  
  # Save summary results for all parameterization
  setwd(simulation_results)
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  filename <- paste0(sim_file_prefix,dim(parms)[1],"_nsims",parms$nsims[1],"_",timestamp,".rds")
  saveRDS(parms, file = filename)
  print("******************************")
  print("Results here:")
  print(filename)
  
 
  return(filename)
}



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


get_coeff_matrix <- function(data, sigma_Y = 1, numcoeffs = 100) {
  sig <- sigma_Y
  out <- matrix(NA_real_, nrow = length(data), ncol = numcoeffs)
  w   <- dnorm(data / sig) / sig
  for (j in 0:(numcoeffs - 1)) {
    out[, j + 1] <- hermite_general(data, j, sigma_Y) * w
  }
  out
}


run_test <- function(data, numcoeffs, sigma_Y = 1, shift = 1.96,
                     L = 6.5, numgrid = 3000, boots = 150, verbose=FALSE,
                     U = NULL, solver=NULL,epsilon_U_pre=NULL) {
  
  #set.seed(seed)
  start_time <- Sys.time()
  
  n <- length(data$t)
  
  if (is.null(U)) {
    U <- create_basis(numcoeffs, L = L, numgrid = numgrid, sigma_Y = sigma_Y)
  }
  if (is.null(solver)) {
    solver <- setup_projection_solver(U)
  }
  if (is.null(epsilon_U_pre)) {
    epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U
  }
  
  epsilon_U <- sqrt(n) * epsilon_U_pre
  
  # Symmetrize and shift once
  ts <- c(abs(data$t), -abs(data$t)) - shift
  
  # Precompute all coefficient contributions once
  coeff_mat   <- get_coeff_matrix(ts, sigma_Y = sigma_Y, numcoeffs = numcoeffs)
  coeffs_orig <- colMeans(coeff_mat)
  
  projection <- compute_residual_fast(coeffs_orig, solver)
  orig_resid <- sqrt(n) * projection$residual
  
  projpoint    <- U %*% projection$alpha_opt

  # Precompute article -> rows mapping in coeff_mat
  article_index <- split(seq_len(n), data$title)
  article_rows  <- lapply(article_index, function(idx) c(idx, idx + n))
  m <- length(article_rows)
  
  tic()
  resids_boot <- rep(NA, boots)
  for (b in 1:boots) {

    sampled_articles <- sample.int(m, size = m, replace = TRUE)
    boot_rows <- unlist(article_rows[sampled_articles], use.names = FALSE)
    
    coeffs_boot <- colMeans(coeff_mat[boot_rows, , drop = FALSE])
    
    estar <- sqrt(n) * (coeffs_boot - coeffs_orig)
    
    sn             <- n^(-1/3)
    proj_pertrubed <- compute_residual_fast(
      coeffs_orig + estar * sn,
      solver,
      alpha_start = projection$alpha_opt
    )$residual
    
    resids_boot[b] <- (proj_pertrubed - orig_resid / sqrt(n)) / sn

    if(verbose){
      print(paste0("boot ", b, " of ", boots))
      print("Resid, 95%, 90%, eps, p, Bhat ")
      print(c( orig_resid,quantile(resids_boot[1:b],0.95),quantile(resids_boot[1:b],0.90),epsilon_U,mean(orig_resid < resids_boot[1:b]+epsilon_U),(orig_resid - quantile(resids_boot[1:b],0.95)-epsilon_U)/sqrt(n)))
    }
  }

  
  pval      <- mean(orig_resid < resids_boot + epsilon_U)
  breakdown <- max(c(0, (orig_resid - quantile(resids_boot, 0.95) - epsilon_U) / sqrt(n)))
  
  return(list(
    n = length(data$t),
    articles = length(unique(data$title)),
    resid = orig_resid,
    epsilon_U = epsilon_U,
    boot95 = quantile(resids_boot, 0.95),
    pval = pval,
    breakdown = breakdown,
    projpoint = projpoint,
    runtime = Sys.time()-start_time
  ))
}



