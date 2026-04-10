
library(osqp)
library(Matrix)
library(tictoc)

#The normalized generalized Hermite polynomials (Carrasco 2011)
hermite_general<- function(x,j,sigma_Y){
  ans <- 0*x
  for (l in 0:floor(j/2)){
    ans <- ans+(-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l)*(x/sigma_Y)^(j-2*l)
    #print((-1)^l*factorial(2*l)/2^l/factorial(l)*choose(j,2*l))
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
  #L <- 6
  pis <- c(seq(-L, L, length.out = numgrid),0,9999)
  U <- matrix(rep(NA,num_coeffs*length(pis)) ,nrow = num_coeffs)
  for (j in 1:num_coeffs){
    for(i in 1:length(pis)) {
      hh <- pis[i]
      U[j,i]= sqrt(sigma_Y^2/(1+sigma_Y^2))^(j-1)*dnorm(hh/sqrt(sigma_Y^2+1))/sqrt(sigma_Y^2+1)*hermite_general(hh,j-1,sqrt(sigma_Y^2+1))
    }
    #print(paste0("Basis row ", j, " of ", num_coeffs," created"))
  }
  return(U)
}

phi <- function(t) dnorm(t)

## For given x,y compute:
## sqrt( ∫ φ(t) (φ(t - x) - φ(t - y))^2 dt )
kernel_diff <- function(x, y, rel.tol = 1e-8) {
  integrand <- function(t) {
    phi(t) * (phi(t - x) - phi(t - y))^2
  }
  val <- integrate(integrand, lower = -Inf, upper = Inf, rel.tol = rel.tol)$value
  sqrt(val)
}

## Approximate max_{|x - y| <= delta/2} of that quantity
## over a search box x,y ∈ [-L, L] with a grid of size nx × ny
max_diff_over_delta <- function(delta, L = 6, nx = 3000,ny=100,
                                rel.tol = 1e-8) {
  xs <- seq(-L, L, length.out = nx)
  max_val <- 0
  
  for (x in xs) {
    ys <- seq(x+delta/2, x-delta/2,length.out = ny)
    for (y in ys) {
      
      val <- kernel_diff(x, y, rel.tol = rel.tol)
      if (val > max_val) {
        max_val <- val
      }
      
    }
  }
  
  max_val
}

draw_exp_matrix <- function(n, k, rate = 1) {
  if (!is.numeric(n) || !is.numeric(k) || n <= 0 || k <= 0) {
    stop("n and k must be positive integers.")
  }
  if (!is.numeric(rate) || rate <= 0) {
    stop("rate must be a positive number.")
  }
  
  # Generate n * k samples from Exponential(rate), and convert to matrix
  matrix(rexp(n * k, rate = rate), nrow = n, ncol = k)
}

draw_unif_matrix <- function(n, k, rate = 1) {
  if (!is.numeric(n) || !is.numeric(k) || n <= 0 || k <= 0) {
    stop("n and k must be positive integers.")
  }
  if (!is.numeric(rate) || rate <= 0) {
    stop("rate must be a positive number.")
  }
  
  # Generate n * k samples from Exponential(rate), and convert to matrix
  matrix(runif(n * k, min = -sqrt(3),max=sqrt(3)), nrow = n, ncol = k)
}

#deprecated on 3/31/2026
compute_residual_old <- function(v,U){
  # Assume U (K×N) and v (length K) are defined...
  P_sp <- as(2 * crossprod(U), "dgCMatrix")
  q    <- -2 * crossprod(U, v)
  n    <- ncol(U)
  
  # Build A_sp via two‑step coercion:
  A    <- rbind(Matrix(1,1,n), Diagonal(n))
  A_sp <- as(as(A, "TsparseMatrix"), "dgCMatrix")
  
  # OR build A_sp directly:
  # eq   <- Matrix(1, 1, n, sparse = TRUE)
  # diag <- Diagonal(n)
  # A_sp <- cBind(eq, diag)
  
  l <- c(1, rep(0, n))
  u <- c(1, rep(Inf, n))
  
  Settings<- osqpSettings(verbose=FALSE)
  
  model  <- osqp(P = P_sp, q = q, A = A_sp, l = l, u = u,pars=Settings)
  result <- model$Solve()
  
  alpha_opt <- result$x
  residual  <- sqrt(sum((U %*% alpha_opt - v)^2))
  out <- list(residual=residual, alpha_opt=alpha_opt)
  return(out)
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
    sigma_Y    <- parms$sigma_Y[parm] #1 
    num_coeffs <- parms$num_coeffs[parm] #20
    nsims      <- parms$nsims[parm] #200 #simulations
    nboots     <- parms$nboots[parm] #400 #bootstrap repetitions
    nu         <- parms$nu[parm] #9999999
    numgrid    <- parms$numgrid[parm] #1500
    L          <- parms$L[parm] #6.5
    h          <- parms$h[parm] # 2
    prob_hack  <- parms$prob_hack[parm] # 2
    
    start_time <- Sys.time()
    
    # Pre-compute important quantities 
    eta_1      <- sqrt(sigma_Y^2/(1+sigma_Y^2))
    
    if(parms$expo[parm]){
      nu_resid <- .0745/nu #2*.378/sqrt(nu)/6
    }
    else{
      integrand <- function(x) { dnorm(x) * (dt(x, df = nu) - dnorm(x))^2 }
      nu_resid  <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value) #this is Delta if we assum t~student-t(nu)
    }
    
    
    #Projection basis
    U          <- create_basis(num_coeffs,L=L,numgrid=numgrid,sigma_Y=sigma_Y)
    delta      <- (2*L)/numgrid
    
    #not quite right
    #integrand <- function(x) {
    #  dnorm(x) * (dnorm(x-1-delta /2) -dnorm(x-1))^2
    #}
    #eps1<-  sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value)  # (dnorm(1-L/numgrid)-dnorm(1+L/numgrid))/( 2*pi )^(1/4) #approx error from grid density
    
    eps1<- max_diff_over_delta(delta,L=L,nx=numgrid)
    
    integrand2 <- function(x) {
      dnorm(x) * (dnorm(x-L))^2
    }
    eps2<-  sqrt(integrate(integrand2, lower = -10, upper = 10)$value)  #sqrt( dnorm(L/sigma_Y)/sigma_Y/sqrt(2*pi) )/2 #approx error from grid support
    c(eps1,eps2) #the approximation error is the maximum of these two
    epsilon_U <- max(eps1,eps2) # approximation error of projection basis U
    
    
    #True dgp: T ~N(h,1)=h+Z
    #Calculate the true coefficients for this dgp
    null_true_coeffs <- rep(NA,num_coeffs)
    for(j in 0:(num_coeffs-1)){
      #no symmeterization
      #null_true_coeffs[j+1] <- eta_1^(j)*hermite_general(h-cv,j,sqrt(sigma_Y^2+1) )*dnorm((h-cv)/sqrt(sigma_Y^2+1))/sqrt(sigma_Y^2+1)
      #with symmeterization
      null_true_coeffs[j+1] <- 0.5*eta_1^(j)*hermite_general(h-cv,j,sqrt(sigma_Y^2+1) )*dnorm((h-cv)/sqrt(sigma_Y^2+1))/sqrt(sigma_Y^2+1)+0.5*eta_1^(j)*hermite_general(-h-cv,j,sqrt(sigma_Y^2+1) )*dnorm((-h-cv)/sqrt(sigma_Y^2+1))/sqrt(sigma_Y^2+1)
    }
    
    # Speed up the bootstrap at the cost of making the cvs bigger
    # (only do this if this isn't a check for size)
    Uboot <- U
    #if (theta <1 & parms$fastboot[parm]==TRUE  ){Uboot <- create_basis(num_coeffs,L=L,numgrid=(numgrid/2),sigma_Y=sigma_Y)}
    
    #For warm-starting
    solver <-  setup_projection_solver(U)
    
    #Set up Elliott et al pval vectors
    lcms_EWK        <- rep(NA,nsims)
    disconts_EWK    <- rep(NA,nsims)
    CS1_EWK      <- rep(NA,nsims)
    Fisher_EWK     <- rep(NA,nsims)
    CS2B_EWK        <- rep(NA,nsims)
    binomial_EWK        <- rep(NA,nsims)
    
    resids    <- rep(NA,nsims)
    boot95    <- rep(NA,nsims)
    rejects   <- rep(NA,nsims)
    true_dist <- rep(NA,nsims)
    proj_error <- rep(NA,nsims) # distance between projection point and true coefficients
    
    for(sim in 1:nsims){
      
      
      tic()
      numsimsper <- 100
      if (TRUE){ #( (sim) /numsimsper == floor( (sim) /numsimsper) ){
        print("")
        print("")
        print(paste0("Parm: ", parm, " of ",num_parameterizations, ", Sim: ", sim, " of ", parms$nsims[parm] ))
        print(Sys.time())
        toc()
        tic()
      }
      
      #Draw sample
      #ts_pre <-  rnorm(n)+h         #draw the sample
      rands   <- runif(n)
      
      hs <- h +rnorm(n)*parms$sigma_h[parm]
      
      #Instead, make ts not exactly normal
      if(parms$expo[parm]){
        #ts1 <-  h - (draw_exp_matrix(n,nu,rate=1)-1)%*%rep(1/sqrt(nu),nu)
        #ts2 <-  h - (draw_exp_matrix(n,nu,rate=1)-1)%*%rep(1/sqrt(nu),nu)
        ts1 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        ts2 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        ts3 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        ts4 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        ts5 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        ts6 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        ts7 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        ts8 <-  hs + (draw_unif_matrix(n,nu))%*%rep(1/sqrt(nu),nu)
        #print(mean((ts1-h)^3) )
      }
      else{
        ts1     <- rt(n,nu)+hs
        ts2     <- rt(n,nu)+hs
        ts3     <- rt(n,nu)+hs
        ts4     <- rt(n,nu)+hs
        ts5     <- rt(n,nu)+hs
        ts6     <- rt(n,nu)+hs
        ts7     <- rt(n,nu)+hs
        ts8     <- rt(n,nu)+hs
      }
      #print(mean((ts1-h)^3) )
      
      #noise  <- rt(n,nu)*(rands>prob_hack)+(rands<= prob_hack)*pmax(rt(n,nu),rt(n,nu) )
      
      
      if(parms$smooth_hack[parm]){ #maximization p-hacking with no threshold and h~N(0,1)
        ts_pre<- pmax(ts1,ts2,ts3,ts4,ts5,ts6,ts7,ts8 )*(1*(rands<= prob_hack))+ts1*(rands>prob_hack)  #maximization p-hacking that doesn't add any discontinuities
      }
      else{
        ts_pre <- ts1*(rands>prob_hack | abs(ts1)>cv)+(rands<= prob_hack & abs(ts1)<=cv)*pmax(ts1,ts2 )
      }
      
      ts_pre <- trunc_population(ts_pre,cv,theta)
      ts <- c(abs(ts_pre),-abs(ts_pre))-cv*parms$shift[parm] #symmeterize and re-center
      ps <- t_to_p_normal_approx(ts_pre)
      
      #Run tests from Elliott et al.
      if(parms$omit_EWK[parm]==FALSE){
        print("Running EWK:")
        
        
        pmax <- 0.15
        pmin <- min(ps[ps>0])
        lcms_EWK[sim] <- LCM(ps,pmin,pmax)
        Fisher_EWK[sim] <- Fisher(ps,pmin,pmax)
        disconts_EWK[sim] <- Discontinuity_test(ps,.05)
        binomial_EWK[sim] <- Binomial(ps,.04, .05, "c")
        
        #Cox-shi tests occaisonally fail to converge. We ignore these simulation iterations.
        #CoxShi(ps, 1:n, pmin, pmax, 30, 1, 0) #for debugging to see why CoxShi sometimes fails to converge
        CS1_EWK_s <- try({CoxShi(ps, 1:n, pmin, pmax, 30, 1, 0)},silent=TRUE)
        if(inherits(CS1_EWK_s,"try-error")){  
          CS1_EWK[sim]<- NA
          print("cs1 NA")
        }else{CS1_EWK[sim] <- CS1_EWK_s}
        CS2_EWK_s <- try({ CoxShi(ps, 1:n, pmin, pmax, 30, 2, 1)},silent=TRUE)
        if(inherits(CS2_EWK_s,"try-error")){  
          CS2B_EWK[sim]<- NA
          print("cs2 NA")
        }else{
          CS2B_EWK[sim] <- CS2_EWK_s
          # print("cs2 ok")
          #print(CS2_EWK[sim])
        }
      }
      
      #CS1_EWK[sim] <- CoxShi(ps, 1:n, pmin, pmax, 30, 1, 0)
      
      #CS2B_EWK[sim] <- CoxShi(ps, 1:n, pmin, pmax, 30, 2, 1)
      #print(paste0("Sim: ", sim," CS2: ", mean(CS2B_EWK[1:sim]), " CS1: ", mean(CS1_EWK[1:sim]) ))
      
      #Projection method
      if(parms$omit_proj[parm]==FALSE){
        print("Running Projection method:")
        
        
        
        #Estimate residuals
        tic()
        coeffs          <- get_coeffs( ts,sigma_Y=sigma_Y,numcoeffs = num_coeffs ) #estimate coefficients
        toc()
        print("Coefficients computed")
        tic()
        projection      <- compute_residual_fast(coeffs,solver) #compute_residual(coeffs,U)
        toc()
        print("Residual computed")
        resids[sim]     <- projection$residual
        true_dist[sim]  <- sqrt(sum((coeffs-null_true_coeffs)^2))
        projpoint       <- U%*%projection$alpha_opt
        proj_error[sim] <- sqrt(sum((projpoint-null_true_coeffs)^2))
        
        
        
        #Bootstrap (skip if this if the test stat is obviously not going to reject)
        
        if(((sim <10 | resids[sim]> 0.9*min(boot95[1:(sim-1)])) & parms$fastboot[parm]==FALSE)  | (sim==1 | parms$fastboot[parm]==FALSE) ){
          
          resids_boot <- rep(NA,nboots)
          for(boot in 1:nboots){
            ts_boot_pre           <- sample(ts_pre,n,replace=TRUE)
            ts_boot               <- c(abs(ts_boot_pre),-abs(ts_boot_pre))-cv*parms$shift[parm] #symmeterize and re-center
            coeffs_boot           <- get_coeffs( ts_boot,sigma_Y=sigma_Y,numcoeffs = num_coeffs ) #estimate coefficients
            #resids_boot[boot]    <- sqrt(sum((coeffs_boot-coeffs)^2 )) #conservative!
            v_boot                <- projpoint+sqrt(n)*(coeffs_boot-coeffs) #
            resids_boot[boot]     <-  compute_residual_fast(v_boot,solver,alpha_start=projection$alpha_opt)$residual /sqrt(n) #compute_residual(v_boot,Uboot)$residual #does noise get you farther away?
            #resids_boot[boot]     <-  compute_residual_fast(coeffs_boot,solver,alpha_start=projection$alpha_opt)$residual - projection$residual
            
            if((boot-1) / 10 == floor((boot-1)/10)){            print(round(c(prob_hack,resids[sim],quantile(resids_boot[1:boot],0.95),boot/nboots,sim/nsims),5))
            }
          }
          boot95[sim]     <- quantile(resids_boot,0.95) 
        }
        else{  boot95[sim]  <- mean(boot95[1:(sim-1)])}
        rejects[sim] <- 1*(resids[sim]>boot95[sim]+epsilon_U+nu_resid) 
      }
      
      
      if(sim / 1 == floor(sim/1)){
        print("Using tangent cone method")
        print(paste0("parm: ", parm, " of ", num_parameterizations,", sim: ", sim, " of ", nsims, ", prob hack: ", prob_hack, " CS1: ", mean(CS1_EWK[1:sim]<.05,na.rm=TRUE), " CS2B: ", mean(CS2B_EWK[1:sim]<.05,na.rm=TRUE), " proj: ", mean(rejects[1:sim]), " pr no adj: ",  mean(resids[1:sim] > boot95[1:sim])))
        print(paste0("epsilon: ",epsilon_U, ", eps1: ", eps1, ", boot 95: ", boot95[sim], ", teststat: ", mean(resids[1:sim]) ))
        print(paste0("rej adj proj: ",mean(resids[1:sim] > (boot95[1:sim]+eps1) ) ))
        #print(round(c(sim/nsims, parm/num_parameterizations,mean(rejects[1:sim]),mean(CS2B_EWK[1:sim]<.05,na.rm=TRUE),prob_hack),5))
      }
      toc()
      
    } # ends sim loop
    
    rejects_no_nu <- resids > boot95+epsilon_U
    
    #Record the results
    parms$epsilon_U[parm]          <- epsilon_U
    parms$nu_resid[parm]           <- nu_resid
    parms$reject_rate_eps1[parm]   <- mean(resids > boot95+eps1)
    parms$reject_rate_no_adj[parm] <- mean(resids > boot95)
    parms$reject_rate_no_nu[parm]  <- mean(rejects_no_nu)
    parms$reject_rate[parm]        <- mean(rejects)
    parms$mean_resid[parm]         <- mean(resids)
    parms$breakdown[parm]          <- mean( resids[rejects_no_nu] - boot95[rejects_no_nu]-epsilon_U )
    parms$boot95[parm]             <- mean(boot95)
    parms$eps2[parm]               <- eps2
    parms$eps1[parm]               <- eps1
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
    setwd(raw_simulations)
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    filename <- sprintf("sim_raw_%s_nu%s_coeffs%d_theta%.2f.rds", 
                        timestamp,nu, num_coeffs, theta)
    results <- list(
      resids = resids,
      boot95 = boot95,
      rejects = rejects,
      proj_error = proj_error,
      true_dist = true_dist
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
  
  toc()
  return(filename)
}

