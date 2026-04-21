
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


setup_projection_solver_tangentcone <- function(Upre,x) {
  U <- sweep(Upre, 1, x, FUN = "-") #Upre-x
  
  nvar <- ncol(U)
  
  P_sp <- as(2 * crossprod(U), "dgCMatrix")
  A    <- rbind(Matrix(1, 1, nvar), Diagonal(nvar))
  A_sp <- as(as(A, "TsparseMatrix"), "dgCMatrix")
  l    <- c(0, rep(0, nvar))
  u    <- c(Inf, rep(Inf, nvar))
  
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

    
    if(parms$expo[parm]){
      nu_resid <- .0745/nu #2*.378/sqrt(nu)/6
    }
    else{
      integrand <- function(x) { dnorm(x) * (dt(x, df = nu) - dnorm(x))^2 }
      nu_resid  <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value) #this is Delta if we assum t~student-t(nu)
    }
    
    
    #Projection basis
    U          <- create_basis(num_coeffs,L=L,numgrid=numgrid,sigma_Y=sigma_Y)
    
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

    for(sim in 1:nsims){
      
      
      tic()
      #numsimsper <- 100
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
      
      
      if(parms$smooth_hack[parm]){ #maximization p-hacking with no threshold and h~N(0,1)
        ts_pre<- pmax(ts1,ts2,ts3,ts4,ts5,ts6,ts7,ts8 )*(1*(rands<= prob_hack))+ts1*(rands>prob_hack)  #maximization p-hacking that doesn't add any discontinuities
      }
      else{
        ts_pre <- ts1*(rands>prob_hack | abs(ts1)>cv)+(rands<= prob_hack & abs(ts1)<=cv)*pmax(ts1,ts2 )
      }
      
      ts_pre <- trunc_population(ts_pre,cv,theta)
      ts <-  ts_pre
      ps <- t_to_p_normal_approx(ts_pre)
      
      #Run tests from Elliott et al.
      if(parms$omit_EWK[parm]==FALSE){
        print("Running EWK:")
        
        
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
          # print("cs2 ok")
          #print(CS2_EWK[sim])
        }
      }
      
      #Projection method
      if(parms$omit_proj[parm]==FALSE){
        print("Running Projection method:")
        
        data        <- data.frame(t=ts,title = as.character(1:length(ts)))
        sim_results <- run_test(data,num_coeffs,sigma_Y=sigma_Y,shift=cv,L=L,numgrid=numgrid,boots=nboots,U=U)
        
        epsilon_U   <- sim_results$epsilon_U
        resids[sim] <- sim_results$resid
        boot95[sim] <- sim_results$boot95
        
        rejects[sim] <- sim_results$pval <= .05  #1*(resids[sim]>boot95[sim]+epsilon_U+nu_resid) 
        print("Pval =", sim_results$pval)
      }
      
      
      if(sim / 1 == floor(sim/1)){
        print(paste0("parm: ", parm, " of ", num_parameterizations,", sim: ", sim, " of ", nsims, ", prob hack: ", prob_hack, " CS1: ", mean(CS1_EWK[1:sim]<.05,na.rm=TRUE), " CS2B: ", mean(CS2B_EWK[1:sim]<.05,na.rm=TRUE), " proj: ", mean(rejects[1:sim])))
      }
      toc()
      
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
  
  toc()
  return(filename)
}




run_test<- function(data,numcoeffs,sigma_Y=1,shift=1.96,L=6.5,numgrid=3000,boots=150,U=NULL){
  
  
  
  #Projection basis
  if(is.null(U)){ #option to pass in U 
    U          <- create_basis(numcoeffs,L=L,numgrid=numgrid,sigma_Y=sigma_Y)
  }
  
  epsilon_U  <- compute_epsilons(L,nx=numgrid)$epsilon_U
  
  
  solver <-  setup_projection_solver(U)
  
  ts                <- c(abs(data$t),-abs(data$t))-shift #symmetrize and then shift
  coeffs_orig       <- get_coeffs( ts,sigma_Y=sigma_Y,numcoeffs=numcoeffs)
  projection        <- compute_residual_fast(coeffs_orig,solver) 
  orig_resid        <- projection$residual 
  
  #The bootstrap
  projpoint           <- U%*%projection$alpha_opt
  solver_tcone        <- setup_projection_solver_tangentcone(U,projpoint)
  resids_boot         <- rep(NA,boots)
  for (b in 1:boots){
    
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
    estar                  <- (coeffs_boot - coeffs_orig)
    resids_boot[b]         <- compute_residual_fast(estar,solver_tcone,alpha_start=warmstart)$residual 
    print(paste0("boot ", b, " of ", boots))
    print(c( orig_resid,quantile(resids_boot[1:b],0.95),quantile(resids_boot[1:b],0.90),mean(orig_resid < resids_boot[1:b]+epsilon_U)))
  }
  pval             <- mean(orig_resid < resids_boot+epsilon_U)
  breakdown        <- orig_resid - quantile(resids_boot,0.95)-epsilon_U
  
  return(list(n=length(data$t),articles=length(unique(data$title)),resid=orig_resid,epsilon_U=epsilon_U,boot95=quantile(resids_boot,0.95),pval=pval,breakdown=breakdown,projpoint =projpoint))
}





