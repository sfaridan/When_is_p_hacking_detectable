
library(osqp)
library(Matrix)
library(tictoc)

# Plateau weights w_j = min(1, (j/J0)^(-p)), j = 0, ..., numcoeffs - 1 (j = 0 gives weight 1)
create_weights <- function(numcoeffs, p = 4, J0=20) {
  if (J0 <= 0) stop("J0 must be positive.")
  if (p <= 0.5) stop("p must exceed 1/2 so that the weights are square-summable.")
  j <- 0:(numcoeffs - 1)
  pmin(1, ( 1/p)^(j-J0) )
}

# Numerically stable phi(u) psi_j(u) / sigma_Y with u = x / sigma_Y, columns j = 0, ..., numcoeffs - 1
phi_psi_stable <- function(x, numcoeffs, sigma_Y = 1) {
  u   <- x / sigma_Y
  sq  <- sqrt(dnorm(u))
  out <- matrix(0, nrow = length(u), ncol = numcoeffs)
  hm2 <- sq                                  # h_0 = psi_0 * sqrt(phi)
  out[, 1] <- hm2 * sq / sigma_Y
  if (numcoeffs >= 2) {
    hm1 <- u * sq                            # h_1 = psi_1 * sqrt(phi)
    out[, 2] <- hm1 * sq / sigma_Y
  }
  if (numcoeffs >= 3) {
    for (k in 2:(numcoeffs - 1)) {
      hk <- (u * hm1 - sqrt(k - 1) * hm2) / sqrt(k)
      out[, k + 1] <- hk * sq / sigma_Y
      hm2 <- hm1
      hm1 <- hk
    }
  }
  out
}

# Stable drop-in for get_coeff_matrix (same signature and output layout)
get_coeff_matrix_stable <- function(data, sigma_Y = 1, numcoeffs = 100) {
  phi_psi_stable(data, numcoeffs, sigma_Y)
}

# Stable drop-in for create_basis (same signature and output layout)
create_basis_stable <- function(num_coeffs, L = 6, numgrid = 1500, sigma_Y = 1) {
  pis <- c(seq(-L, L, length.out = numgrid), 0, 9999)
  s   <- sqrt(sigma_Y^2 + 1)
  t(phi_psi_stable(pis, num_coeffs, s)) * sqrt(sigma_Y^2 / (1 + sigma_Y^2))^(0:(num_coeffs - 1))
}

# TRUE if a solver's stored basis A equals the weighted basis B
same_basis <- function(A, B) {
  !is.null(A) && identical(dim(A), dim(B)) && isTRUE(max(abs(as.matrix(A) - B)) < 1e-12)
}


# Exact projection of v onto conv(columns of Uw): NNLS with the simplex constraint as a heavy row
project_simplex_nnls <- function(v, A_aug, Uw, M) {
  fit   <- nnls::nnls(A_aug, c(v, M))
  alpha <- fit$x
  r     <- v - as.numeric(Uw %*% alpha)
  list(residual = sqrt(sum(r^2)), alpha_opt = alpha)
}

# Builds the per-draw computation as a self-contained closure, so it can be shipped to workers
make_worker <- function(Asum, Rg, coeffs_orig, coeffs_w, p_hat_w, w, n, sn,
                        orig_resid, A_aug, Uw, M, project) {
  force(Asum); force(Rg); force(coeffs_orig); force(coeffs_w); force(p_hat_w); force(w)
  force(n); force(sn); force(orig_resid); force(A_aug); force(Uw); force(M); force(project)
  function(sampled_articles) {
    cnt         <- tabulate(sampled_articles, nbins = nrow(Asum))
    coeffs_boot <- as.numeric(crossprod(cnt, Asum)) / sum(cnt * Rg)
    estar_w     <- w * sqrt(n) * (coeffs_boot - coeffs_orig)          # weighted bootstrap draw
    d_fs <- project(coeffs_w + estar_w * sn, A_aug, Uw, M)$residual    # FS perturbation
    c(fs   = (d_fs - orig_resid / sqrt(n)) / sn)
  }
}


run_test_weight_fast <- function(data, numcoeffs, sigma_Y = 1, shift = 1.96,
                                 L = 6.5, numgrid = 3000, boots = 150, verbose = FALSE,
                                 U = NULL, solver = NULL, solver_boot = NULL, epsilon_U_pre = NULL,
                                 p = 4, J0 = (numcoeffs - 12), cl = NULL, M = 1e5) {
  
  start_time <- Sys.time()
  
  n <- length(data$t)
  
  # Weights and weighted basis
  w <- create_weights(numcoeffs, p = p, J0 = J0)
  if (is.null(U)) {
    U <- create_basis_stable(numcoeffs, L = L, numgrid = numgrid, sigma_Y = sigma_Y)
  } else if (any(!is.finite(U))) {
    warning("U has non-finite entries (create_basis overflows for numcoeffs >= 82); rebuilding it stably.")
    U <- create_basis_stable(numcoeffs, L = L, numgrid = numgrid, sigma_Y = sigma_Y)
  }
  Uw    <- U * w             # row j scaled by w_j
  A_aug <- rbind(Uw, M)      # last row enforces sum(alpha) = 1
  
  if (is.null(epsilon_U_pre)) {
    epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U
  }
  epsilon_U <- sqrt(n) * epsilon_U_pre   # unchanged by the weights, since w_j <= 1
  initialization_finished_time <- Sys.time()
  
  # Symmetrize and shift once
  ts <- c(abs(data$t), -abs(data$t)) - shift
  
  coeff_mat   <- get_coeff_matrix_stable(ts, sigma_Y = sigma_Y, numcoeffs = numcoeffs)
  coeffs_orig <- colMeans(coeff_mat)
  coeffs_w    <- w * coeffs_orig
  
  projection <- project_simplex_nnls(coeffs_w, A_aug, Uw, M)
  orig_resid <- sqrt(n) * projection$residual
  
  projpoint <- U %*% projection$alpha_opt                   # fitted honest t-curve, unweighted coefficients
  p_hat_w   <- as.numeric(Uw %*% projection$alpha_opt)      # the same point, weighted coordinates
  
  # Article-level sums. Articles are indexed as in run_test(): split() orders them by factor(title).
  art   <- as.integer(factor(data$title))
  m     <- max(art)
  Asum  <- rowsum(coeff_mat, c(art, art))                   # m x numcoeffs; rows i and i + n belong to art[i]
  Rg    <- 2 * tabulate(art, nbins = m)                     # rows per article
  
  # Draw all article samples first, in the same order as run_test(), then evaluate them
  draws  <- lapply(seq_len(boots), function(b) sample.int(m, size = m, replace = TRUE))
  sn     <- n^(-1/3)
  worker <- make_worker(Asum, Rg, coeffs_orig, coeffs_w, p_hat_w, w, n, sn,
                        orig_resid, A_aug, Uw, M, project_simplex_nnls)
  if (is.null(cl)) {
    out <- lapply(draws, worker)
  } else {
    out <- parallel::parLapply(cl, draws, worker)
  }
  out <- do.call(rbind, out)
  resids_boot <- out[, "fs"]
  
  # 1. FS test
  boot95    <- quantile(resids_boot, 0.95)
  pval      <- mean(orig_resid < resids_boot + epsilon_U)
  breakdown <- max(c(0, (orig_resid - quantile(resids_boot, 0.95) - epsilon_U) / sqrt(n)))
  
  
  return(list(
    n = length(data$t),
    articles = m,
    resid = orig_resid,
    epsilon_U = epsilon_U,
    boot95 = boot95,
    pval = pval,
    breakdown = breakdown,
    projpoint = projpoint,
    runtime = Sys.time() - start_time,
    init_time = initialization_finished_time - start_time,
    weights = w,
    p = p,
    J0 = J0,
    numcoeffs=numcoeffs,
    resids_boot = resids_boot
  ))
}


#takes in population without publication bias and truncates it
# if absolute value of t-score is above cv, you stay
# if |T|<cv then stay with probability theta
trunc_population <- function(untrunc_population, cv, theta){
  randkeep <- runif(length(untrunc_population))
  truncated <- untrunc_population[ abs(untrunc_population)>= cv | randkeep<=theta  ]
  return(truncated)
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
    pi0_shape  <- parms$pi0_shape[parm]
    prob_hack  <- parms$prob_hack[parm] # 2
    hack_type  <- parms$hack_type[parm]
    J0         <- parms$J0[parm]
    p          <- parms$p[parm]
    
    start_time <- Sys.time()
    
    
    
    integrand <- function(x) { dnorm(x) * (dt(x, df = nu) - dnorm(x))^2 }
    nu_resid  <- sqrt(integrate(integrand, lower = -Inf, upper = Inf)$value) #this is Delta if we assume t~student-t(nu)
    
    
    
    #Precompute to accelerate the projection test
    epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U
    
    #Set up Elliott et al pval vectors
    lcms_EWK        <- rep(NA,nsims)
    disconts_EWK    <- rep(NA,nsims)
    CS1_EWK         <- rep(NA,nsims)
    Fisher_EWK      <- rep(NA,nsims)
    CS2B_EWK        <- rep(NA,nsims)
    binomial_EWK    <- rep(NA,nsims)
    
    #set up kudrin tests
    kt_pvals        <- matrix(NA_real_, nsims, 3, dimnames = list(NULL, c("T_inf", "T_H", "T_H_hybrid")))  # Kudrin (2025) tests
    
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
      
      
      if(pi0_shape == "normal") { hnoise <- rnorm(n)*sigma_h + h_center}
      if(pi0_shape == "uniform"){ hnoise <- (runif(n)-0.5)*sigma_h  + h_center}
      if(pi0_shape == "point")  { hnoise <- rep(0,n)*sigma_h +h_center}
      if(pi0_shape == "chi2")   { hnoise <- rchisq(n,h_center)}  # ignores other paramters
      if(pi0_shape == "double_normal") { 
        hnoise <- rnorm(n)*sigma_h + h_center
        rr <- runif(n)
        hnoise[rr >= 0.5]  <- hnoise[rr >= 0.5] - 2
      }
      if(pi0_shape == "cauchy") { hnoise <- rcauchy(n)}
      if(pi0_shape == "poisson") { hnoise <- rpois(n,h_center)}
      if(pi0_shape == "kudrin_gamma") { hnoise <- rgamma(n, shape = 1.4134, rate = 1.2980) }
      if(pi0_shape == "normal_11") { hnoise <- rnorm(n)+1 }
      
      
      if(pi0_shape == "poisson(2)") { hnoise <- rpois(n,2)}
      if(pi0_shape == "chi2(2)")   { hnoise <- rchisq(n,2)}  # ignores other paramters
      if(pi0_shape == "null")   { hnoise <- rep(0,n)}
      if(pi0_shape == "onehalf")   { hnoise <- rep(0.5,n)}
      if(pi0_shape == "two")   { hnoise <- rep(2,n)}
      
      if(pi0_shape == "fitted")   { #Fitted by zcurve from Faridani (2026)
        vv <- c(0, 1, 2, 3, 4, 5, 6)
        p_vv  <- c(0.02, 0.47, 0.01, 0.27, 0.14, 0.00, 0.09)
        hnoise <-   sample(vv, size = n, replace = TRUE, prob = p_vv)
      }
      if(pi0_shape == "unif_half"){ hnoise <- runif(n)-0.5  }
      if(pi0_shape == "standard_normal") { hnoise <- rnorm(n)}
      if(pi0_shape == "hard") { hnoise <- rnorm(n)*.7+2 }
      
      if(pi0_shape == "smooth_pub") { hnoise <- rnorm(n)*sqrt(smooth_pub_par(parms, parm)$tau2) } # smooth_pub: latent effects (printed only)
      
      
      hs      <- hnoise #h_center + hnoise*sigma_h # latent true effects
      
      #Draw individual t-scores
      ts1     <- rt(n,nu)+hs
      ts2     <- rt(n,nu)+hs
      ts3     <- rt(n,nu)+hs
      
      if(hack_type == "max"){ #maximization p-hacking with no threshold 
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
      if(startsWith(as.character(hack_type), "kudrin")){
        ts_pre <- kudrin_tstats(hs, hacker = (rands <= prob_hack), hack_type = hack_type, parms_row = parms[parm, ])
      }
      
      if(theta != 1){  #publication bias if theta < 1
        ts_pre <- trunc_population(ts_pre,cv,theta) 
      }
      
      if(pi0_shape == "smooth_pub"){ # smooth_pub DGP (smooth_pub_dgp.R): exactly n published t-stats; replaces hacking/truncation above
        sp     <- smooth_pub_par(parms, parm)
        ts_pre <- r_smooth_pub(n, tau2 = sp$tau2, sw2 = sp$sw2, theta = sp$theta, nu = nu)
        print(paste0("smooth_pub: tau2=", sp$tau2, ", sw2=", sp$sw2, ", theta=", sp$theta,
                     ", P(R)=", round(smooth_pub_prob(sp$tau2, sp$sw2, sp$theta), 3),
                     ", detectable=", sp$theta > smooth_pub_threshold(sp$tau2, sp$sw2)))
      }
      
      ts <-  ts_pre # finalize observed t-scores
      
      print(paste0("H: ", pi0_shape,"frac hacking: ", prob_hack,", hacktype: ", hack_type))
      #print(paste0("frac hacking: ", prob_hack,", hacktype: ", hack_type))
      #print(summary(ts))
      
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
      
      #Kudrin (2025) tests T_inf, T_H, T_H_hybrid; requires kudrin_tests.R
      if(isTRUE(parms[["kudrin_tests"]][parm])){
        kt_start <- Sys.time()
        kt_s <- try(kudrin_tests(ts_pre, boots = nboots)$pval, silent = TRUE)
        if(inherits(kt_s, "try-error")){ print(paste0("Kudrin tests NA: ", kt_s)) }else{ kt_pvals[sim, ] <- kt_s }
        print(paste0("Kudrin tests time: ", Sys.time() - kt_start))
      }
      
      #Projection method
      if(parms$omit_proj[parm]==FALSE){
        proj_start  <- Sys.time()
        
        data        <- data.frame(t=ts,title = as.character(1:length(ts)))
        sim_results <- run_test_weight_fast(data,num_coeffs,sigma_Y=sigma_Y,shift=cv,
                                 L=L,numgrid=numgrid,boots=nboots,
                                 J0=J0,p=p,epsilon_U_pre=epsilon_U_pre)
        
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
        if(isTRUE(parms[["kudrin_tests"]][parm])){   # Kudrin (2025) tests
          for(nm in colnames(kt_pvals)) print(paste0(nm,": ", round(mean(kt_pvals[, nm] <= .05, na.rm = TRUE),3)))
        }
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
    
    #Record Kudrin test results
    if(isTRUE(parms[["kudrin_tests"]][parm])){   # Kudrin (2025) tests
      for(nm in colnames(kt_pvals)) parms[parm, paste0("Kudrin_", nm)] <- mean(kt_pvals[, nm] <= .05, na.rm = TRUE)
    }
    
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




