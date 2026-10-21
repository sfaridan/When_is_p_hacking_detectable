#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_all.R"))

# Plateau weights w_j = min(1, (j/J0)^(-p)), j = 0, ..., numcoeffs - 1 (j = 0 gives weight 1)
rescue_weights <- function(numcoeffs, p = 4, J0 = (numcoeffs - 1) / 2) {
  if (J0 <= 0) stop("J0 must be positive.")
  if (p <= 0.5) stop("p must exceed 1/2 so that the weights are square-summable.")
  j <- 0:(numcoeffs - 1)
  pmin(1, (j / J0)^(-p))
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
make_rescue_worker <- function(Asum, Rg, coeffs_orig, coeffs_w, p_hat_w, w, n, sn,
                               orig_resid, A_aug, Uw, M, project) {
  force(Asum); force(Rg); force(coeffs_orig); force(coeffs_w); force(p_hat_w); force(w)
  force(n); force(sn); force(orig_resid); force(A_aug); force(Uw); force(M); force(project)
  function(sampled_articles) {
    cnt         <- tabulate(sampled_articles, nbins = nrow(Asum))
    coeffs_boot <- as.numeric(crossprod(cnt, Asum)) / sum(cnt * Rg)
    estar_w     <- w * sqrt(n) * (coeffs_boot - coeffs_orig)          # weighted bootstrap draw
    d_fs <- project(coeffs_w + estar_w * sn, A_aug, Uw, M)$residual    # FS perturbation
    d_rc <- project(p_hat_w + estar_w / sqrt(n), A_aug, Uw, M)$residual # recentered
    c(fs   = (d_fs - orig_resid / sqrt(n)) / sn,
      rc   = sqrt(n) * d_rc,
      norm = sqrt(sum(estar_w^2)))
  }
}


run_test_rescue_fast <- function(data, numcoeffs, sigma_Y = 1, shift = 1.96,
                                 L = 6.5, numgrid = 3000, boots = 150, verbose = FALSE,
                                 U = NULL, solver = NULL, solver_boot = NULL, epsilon_U_pre = NULL,
                                 p = 4, J0 = (numcoeffs - 1) / 2, cl = NULL, M = 1e5) {
  
  start_time <- Sys.time()

  n <- length(data$t)
  
  # Weights and weighted basis
  w <- rescue_weights(numcoeffs, p = p, J0 = J0)
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
  worker <- make_rescue_worker(Asum, Rg, coeffs_orig, coeffs_w, p_hat_w, w, n, sn,
                               orig_resid, A_aug, Uw, M, project_simplex_nnls)
  if (is.null(cl)) {
    out <- lapply(draws, worker)
  } else {
    out <- parallel::parLapply(cl, draws, worker)
  }
  out <- do.call(rbind, out)
  resids_boot <- out[, "fs"]
  recent_boot <- out[, "rc"]
  norms_boot  <- out[, "norm"]
  
  # 1. FS test
  boot95    <- quantile(resids_boot, 0.95)
  pval      <- mean(orig_resid < resids_boot + epsilon_U)
  breakdown <- max(c(0, (orig_resid - quantile(resids_boot, 0.95) - epsilon_U) / sqrt(n)))
  
  # 2. Recentered test
  boot95_recent    <- quantile(recent_boot, 0.95)
  pval_recent      <- mean(orig_resid < recent_boot + epsilon_U)
  breakdown_recent <- max(c(0, (orig_resid - boot95_recent - epsilon_U) / sqrt(n)))
  
  # 3. Max test (recommended)
  boot95_max    <- max(unname(boot95), unname(boot95_recent))
  pval_max      <- max(pval, pval_recent)
  breakdown_max <- max(c(0, (orig_resid - boot95_max - epsilon_U) / sqrt(n)))
  
  # 4. Norm test
  boot95_norm    <- quantile(norms_boot, 0.95)
  pval_norm      <- mean(orig_resid < norms_boot + epsilon_U)
  breakdown_norm <- max(c(0, (orig_resid - boot95_norm - epsilon_U) / sqrt(n)))
  
  q95 <- c(unname(boot95), unname(boot95_recent), boot95_max, unname(boot95_norm))
  pv  <- c(pval, pval_recent, pval_max, pval_norm)
  comparison <- data.frame(
    test       = c("FS numerical derivative", "recentered at projection",
                   "max(FS, recentered) [recommended]", "norm (conservative)"),
    resid      = orig_resid,
    boot95     = q95,
    crit_value = q95 + epsilon_U,
    pval       = pv,
    reject_05  = pv <= 0.05,
    breakdown  = c(breakdown, breakdown_recent, breakdown_max, breakdown_norm)
  )
  
  # Crude truncation bound (enters only the class of alternatives guaranteed to be detected)
  jt          <- numcoeffs:(numcoeffs + 2e5)
  tail_w2     <- sum(pmin(1, (jt / J0)^(-p))^2) + J0^(2 * p) * max(jt)^(1 - 2 * p) / (2 * p - 1)
  trunc_bound <- sqrt(tail_w2 / (2 * pi)) + 2^(-(numcoeffs - 1) / 2) / sqrt(4 * pi)
  
  # Theory: FS and recentered draws are each <= the norm draw (exact here, up to ~1e-9)
  diagnostics <- list(
    frac_fs_above_norm     = mean(resids_boot > norms_boot + 1e-8),
    frac_recent_above_norm = mean(recent_boot > norms_boot + 1e-8),
    max_recent_minus_norm  = max(recent_boot - norms_boot),
    trunc_bound_crude      = trunc_bound
  )
  
  if (verbose) {
    print(comparison)
    print(unlist(diagnostics))
    print(Sys.time() - start_time)
  }
  
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
    boot95_recent = boot95_recent,
    pval_recent = pval_recent,
    breakdown_recent = breakdown_recent,
    boot95_max = boot95_max,
    pval_max = pval_max,
    breakdown_max = breakdown_max,
    boot95_norm = boot95_norm,
    pval_norm = pval_norm,
    breakdown_norm = breakdown_norm,
    comparison = comparison,
    diagnostics = diagnostics,
    weights = w,
    p = p,
    J0 = J0,
    coeffs=coeffs,
    resids_boot = resids_boot,
    recent_boot = recent_boot,
    norms_boot = norms_boot
  ))
}

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
boots      <- 1000 #bootstrap repetitions
J          <- 60
coeffs     <- J+1


#Run the projection test with more coefficients
for (J0 in c(20,30)){
  rct_results <- run_test_rescue_fast(MM_data[MM_data$method=="RCT",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0)
  did_results <- run_test_rescue_fast(MM_data[MM_data$method=="DID",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0)
  rdd_results <- run_test_rescue_fast(MM_data[MM_data$method=="RDD",],numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0)
  iv_results  <- run_test_rescue_fast(MM_data[MM_data$method=="IV",] ,numcoeffs=coeffs,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0)
  
  #Save results
  empirical_results <- as.data.frame(do.call(rbind, list(rct_results,iv_results, did_results,rdd_results )))
  empirical_results$methods <- c("RCT","IV","DID","RDD")
  empirical_results<- empirical_results[,c("methods","n","articles","pval","breakdown","coeffs","J0")]
  setwd(paste0(root,"/results"))
  name <- paste0("MM_weights_results_",J0,"_J0")
  saveRDS(empirical_results, file = name)
  print(readRDS(name))
}
