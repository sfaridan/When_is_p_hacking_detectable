# -----------------------------------------------------------------------------------
# Drop-in replacement for compute_epsilons(), based on Lemma 6' (second-order grid bound).
#
# Same call and same output fields as the original:
#   compute_epsilons(L, nx)$epsilon_U, with
#     eps1      = grid-spacing term   (was max_diff_over_delta(L, nx));  E_2(L, delta) in Lemma 6'
#     eps2      = finite-L tail term  (was compute_eps2(L));             E_1(L)        in Lemma 6'
#     epsilon_U = eps1 + eps2                                             E(L, delta)   in Lemma 6'
#
# Assumptions, all met by create_basis_stable() and create_weights():
#   * the grid is seq(-L, L, length.out = nx) plus a far point x_inf >= 2L
#     (create_basis_stable() appends 9999; its extra point 0 can only help);
#   * weights satisfy 0 <= w_j <= 1, for any number of coefficients;
#   * sigma_Y = 1, as in the original compute_epsilons().
#
# eps1: a point mass at y is replaced by the chord between its two neighbouring grid points.
#       The linear-interpolation error is at most delta^2/8 * max|phi''| = delta^2/(8*sqrt(2*pi)).
# eps2: a point mass beyond +-L is replaced by the chord between +-L and x_inf.
# Both terms are closed forms, so no numerical integration is needed.
# -----------------------------------------------------------------------------------
compute_epsilons <- function(L, nx = 3000, x_inf = 9999) {
  if (!is.numeric(L) || length(L) != 1 || !is.finite(L) || L <= 0) stop("L must be a positive number.")
  if (!is.numeric(nx) || length(nx) != 1 || nx < 2) stop("nx must be at least 2.")
  if (x_inf < 2 * L) stop("Lemma 6' needs a far grid point x_inf >= 2L (create_basis_stable() uses 9999).")
  
  delta <- 2 * L / (nx - 1)                                    # spacing of seq(-L, L, length.out = nx)
  eps1  <- delta^2 / (8 * sqrt(2 * pi))                        # interior: chord between neighbours
  eps2  <- 3^(1/4) * exp(-1 - L^2 / 6) / (sqrt(2 * pi) * L) +  # tail: chord between L and x_inf
    exp(-2 * L^2 / 3) / sqrt(2 * pi * sqrt(3))
  
  list(eps1 = eps1, eps2 = eps2, epsilon_U = eps1 + eps2)
}

L<-7
numgrid<-1000
epsilon_U_pre <- compute_epsilons(L=L, nx = numgrid)$epsilon_U #pre-compute for speed
J0<- 20
shift<-1.96
rct_results <- run_test_weight_fast(MM_data[MM_data$method=="RCT",],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
did_results <- run_test_weight_fast(MM_data[MM_data$method=="DID",],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
rdd_results <- run_test_weight_fast(MM_data[MM_data$method=="RDD",],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
iv_results  <- run_test_weight_fast(MM_data[MM_data$method=="IV",] ,numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)

#Save results
empirical_results <- as.data.frame(do.call(rbind, list(rct_results,iv_results, did_results,rdd_results )))
empirical_results$methods <- c("RCT","IV","DID","RDD")
empirical_results<- empirical_results[,c("methods","n","articles","pval","breakdown","numcoeffs","J0")]
empirical_results

# Expected output at the application settings:
#   compute_epsilons(6.5, nx = 3000)
#   $eps1      9.370e-07
#   $eps2      2.599e-05
#   $epsilon_U 2.693e-05
# The original version gives about 3.79e-04, 2.65e-04 and 6.44e-04 at the same settings.