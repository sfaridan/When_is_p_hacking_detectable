################################################################################
#  simulations_brodeur_rct.R
#
#  Power curves of the projection test on meta-samples built to look like the RCTs (or any
#  subsample: CFG$method; "ALL" pools all four methods) of Brodeur,
#  Cook and Heyes (2020) as application_MM_all_four_rounding.R analyzes them. Every simulated
#  meta-sample has the application's subsample (Kudrin's drop of mu = 0, no standard error
#  printed as zero, kappa = a / s_min <= 0.5): the same t-statistics, the same articles with the
#  same number of t-statistics each, the same printed precision of every coefficient, standard
#  error and reported t, effects fitted to their t-curve, and their within-article dependence.
#  Each meta-sample is p-hacked to extent tau by one of two methods of run_sims() and
#  de-rounded by each of the application's four options: 2 x 4 = 8 power curves, whose
#  leftmost points (tau = 0) are sizes.
#
#  The DGP, for t-statistic i of article g(i); everything but the draws is fixed by the data:
#    Effect    H_i = Q(Phi(sqrt(rho) A_g + sqrt(1 - rho) E_i)), A_g, E_i ~ N(0, 1). Q is the
#              quantile function of Pi-hat, by default a Gamma, Kudrin's family, fitted by
#              maximum likelihood to the subsample's |t| (de-rounded as CFG$fit_option) under the
#              honest model |T| = |H + Z| (CFG$pi_fit), so H_i ~ Pi-hat; rho correlates
#              effects within an article and is set so that the application's design effect
#              deff() of the projection test's features matches the data (CFG$deff_option).
#              Effects do not depend on how a result is printed.
#    Honest    T1 = H + Z1, with Z_k ~ N(0, 1) (Student-t(nu) if CFG$nu is finite) independent;
#              with corr_in = "noise" the dependence is in Z instead: Z_k = sqrt(rho) W_gk +
#              sqrt(1 - rho) e_ik, still N(0, 1) for each result, and effects are independent
#    Hacked    with probability tau (each result, or each article: CFG$hack_level), as run_sims():
#                "max"        max(H + Z1, ..., H + Z8)
#                "threshold"  T1 if T1 > 1.96, else max(T1, H + Z2)
#    Printed   A coefficient/SE result keeps the data's printed standard error and its table's
#              coefficient unit a_tab. With u = min(b, a_tab), the finest unit its table prints
#              (so the trailing zeros the data lost come back, as Kranz and Puetz restore them),
#              the true se is x ~ U[s/u - 1/2, s/u + 1/2) in units of u and the coefficient T x
#              is printed to a_tab. A reported t-statistic keeps its printed decimals.
#    Analyst   The application's units and options: b read from the stored se, a = max(a_tab, b),
#              s_min = s - b/2, s_max = s + b/2, and
#                "none"        m / s                     "smin"        m / s_min
#                "uniform"     (m + V) / W               "smin_noise"  (m + V) / s_min
#              with V ~ U(-a/2, a/2) the same in "uniform" and "smin_noise", W ~ U(s_min, s_max);
#              a reported t is used as printed, plus V in its units under the two noisy options.
#              The true se is at least s_min, so the stretch se / s_min of Proposition 4 is >= 1.
#
#  Tests: your run_test_weight_fast() with the application's settings (J = 40, J0 = 20, p = 2,
#  clustered by article); the application's run_ewk() (EKW tests, each reported separately);
#  the p-curve test for left skew that Simonsohn's report asks for (point 5: each p < .05 times
#  20, normal quantiles, combined by Stouffer's method); optionally your kudrin_tests().
#  Under "smin" and "smin_noise" the projection test rejects only when its breakdown statistic
#  B-hat = max(0, (resid - cv - epsilon_U) / sqrt(n)) exceeds the rounding allowance C mean(kappa^6),
#  C = 0.13 for J0 = 20 and 0.31 for J0 = 30 (Appendix D.1, step 3; Propositions 3 and 4), as the
#  paper's tables do; under "none" and "uniform", which have no such bound, when p <= alpha.
#  "full" mode bootstraps every test (proj_boots draws); "warp" draws once per test and pools
#  the critical value over repetitions (Giacomini, Politis and White 2013), about 100x faster.
#  The application averages p-values over 200 de-roundings for "uniform" and "smin_noise";
#  here each repetition uses CFG$noise_draws of them (1 by default; more needs "full").
#
#  Outputs in OUT_DIR (NAME = CFG$name, by default brodeur_<method>_power: brodeur_rct_power,
#  brodeur_all_power, ..., so pooled and single-method runs keep separate caches):
#    NAME_proj.pdf/.png       the 8 power curves of the projection test (one panel per method)
#    NAME_<test>.pdf/.png     the same for the Stouffer test, the smallest EKW p-value, T_inf, T_H
#                             and T_H_hybrid (those that are run)
#    NAME_<test>.tex          their tables; NAME_EKW.tex reports every EKW test separately
#    NAME_rates.csv           all rejection rates by method, option and tau
#    NAME_dw.csv / .tex       exact distortion d_w of each simulated form of p-hacking at each tau,
#                             and the tau at which it equals each B-hat of the application
#    NAME_calibration.txt/.rds  the data, Pi-hat, rho and design-effect fit
#    reps/rep_NNNNN.rds       one file per repetition, so an interrupted run resumes
#
#  Run once with QUICK <- TRUE (about an hour on 15 cores) to check that everything works.
################################################################################

rm(list = ls())

## ---- 0. Settings -------------------------------------------------------------------------

ROOT      <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
CODE_DIR  <- file.path(ROOT, "code")
DATA_FILE <- "C:/Users/sfaridani6/Documents/Data/120246-V1/MM data.dta"   # as in the application

QUICK <- FALSE

CFG <- list(
  name          = NULL,                  # NULL: brodeur_<method>_power
  method        = "ALL",                 # "RCT", "IV", "DID", "RDD", several (c("DID", "IV")), or
  # "ALL": all four pooled (one Pi-hat, one rho, clustered
  # by article across methods)
  kappa_cap     = 0.5,                   # as in the application script (the paper says 0.3)
  se_rounded_up = FALSE,                 # as in the application
  
  # p-hacking
  taus          = seq(0, 1, by = 0.1),   # share of p-hacked results; tau = 0 gives size
  hack_methods  = c("max", "threshold"),
  max_draws     = 8,                     # "max": best of 8
  cv            = 1.96,                  # "threshold": second draw when T1 <= cv
  hack_level    = "result",              # "result" or "article" (an article hacks all or none)
  nu            = Inf,                   # noise: Inf = normal, else Student-t degrees of freedom
  
  # De-rounding: the application's options 1-4
  options       = c("none", "uniform", "smin", "smin_noise"),
  noise_draws   = 1,                     # de-roundings averaged for "uniform" and "smin_noise"
  # (p-values; application: 200; > 1 needs "full")
  
  # Calibration
  fit_option    = "smin",                # t-statistics Pi-hat is fitted to
  deff_option   = "smin",                # t-statistics whose design effect rho matches
  corr_in       = "effects",             # what carries the within-article dependence: "effects"
  # (a copula on H) or "noise" (Z equicorrelated within an
  # article, still N(0, 1) for each result)
  n_cal         = 20,                    # simulated meta-samples per candidate rho
  rho_max       = 0.99,
  pi_fit        = "kudrin_scale",               # effects: "gamma" (maximum likelihood, Kudrin's family),
  # "npmle", or Kudrin's published Gamma(1.4134, 1.2980)
  # read with 1.2980 a rate ("kudrin_rate", mean 1.09, as
  # his footnote 13) or a scale ("kudrin_scale", mean 1.83,
  # as kudrin_pi_hat in the companion scripts); with
  # these two nothing is estimated
  compare_kudrin = TRUE,                 # report the log-likelihood of Kudrin's two readings on
  # these |t| (an evaluation for the report, not a fit)
  npmle_iters   = 3000,
  fit_cap       = 20,                    # in the effect fit, |t| above this is censored: it only
  # says "very large" (the data hold |t| in the millions)
  se_cell       = "uniform",             # true se within its printed cell: "uniform" or "benford"
  
  # Projection test, with the application's settings
  num_coeffs    = 41,                    # J = 40
  J0            = 20,
  p_weight      = 2,
  shift         = 1.96,
  L             = 6.5,
  numgrid       = 3000,
  proj_mode     = "full",                # "full" or "warp"
  proj_boots    = 200,                   # application: 1000
  epsilon_U_pre = NULL,                  # NULL: compute_epsilons() once (slow); paper: 0.00063
  allowance_C   = NULL,                  # rounding allowance constant; NULL: 0.13 (J0 = 20) or 0.31
  # (J0 = 30), the paper's values for p = 2
  
  # Other tests (each adds its own run time)
  run_kudrin    = TRUE,
  kudrin_boots  = 199,
  run_ekw       = TRUE,                  # each EKW test separately (adds a few s per test call)
  run_stouffer  = TRUE,                  # p-curve test for left skew (Simonsohn); instantaneous
  table_taus    = NULL,                  # taus in the tables (NULL: all), e.g. c(0, .1, .25, .5, 1)
  
  # Exact distortions d_w (no sampling), to read the application's B-hat (needs the nnls package)
  dw_taus       = seq(0, 1, by = 0.1),   # rows of the d_w table
  bhat          = c(RCT = 2.8e-3, IV = 3.1e-3, DID = 0.9e-3, RDD = 1.3e-3),   # B-hat by design in
  # the application (Table 3 of the 9/20/2026 draft; update
  # after rerunning it), or NULL
  
  # Computing
  reps          = 500,
  alpha         = 0.05,
  n_cores       = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  seed          = 20261029,
  resume        = TRUE
)

# The subsample, its label in figures and tables, and the output folder
METHODS_ALL <- c("RCT", "IV", "DID", "RDD")
POOL_ALL    <- length(CFG$method) == 1 && toupper(CFG$method) == "ALL"
SUBSAMPLE   <- if (POOL_ALL) METHODS_ALL else CFG$method
SUB_TAG     <- if (POOL_ALL) "all" else tolower(paste(SUBSAMPLE, collapse = "_"))
SUB_LABEL   <- if (POOL_ALL) "all four methods pooled" else if (identical(SUBSAMPLE, "RCT")) "RCTs" else
  paste(SUBSAMPLE, collapse = " + ")
if (is.null(CFG$name)) CFG$name <- sprintf("brodeur_%s_power", SUB_TAG)
OUT_DIR <- file.path(ROOT, "simulations", CFG$name)
if (QUICK) {
  CFG$reps          <- 40
  CFG$taus          <- c(0, 0.25, 0.5, 0.75, 1)
  CFG$proj_boots    <- 99
  CFG$kudrin_boots  <- 99
  CFG$n_cal         <- 6
  CFG$epsilon_U_pre <- 0.00063
  OUT_DIR           <- paste0(OUT_DIR, "_quick")
}
dir.create(file.path(OUT_DIR, "reps"), recursive = TRUE, showWarnings = FALSE)

stopifnot(all(CFG$options %in% c("none", "uniform", "smin", "smin_noise")),
          all(CFG$hack_methods %in% c("max", "threshold")), CFG$max_draws >= 2,
          CFG$hack_level %in% c("result", "article"), CFG$proj_mode %in% c("full", "warp"),
          CFG$se_cell %in% c("uniform", "benford"), CFG$noise_draws >= 1,
          CFG$corr_in %in% c("effects", "noise"),
          CFG$pi_fit %in% c("gamma", "npmle", "kudrin_rate", "kudrin_scale"),
          CFG$fit_option %in% c("none", "uniform", "smin", "smin_noise"),
          CFG$deff_option %in% c("none", "uniform", "smin", "smin_noise"),
          any(CFG$taus == 0), all(CFG$taus >= 0 & CFG$taus <= 1))
if (CFG$noise_draws > 1 && CFG$proj_mode == "warp") stop("noise_draws > 1 needs proj_mode = \"full\".")


## ---- 1. Your code ------------------------------------------------------------------------

CODE_FILES <- file.path(CODE_DIR, c("functions_all_weighted.R",
                                    if (CFG$run_ekw) "Tests.R",
                                    if (CFG$run_kudrin) "kudrin_tests.R"))
missing_files <- CODE_FILES[!file.exists(CODE_FILES)]
if (length(missing_files) > 0) stop("Not found: ", paste(missing_files, collapse = ", "))
for (f in CODE_FILES) source(f)
stopifnot(exists("run_test_weight_fast"), exists("create_basis_stable"), exists("compute_epsilons"))

# EKW tests, called as in the application's run_ewk(); each is guarded, so one failure leaves
# the others (and min is the smallest p-value among those that ran)
run_ewk <- function(ps, ids = 1:length(ps), pmax = 0.15) {
  pmin <- min(ps[ps > 0])
  g    <- function(expr) tryCatch(as.numeric(expr)[1], error = function(e) NA_real_)
  out  <- list(lcms_EWK     = g(LCM(ps, pmin, pmax)),
               Fisher_EWK   = g(Fisher(ps, pmin, pmax)),
               disconts_EWK = g(Discontinuity_test(ps, .05)),
               binomial_EWK = g(Binomial(ps, .04, .05, "c")),
               CS1_EWK_s    = g(CoxShi(ps, ids, pmin, pmax, 15, 1, 0)),
               CS2_EWK_s    = g(CoxShi(ps, ids, pmin, pmax, 15, 2, 1)))
  v       <- unlist(out)
  out$min <- if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
  out
}

# Result columns of the EKW tests, and the run_ewk() element each comes from
EKW_TESTS <- c(EKW_LCM = "lcms_EWK", EKW_Fisher = "Fisher_EWK", EKW_Discontinuity = "disconts_EWK",
               EKW_Binomial = "binomial_EWK", EKW_CS1 = "CS1_EWK_s", EKW_CS2B = "CS2_EWK_s",
               EKW_min = "min")


## ---- 2. The Brodeur et al. subsample, as the application prepares it ----------------------

n_dec <- function(x, max_dec = 12) {     # decimals as printed (the application's n_dec())
  x   <- signif(x, 6)
  out <- rep(NA_integer_, length(x))
  for (k in 0:max_dec) {
    y   <- x * 10^k
    hit <- is.na(out) & is.finite(x) & abs(y - round(y)) <= 1e-9 * pmax(1, abs(y))
    out[hit] <- k
  }
  out[is.na(out) & is.finite(x)] <- max_dec
  out
}

plain <- function(v) {
  if (inherits(v, "haven_labelled")) v <- haven::zap_labels(v)
  as.vector(v)
}

# One row per t-statistic of the subsample, with everything the DGP and the analyst need. Unit
# columns are in units of u = min(b, a_tab) for coefficient/SE results: S the printed se, c_rel
# the coefficient's printing unit, a_rel and b_rel the analyst's a and b, smin_rel and smax_rel.
build_design <- function() {
  if (!file.exists(DATA_FILE)) stop("Brodeur et al. data not found: ", DATA_FILE)
  raw  <- haven::read_dta(DATA_FILE)
  miss <- setdiff(c("mu", "sd", "t", "title", "table", "method"), names(raw))
  if (length(miss) > 0) stop("Columns not in the data: ", paste(miss, collapse = ", "))
  d <- data.frame(mu = as.numeric(plain(raw$mu)), sd = as.numeric(plain(raw$sd)),
                  t = as.numeric(plain(raw$t)), title = as.character(plain(raw$title)),
                  table = as.character(plain(raw$table)), method = as.character(plain(raw$method)),
                  stringsAsFactors = FALSE)
  d <- d[!(is.na(d$mu) | d$mu == 0), ]                    # Kudrin's (2025) drop
  dec_mu <- n_dec(d$mu)
  dec_sd <- n_dec(d$sd)
  dec_t  <- n_dec(d$t)
  mu_p   <- round(d$mu, dec_mu)
  sd_p   <- round(d$sd, dec_sd)
  t_p    <- round(d$t,  dec_t)
  tab_id <- paste(d$title, d$table, sep = " | ")
  b_sd   <- 10^(-dec_sd)
  a_tab  <- 10^(-ave(dec_mu, tab_id, FUN = function(x) if (all(is.na(x))) NA_integer_ else max(x, na.rm = TRUE)))
  a_mu   <- pmax(a_tab, b_sd)
  s_min  <- sd_p - if (CFG$se_rounded_up) b_sd else b_sd / 2
  s_max  <- sd_p + if (CFG$se_rounded_up) 0 else b_sd / 2
  a_t    <- 10^(-dec_t)
  has_ratio  <- !is.na(mu_p) & !is.na(sd_p)
  drop_se0   <- has_ratio & sd_p <= 0
  drop_kappa <- has_ratio & !drop_se0 & a_mu / s_min > CFG$kappa_cap
  use_ratio  <- has_ratio & !drop_se0 & !drop_kappa
  direct     <- !has_ratio & !is.na(t_p)
  in_method  <- d$method %in% SUBSAMPLE
  if (!any(in_method)) {
    stop("No t-statistics with method ", paste(SUBSAMPLE, collapse = ", "), "; the data have: ",
         paste(sort(unique(d$method)), collapse = ", "), ". Use one of them, several, or \"ALL\".")
  }
  message(sprintf("%s: %d t-statistics after Kudrin's drop; %d dropped (se printed as zero), %d dropped (kappa > %s), %d with neither a coefficient/SE nor a printed t.",
                  SUB_LABEL, sum(in_method), sum(in_method & drop_se0),
                  sum(in_method & drop_kappa), format(CFG$kappa_cap),
                  sum(in_method & !has_ratio & !direct)))
  k  <- which(in_method & (use_ratio | direct))
  r  <- use_ratio[k]
  u  <- ifelse(r, pmin(b_sd[k], a_tab[k]), NA)
  ds <- data.frame(title = d$title[k], gid = as.integer(factor(d$title[k])), method = d$method[k], ratio = r,
                   u_se = u, S = round(sd_p[k] / u), c_rel = round(a_tab[k] / u),
                   a_rel = round(a_mu[k] / u), b_rel = round(b_sd[k] / u),
                   a_t = ifelse(r, NA, a_t[k]), kappa = ifelse(r, a_mu[k] / s_min[k], a_t[k]),
                   mu_p = mu_p[k], sd_p = sd_p[k], s_min = s_min[k], s_max = s_max[k],
                   a_mu = a_mu[k], t_p = t_p[k], stringsAsFactors = FALSE)
  ds$smin_rel <- ds$S - if (CFG$se_rounded_up) ds$b_rel else ds$b_rel / 2
  ds$smax_rel <- ds$S + if (CFG$se_rounded_up) 0 else ds$b_rel / 2
  ds
}

# The application's deround_draw() on the data (u1, u2: uniforms on (-1/2, 1/2))
data_t <- function(option, u1 = 0, u2 = 0) {
  D     <- DESIGN
  noisy <- option %in% c("uniform", "smin_noise")
  v     <- if (noisy) u1 else 0
  num   <- D$mu_p + D$a_mu * v
  den   <- switch(option,
                  none       = D$sd_p,
                  uniform    = (D$s_min + D$s_max) / 2 + (D$s_max - D$s_min) * u2,
                  smin       = D$s_min,
                  smin_noise = D$s_min)
  ifelse(D$ratio, num / den, D$t_p + D$a_t * v)
}


## ---- 3. The simulated tables and the analyst -----------------------------------------------

draw_noise   <- function(n) if (is.finite(CFG$nu)) rt(n, CFG$nu) else rnorm(n)
draw_effects <- function(U) {
  if (identical(PI_HAT$type, "gamma")) {
    return(qgamma(pmin(U, 1 - 1e-12), shape = PI_HAT$shape, rate = PI_HAT$rate))
  }
  PI_HAT$h[pmin(findInterval(U, PI_HAT$cw) + 1L, length(PI_HAT$h))]
}

# One simulated world: effects, K noise draws (each correlated within articles through its own
# article shocks when corr_in = "noise"), the true se of every coefficient/SE result,
# the hacking uniforms and the analyst's de-rounding noise. Consumes the same random numbers
# whatever rho is, so calibration runs on common random numbers.
draw_world <- function(rho, K) {
  n  <- nrow(DESIGN)
  A  <- rnorm(N_ART)
  E  <- rnorm(n)
  W  <- matrix(rnorm(N_ART * K), N_ART, K)               # article shocks, used if corr_in = "noise"
  ce <- identical(CFG$corr_in, "effects")
  H  <- draw_effects(pnorm(if (ce) sqrt(rho) * A[DESIGN$gid] + sqrt(1 - rho) * E else E))
  Z  <- matrix(draw_noise(n * K), n, K)
  if (!ce) Z <- sqrt(rho) * W[DESIGN$gid, , drop = FALSE] + sqrt(1 - rho) * Z
  v  <- runif(n)
  lo <- DESIGN$S - if (CFG$se_rounded_up) 1 else 0.5      # the true se's cell, in units of u
  x  <- if (CFG$se_cell == "benford") lo * ((lo + 1) / lo)^v else lo + v
  hu <- if (CFG$hack_level == "article") runif(N_ART)[DESIGN$gid] else runif(n)
  list(H = H, Z = Z, T1 = H + Z[, 1], x = x, hu = hu,
       u1 = runif(n, -0.5, 0.5), u2 = runif(n, -0.5, 0.5))
}

# What the tables print for true t-statistics tstat: the coefficient in units of u, and the t
print_tables <- function(tstat, x) {
  list(m  = DESIGN$c_rel * round(tstat * x / DESIGN$c_rel),
       tp = DESIGN$a_t * round(tstat / DESIGN$a_t))
}

# The analyst's t-statistics under an option (as data_t(), in units of u)
analyst_t <- function(option, tabs, u1, u2) {
  D     <- DESIGN
  noisy <- option %in% c("uniform", "smin_noise")
  v     <- if (noisy) u1 else 0
  num   <- tabs$m + D$a_rel * v
  den   <- switch(option,
                  none       = D$S,
                  uniform    = (D$smin_rel + D$smax_rel) / 2 + (D$smax_rel - D$smin_rel) * u2,
                  smin       = D$smin_rel,
                  smin_noise = D$smin_rel,
                  stop("Unknown option: ", option))
  ifelse(D$ratio, num / den, tabs$tp + D$a_t * v)
}


## ---- 4. Calibration: Pi-hat and rho --------------------------------------------------------

# The effect fits share a grid for H and a likelihood in which |t| above CFG$fit_cap is censored:
# it contributes P(|T| > cap). The grid stops at cap + 10, its last bin running to Inf, so its
# size no longer depends on the largest |t|.
effect_grid <- function() {
  edges <- unique(c(seq(0, 10, by = 0.02), seq(10.1, CFG$fit_cap + 10, by = 0.1)))
  list(edges = edges, mids = c((edges[-length(edges)] + edges[-1]) / 2, edges[length(edges)]))
}

# Densities of the uncensored |t| (rows) for each grid point (columns), built column by column,
# and P(|h + Z| > cap) for each grid point
fit_kernel <- function(absT, mids) {
  absT <- absT[!is.na(absT)]
  dens <- if (is.finite(CFG$nu)) function(z) dt(z, CFG$nu) else dnorm
  surv <- if (is.finite(CFG$nu)) function(z) pt(z, CFG$nu, lower.tail = FALSE) else
    function(z) pnorm(z, lower.tail = FALSE)
  cen <- absT > CFG$fit_cap
  y   <- absT[!cen]
  K   <- matrix(vapply(mids, function(h) dens(y - h) + dens(y + h), numeric(length(y))), nrow = length(y))
  list(K = K, cens = surv(CFG$fit_cap - mids) + surv(CFG$fit_cap + mids), n_cens = sum(cen),
       n = length(absT))
}

# Log-likelihood of grid weights w
fit_loglik <- function(fk, w) {
  sum(log(pmax(as.vector(fk$K %*% w), 1e-300))) + fk$n_cens * log(max(sum(fk$cens * w), 1e-300))
}

# Gamma(shape, rate) mass in each bin of the effect grid
gamma_bins <- function(shape, rate, edges) {
  cdf <- pgamma(edges, shape = shape, rate = rate)
  c(diff(cdf), 1 - cdf[length(cdf)])
}

# A Gamma(shape, rate) used as given (Kudrin's values): nothing is estimated and the data are not
# touched. sig_share: P(|T| > 1.96) it implies with normal noise.
gamma_fixed <- function(shape, rate) {
  g <- effect_grid()
  w <- gamma_bins(shape, rate, g$edges)
  list(type = "gamma", shape = shape, rate = rate, loglik = NA_real_, n_censored = NA_integer_,
       sig_share = sum(w * (pnorm(-1.96 - g$mids) + pnorm(g$mids - 1.96))))
}

# Maximum likelihood Gamma(shape, rate) for the effects from |t|, under |T| = |H + Z| with Z
# normal (or Student-t(nu)): the likelihood integrates over H on the grid of gamma-CDF bins
fit_gamma <- function(absT) {
  g   <- effect_grid()
  fk  <- fit_kernel(absT, g$mids)
  nll <- function(par) -fit_loglik(fk, gamma_bins(exp(par[1]), exp(par[2]), g$edges))
  fit <- optim(log(c(1.4134, 1.2980)), nll, control = list(reltol = 1e-12, maxit = 4000))
  fit <- optim(fit$par, nll, method = "BFGS", control = list(reltol = 1e-12))
  out <- gamma_fixed(exp(fit$par[1]), exp(fit$par[2]))
  out$loglik     <- -fit$value
  out$n_censored <- fk$n_cens
  out
}

# Log-likelihoods of given gammas (a named list of c(shape, rate)) on |t|, with one kernel:
# an evaluation for the calibration report, not a fit
gamma_logliks <- function(absT, pars) {
  g  <- effect_grid()
  fk <- fit_kernel(absT, g$mids)
  vapply(pars, function(p) fit_loglik(fk, gamma_bins(p[1], p[2], g$edges)), 0)
}

# NPMLE of the effect distribution from |t|: EM for the weights of the grid's atoms, with the
# same censored likelihood as fit_gamma()
fit_pi_hat <- function(absT) {
  h  <- effect_grid()$mids
  fk <- fit_kernel(absT, h)
  w  <- rep(1 / length(h), length(h))
  for (it in seq_len(CFG$npmle_iters)) {
    post <- as.vector(crossprod(fk$K, 1 / as.vector(fk$K %*% w))) + fk$n_cens * fk$cens / sum(fk$cens * w)
    wn   <- w * post / fk$n
    done <- max(abs(wn - w)) < 1e-10
    w    <- wn
    if (done) break
  }
  loglik <- fit_loglik(fk, w)
  keep   <- w > 1e-9
  h      <- h[keep]
  w      <- w[keep] / sum(w[keep])
  list(type = "npmle", h = h, w = w, cw = c(cumsum(w)[-length(w)], 1), loglik = loglik, iters = it,
       n_censored = fk$n_cens, sig_share = sum(w * (pnorm(-1.96 - h) + pnorm(h - 1.96))))
}

# The application's design effect: paper-level clustering of the projection test's weighted
# Hermite features (trace of the cluster-sum over the plain sums of squares)
hf <- function(x, J) {                     # psi_j(x) * phi(x), j = 0..J
  G <- matrix(0, length(x), J + 1)
  G[, 1] <- dnorm(x)
  G[, 2] <- x * G[, 1]
  for (j in 1:(J - 1)) G[, j + 2] <- (x * G[, j + 1] - sqrt(j) * G[, j]) / sqrt(j + 1)
  G
}
DEFF_W <- pmin(1, CFG$p_weight^(CFG$J0 - (0:(CFG$num_coeffs - 1))))
deff_trace <- function(t, gid) {
  t <- ifelse(is.finite(t), t, 1e6)        # features vanish beyond |t| ~ 40 anyway
  X <- sweep(0.5 * (hf(t - CFG$shift, CFG$num_coeffs - 1) + hf(-t - CFG$shift, CFG$num_coeffs - 1)),
             2, DEFF_W, `*`)
  R <- sweep(X, 2, colMeans(X))
  S <- rowsum(R, gid)
  sum(S^2) / sum(R^2)
}

# rho such that honest meta-samples, de-rounded as CFG$deff_option, have the data's design effect
calibrate_rho <- function(target) {
  f <- function(rho) {
    set.seed(CFG$seed + 31)                # common random numbers across candidates
    mean(replicate(CFG$n_cal, {
      w <- draw_world(rho, 1)
      deff_trace(analyst_t(CFG$deff_option, print_tables(w$T1, w$x), w$u1, w$u2), DESIGN$gid)
    }))
  }
  lo <- 0
  hi <- CFG$rho_max
  f0 <- f(lo)
  f1 <- f(hi)
  if (target <= f0) return(list(rho = 0, deff = f0, deff_0 = f0, deff_max = f1, target = target))
  if (target >= f1) {
    warning(sprintf("The data's design effect %.2f exceeds %.2f, the most that rho = %s gives with corr_in = \"%s\"; try corr_in = \"%s\".",
                    target, f1, format(hi), CFG$corr_in, if (CFG$corr_in == "effects") "noise" else "effects"))
    return(list(rho = hi, deff = f1, deff_0 = f0, deff_max = f1, target = target))
  }
  for (it in 1:14) {
    mid <- (lo + hi) / 2
    if (f(mid) < target) lo <- mid else hi <- mid
  }
  rho <- (lo + hi) / 2
  list(rho = rho, deff = f(rho), deff_0 = f0, deff_max = f1, target = target)
}


## ---- 5. One repetition ---------------------------------------------------------------------

proj_test <- function(t) {
  pr <- tryCatch(run_test_weight_fast(data.frame(t = t, title = DESIGN$title, stringsAsFactors = FALSE),
                                      numcoeffs = CFG$num_coeffs, sigma_Y = 1, shift = CFG$shift,
                                      L = CFG$L, numgrid = CFG$numgrid,
                                      boots = if (CFG$proj_mode == "warp") 1 else CFG$proj_boots,
                                      U = PROJ$U, epsilon_U_pre = PROJ$epsilon_U_pre,
                                      p = CFG$p_weight, J0 = CFG$J0),
                 error = function(e) NULL)
  one <- function(v) if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)[1]
  if (is.null(pr)) {
    return(c(pval = NA_real_, resid = NA_real_, boot1 = NA_real_, epsU = NA_real_, breakdown = NA_real_))
  }
  c(pval = one(pr$pval), resid = one(pr$resid), boot1 = one(pr$resids_boot),
    epsU = one(pr$epsilon_U), breakdown = one(pr$breakdown))
}

kudrin_p <- function(t) {
  nm <- c("T_inf", "T_H", "T_H_hybrid")
  kt <- tryCatch(kudrin_tests(t, cluster = DESIGN$title, boots = CFG$kudrin_boots),
                 error = function(e) NULL)
  if (is.null(kt)) return(setNames(rep(NA_real_, 3), nm))
  vapply(nm, function(k) if (k %in% names(kt$pval)) as.numeric(kt$pval[[k]])[1] else NA_real_, 0)
}

ekw_p <- function(t) {
  r <- tryCatch(run_ewk(2 * pnorm(-abs(t)), ids = DESIGN$gid), error = function(e) NULL)
  vapply(EKW_TESTS, function(k) if (is.null(r) || is.null(r[[k]])) NA_real_ else as.numeric(r[[k]])[1], 0)
}

# The p-curve test for left skew (Simonsohn, Nelson and Simmons 2014; Simonsohn's report, point 5):
# every p < .05 times 20, through the normal quantile function, combined by Stouffer's method,
# Z = sum(qnorm(20 p)) / sqrt(k); one-sided p-value 1 - Phi(Z) for left skew. Like p-curve, it
# treats the p-values as independent. qnorm(20 p) is computed on the log scale, so huge |t| do
# not underflow to -Inf.
stouffer_p <- function(t) {
  a   <- abs(t)
  sig <- is.finite(a) & 2 * pnorm(-a) < 0.05
  k   <- sum(sig)
  if (k == 0) return(c(stouffer = NA_real_, stouffer_k = 0))
  z   <- qnorm(log(40) + pnorm(-a[sig], log.p = TRUE), log.p = TRUE)
  c(stouffer = pnorm(sum(z) / sqrt(k), lower.tail = FALSE), stouffer_k = k)
}

# Every test on the analyst's t-statistics under one option; the noisy options average the
# projection test's p-value over CFG$noise_draws de-roundings (the first one is u1, u2)
run_tests <- function(option, tabs, u1, u2, extra) {
  t1  <- analyst_t(option, tabs, u1, u2)
  out <- proj_test(t1)
  if (CFG$noise_draws > 1 && option %in% c("uniform", "smin_noise")) {
    pv <- out[["pval"]]
    bd <- out[["breakdown"]]
    for (j in seq_len(CFG$noise_draws - 1)) {
      pj <- proj_test(analyst_t(option, tabs, extra$u1[, j], extra$u2[, j]))
      pv <- c(pv, pj[["pval"]])
      bd <- c(bd, pj[["breakdown"]])
    }
    out[["pval"]]      <- mean(pv, na.rm = TRUE)   # the application averages p-values and B-hat
    out[["breakdown"]] <- mean(bd, na.rm = TRUE)   # over de-roundings
  }
  if (CFG$run_kudrin)   out <- c(out, kudrin_p(t1))
  if (CFG$run_ekw)      out <- c(out, ekw_p(t1))
  if (CFG$run_stouffer) out <- c(out, stouffer_p(t1))
  out
}

# Repetition r: one world, p-hacked at every tau by every method (tau = 0 once, as "honest"),
# each printed and then de-rounded by every option, all on common random numbers
one_rep <- function(r) {
  set.seed(CFG$seed + 1000 * r)
  w  <- draw_world(RHO, CFG$max_draws)
  mx <- w$T1
  for (k in 2:CFG$max_draws) mx <- pmax(mx, w$H + w$Z[, k])
  hacked <- list(max = mx, threshold = ifelse(w$T1 > CFG$cv, w$T1, pmax(w$T1, w$H + w$Z[, 2])))
  n      <- nrow(DESIGN)
  extra  <- if (CFG$noise_draws > 1) {
    list(u1 = matrix(runif(n * (CFG$noise_draws - 1), -0.5, 0.5), n),
         u2 = matrix(runif(n * (CFG$noise_draws - 1), -0.5, 0.5), n))
  } else NULL
  rows <- list()
  for (tau in CFG$taus) {
    for (m in if (tau == 0) "honest" else CFG$hack_methods) {
      tt   <- if (m == "honest") w$T1 else ifelse(w$hu < tau, hacked[[m]], w$T1)
      tabs <- print_tables(tt, w$x)
      for (o in CFG$options) {
        rows[[length(rows) + 1]] <- data.frame(rep = r, method = m, tau = tau, option = o,
                                               as.list(run_tests(o, tabs, w$u1, w$u2, extra)),
                                               stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, rows)
}

rep_file <- function(r) file.path(OUT_DIR, "reps", sprintf("rep_%05d.rds", as.integer(r)))

one_rep_cached <- function(r) {
  f <- rep_file(r)
  if (CFG$resume && file.exists(f)) {
    z <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(z) && identical(z$key, RUN_KEY)) return(z$res)
  }
  res <- one_rep(r)
  saveRDS(list(key = RUN_KEY, res = res), f)
  res
}


## ---- 6. Summaries, figures and tables ------------------------------------------------------

TESTS     <- c("proj", if (CFG$run_kudrin) c("T_inf", "T_H", "T_H_hybrid"),
               if (CFG$run_ekw) names(EKW_TESTS), if (CFG$run_stouffer) "stouffer")
FIG_TESTS <- intersect(c("proj", "stouffer", "EKW_min", "T_inf", "T_H", "T_H_hybrid"), TESTS)   # figures

# Rejection rates by method, option and tau; the tau = 0 ("honest") results count for every method
summarise_rates <- function(raw) {
  h   <- raw[raw$method == "honest", , drop = FALSE]
  raw <- rbind(raw[raw$method != "honest", , drop = FALSE],
               do.call(rbind, lapply(CFG$hack_methods, function(m) { z <- h; z$method <- m; z })))
  keys <- unique(raw[c("method", "tau", "option")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    d <- raw[raw$method == keys$method[i] & raw$tau == keys$tau[i] & raw$option == keys$option[i], ,
             drop = FALSE]
    # proj_pval: p <= alpha; proj_bhat: B-hat > ALLOWANCE (with warp, B-hat from the pooled cv)
    if (CFG$proj_mode == "warp") {
      ok        <- is.finite(d$resid) & is.finite(d$boot1)
      crit      <- quantile(d$boot1[ok], 1 - CFG$alpha, names = FALSE)
      proj_pval <- mean(d$resid[ok] > crit + d$epsU[ok])
      proj_bhat <- mean(d$resid[ok] > crit + d$epsU[ok] + sqrt(nrow(DESIGN)) * ALLOWANCE)
    } else {
      proj_pval <- mean(d$pval <= CFG$alpha, na.rm = TRUE)
      proj_bhat <- mean(d$breakdown[is.finite(d$breakdown)] > ALLOWANCE)
    }
    smin <- keys$option[i] %in% c("smin", "smin_noise")      # the paper's rule under s_min
    s <- data.frame(method = keys$method[i], tau = keys$tau[i], option = keys$option[i],
                    reps = nrow(d), proj = if (smin) proj_bhat else proj_pval,
                    proj_pval = proj_pval, proj_bhat = proj_bhat,
                    proj_failed = mean(!is.finite(d$resid)), stringsAsFactors = FALSE)
    for (tn in setdiff(TESTS, "proj")) s[[tn]] <- mean(d[[tn]] <= CFG$alpha, na.rm = TRUE)
    if (CFG$run_stouffer) s$stouffer_k <- mean(d$stouffer_k, na.rm = TRUE)   # p-values below .05
    s
  }))
}

OPT_STYLE <- list(
  none       = list(col = "#8C8C8C", lty = "longdash", lwd = 1.3, pch = 4,  label = "No de-rounding"),
  uniform    = list(col = "#E69F00", lty = "dashed",   lwd = 1.3, pch = 2,  label = "Uniform de-rounding"),
  smin       = list(col = "#0072B2", lty = "solid",    lwd = 1.3, pch = 1,  label = "Smallest consistent SE"),
  smin_noise = list(col = "#000000", lty = "solid",    lwd = 1.8, pch = 16, label = "Smallest consistent SE + noise"))
OPT_TEX      <- c(none = "None", uniform = "Unif.", smin = "$s_{\\min}$", smin_noise = "$+$noise")
METHOD_TITLE <- c(max = sprintf("Maximum of %d specifications", as.integer(CFG$max_draws)),
                  threshold = sprintf("Threshold: second draw if t <= %s", format(CFG$cv)))
METHOD_TEX   <- c(max = METHOD_TITLE[["max"]],
                  threshold = sprintf("Threshold: second draw if $t \\le %s$", format(CFG$cv)))
TEST_TITLE   <- c(proj = "Projection test", T_inf = "Kudrin's T_inf", T_H = "Kudrin's T_H",
                  T_H_hybrid = "Kudrin's T_H hybrid", EKW_min = "Any EKW test (smallest p-value)",
                  EKW_LCM = "EKW: LCM", EKW_Fisher = "EKW: Fisher", EKW_Discontinuity = "EKW: Discontinuity",
                  EKW_Binomial = "EKW: Binomial", EKW_CS1 = "EKW: CS1", EKW_CS2B = "EKW: CS2B",
                  stouffer = "p-curve test for left skew (Stouffer)")
TEST_TEX     <- c(proj = "Projection test", T_inf = "$T_\\infty$", T_H = "$T_H$",
                  T_H_hybrid = "$T_H^{\\mathrm{hyb}}$", EKW_min = "Any EKW test (smallest $p$-value)",
                  EKW_LCM = "LCM", EKW_Fisher = "Fisher", EKW_Discontinuity = "Discontinuity",
                  EKW_Binomial = "Binomial", EKW_CS1 = "CS1", EKW_CS2B = "CS2B",
                  stouffer = "$p$-curve test for left skew (Stouffer)")

plot_power <- function(sm, test, stem, note) {
  ms <- CFG$hack_methods
  nl <- length(note)
  for (ext in c("pdf", "png")) {
    if (ext == "pdf") {
      grDevices::pdf(paste0(stem, ".pdf"), width = 3.3 * length(ms), height = 3.58 + 0.12 * nl, pointsize = 8)
    } else {
      grDevices::png(paste0(stem, ".png"), width = 3.3 * length(ms), height = 3.58 + 0.12 * nl, units = "in",
                     res = 300, pointsize = 8)
    }
    par(mfrow = c(1, length(ms)), mar = c(3.1, 3.3, 2.2, 0.8), oma = c(3.2, 0, 1.6 + nl, 0),
        mgp = c(1.9, 0.4, 0), tcl = -0.25, las = 1)
    for (m in ms) {
      d <- sm[sm$method == m, , drop = FALSE]
      plot.new()
      plot.window(xlim = c(-0.03, 1.03), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      half <- qnorm(0.975) * sqrt(CFG$alpha * (1 - CFG$alpha) / max(1, min(d$reps)))
      rect(-0.03, CFG$alpha - half, 1.03, CFG$alpha + half, col = "grey92", border = NA)
      abline(h = c(0.25, 0.5, 0.75), v = c(0.25, 0.5, 0.75), col = "grey88", lwd = 0.5)
      abline(h = CFG$alpha, lty = 3)
      for (o in CFG$options) {
        z  <- d[d$option == o, , drop = FALSE]
        z  <- z[order(z$tau), , drop = FALSE]
        st <- OPT_STYLE[[o]]
        ok <- is.finite(z[[test]])
        lines(z$tau[ok], z[[test]][ok], col = st$col, lty = st$lty, lwd = st$lwd)
        points(z$tau[ok], z[[test]][ok], col = st$col, pch = st$pch, cex = 0.75, xpd = TRUE)
      }
      axis(1, at = seq(0, 1, 0.25), labels = c("0", ".25", ".5", ".75", "1"))
      axis(2, at = seq(0, 1, 0.25), labels = c("0", ".25", ".5", ".75", "1"))
      box()
      mtext(expression("Share of p-hacked results ("*tau*"); at "*tau == 0*": size"), side = 1,
            line = 1.9, cex = 0.85)
      mtext("Rejection rate", side = 2, line = 2.3, las = 0, cex = 0.85)
      mtext(METHOD_TITLE[[m]], side = 3, line = 0.5, font = 2, cex = 0.95)
    }
    mtext(TEST_TITLE[[test]], side = 3, outer = TRUE, line = 0.3 + nl, font = 2, cex = 1)
    for (j in seq_len(nl)) mtext(note[j], side = 3, outer = TRUE, line = 0.2 + (nl - j), cex = 0.7)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    sty <- OPT_STYLE[CFG$options]
    legend("bottom", legend = vapply(sty, function(s) s$label, ""), ncol = 2, bty = "n",
           col = vapply(sty, function(s) s$col, ""), lty = vapply(sty, function(s) s$lty, ""),
           lwd = vapply(sty, function(s) s$lwd, 0), pch = vapply(sty, function(s) s$pch, 0),
           cex = 0.95, seg.len = 2.6, inset = c(0, 0.01), xpd = NA)
    grDevices::dev.off()
  }
  invisible(paste0(stem, c(".pdf", ".png")))
}

# One LaTeX table with notes: rows tau (a block of rows per test when there are several),
# columns p-hacking method x de-rounding option
write_rates_tex <- function(sm, tests, file, caption, notes) {
  ms   <- CFG$hack_methods
  os   <- CFG$options
  no   <- length(os)
  nc   <- 1 + no * length(ms)
  taus <- sort(unique(sm$tau))
  if (!is.null(CFG$table_taus)) taus <- taus[round(taus, 8) %in% round(CFG$table_taus, 8)]
  head <- c(paste0(" & ", paste(sprintf("\\multicolumn{%d}{c}{%s}", no, METHOD_TEX[ms]), collapse = " & "),
                   " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + (seq_along(ms) - 1) * no, 1 + seq_along(ms) * no),
                  collapse = " "),
            paste0("$\\tau$ & ", paste(rep(OPT_TEX[os], length(ms)), collapse = " & "), " \\\\"))
  body <- character(0)
  for (j in seq_along(tests)) {
    tn <- tests[j]
    if (length(tests) > 1) {
      if (j > 1) body <- c(body, "\\midrule")
      body <- c(body, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", nc, TEST_TEX[[tn]]))
    }
    body <- c(body, vapply(taus, function(tau) {
      v <- unlist(lapply(ms, function(m) vapply(os, function(o) {
        z <- sm[sm$method == m & sm$option == o & abs(sm$tau - tau) < 1e-9, tn]
        if (length(z) == 0 || !is.finite(z[1])) "--" else sprintf("%.3f", z[1])
      }, "")))
      paste0(format(tau), " & ", paste(v, collapse = " & "), " \\\\")
    }, ""))
  }
  writeLines(c("\\begin{table}[htbp]", "\\centering", sprintf("\\caption{%s}", caption),
               sprintf("\\label{tab:%s}", tools::file_path_sans_ext(basename(file))),
               "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
               sprintf("\\begin{tabular}{l%s}", strrep("c", no * length(ms))), "\\toprule", head,
               "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\par\\smallskip",
               sprintf("\\begin{minipage}{\\linewidth}\\footnotesize\\textit{Notes:} %s\\end{minipage}", notes),
               "\\end{table}"), file)
}

# The notes every table shares (so each is self-contained), followed by its own
table_notes <- function(extra) {
  big <- function(x) gsub(",", "{,}", format(x, big.mark = ","), fixed = TRUE)
  pit <- if (identical(PI_HAT$type, "gamma")) {
    sprintf("a Gamma with shape %.3f and rate %.3f%s", PI_HAT$shape, PI_HAT$rate,
            if (CFG$pi_fit == "gamma") " fitted to their $|t|$" else ", Kudrin's")
  } else "the NPMLE of the effect distribution fitted to their $|t|$"
  paste0(sprintf("Rejection rates of %s\\%%-level tests over %d repetitions. Each meta-sample has the %s $t$-statistics in %s articles of Brodeur et al.\\ (2020), %s, that the application keeps ($\\kappa \\le %s$), their printed precision and their within-article dependence ($\\rho = %.2f$, in the %s); effects are drawn from %s. ",
                 format(100 * CFG$alpha), as.integer(CFG$reps), big(nrow(DESIGN)), big(N_ART),
                 SUB_LABEL, format(CFG$kappa_cap), RHO, CFG$corr_in, pit),
         "$\\tau$: share of $p$-hacked results ($\\tau = 0$: size). Columns: the $p$-hacking method, then the de-rounding: None ($m/s$), Unif.\\ (uniform draws for the coefficient's last digit and the standard error), $s_{\\min}$ ($m/s_{\\min}$, the smallest consistent standard error), $+$noise ($s_{\\min}$ with uniform noise at the coefficient's last digit). ",
         extra)
}
NOTE_PROJ  <- "Projection test: run\\_test\\_weight\\_fast() with $J = 40$, $J_0 = 20$ and $p = 2$, clustered by article."
NOTE_EKW   <- "EKW tests as in the application, on the two-sided $p$-values below 0.15: LCM, the least concave majorant test, and Fisher, Fisher's test, both of a concave $p$-curve CDF; Discontinuity, the density discontinuity test at 0.05 (Cattaneo et al.\\ 2020); Binomial, the exact binomial test of a non-increasing $p$-curve on $[0.04, 0.05]$; CS1 and CS2B, the Cox and Shi (2023) tests of a non-increasing histogram and of a 2-monotone histogram with bounds (15 bins, clustered by article). Last block: any of the six rejects at 5\\% (smallest $p$-value, no multiple-testing correction)."
NOTE_STOUF <- "$p$-curve test for left skew (Simonsohn, Nelson and Simmons 2014), as suggested in Simonsohn's report: each $p < .05$ is multiplied by 20, transformed by the normal quantile function and combined by Stouffer's method; it rejects for left skew and, like $p$-curve, treats the $p$-values as independent."
NOTE_KUD   <- "Kudrin (2025) tests $T_\\infty$, $T_H$ and $T_H^{\\mathrm{hyb}}$, clustered by article."


## ---- 6b. Exact distortions d_w, to read B-hat in the application ---------------------------
## For each simulated form of p-hacking and each share tau, the d_w of the t-curve it produces in
## the calibrated DGP (effects from Pi-hat, no rounding), computed exactly by integration: the
## weighted distance from its projection-test coefficients (symmetrized, shifted, J0 and p as in
## CFG) to the test's own honest set, the hull of the point-mass vectors b_x on its grid of numgrid
## points in [-L, L] plus the far point 9999 (Lemma 3, eq. (34)). B-hat is a lower 95% confidence
## bound on the same distance for the reported t-curve, so a design whose B-hat exceeds d_w at tau
## is at least as far from honest reporting as that p-hacking by a share tau of studies makes it.

# b_x,j = 2^(-j/2) chi_j(x) phi_2(x) = 2^(-j/2) psi_j(x / sqrt(2)) phi(x / sqrt(2)) / sqrt(2), one
# column per x
point_mass_coeffs <- function(x) {
  J <- CFG$num_coeffs - 1
  t(hf(x / sqrt(2), J) / sqrt(2)) * 2^(-(0:J) / 2)
}

# Weighted distance from theta to the convex hull of the columns of B (nonnegative weights, with
# their sum held at one by a heavy extra row)
hull_distance <- function(theta, B, w) {
  fit <- nnls::nnls(rbind(w * B, 1e3), c(w * theta, 1e3))
  sqrt(sum((w * (B %*% fit$x - theta))^2))
}

# d_w at each tau for "max" and "threshold"; at tau = 0 the t-curve is honest, so that row shows
# the numerical error (zero to many decimals with normal noise)
dw_distortions <- function(taus) {
  J    <- CFG$num_coeffs - 1
  step <- 0.004
  tg   <- seq(-15, 30, by = step)
  XI   <- 0.5 * (hf(tg - CFG$shift, J) + hf(-tg - CFG$shift, J))   # features of a reported t
  if (identical(PI_HAT$type, "gamma")) {
    g  <- effect_grid()
    hs <- g$mids
    ws <- gamma_bins(PI_HAT$shape, PI_HAT$rate, g$edges)
  } else {
    hs <- PI_HAT$h
    ws <- PI_HAT$w
  }
  dn <- function(z) if (is.finite(CFG$nu)) dt(z, CFG$nu) else dnorm(z)
  pn <- function(z) if (is.finite(CFG$nu)) pt(z, CFG$nu) else pnorm(z)
  K  <- CFG$max_draws
  fh <- fm <- ft <- numeric(length(tg))                  # t-curves: honest, max of K, threshold
  for (k in seq_along(hs)) {
    if (ws[k] < 1e-14) next
    d  <- dn(tg - hs[k])
    P  <- pn(tg - hs[k])
    fh <- fh + ws[k] * d
    fm <- fm + ws[k] * K * d * P^(K - 1)
    ft <- ft + ws[k] * ifelse(tg <= CFG$cv, 2 * P * d, d * (1 + pn(CFG$cv - hs[k])))
  }
  th <- lapply(list(honest = fh, max = fm, threshold = ft), function(f) drop(crossprod(XI, f)) * step)
  B  <- cbind(point_mass_coeffs(seq(-CFG$L, CFG$L, length.out = CFG$numgrid)), 0)   # 0: the point 9999
  data.frame(tau = taus,
             max       = vapply(taus, function(a) hull_distance((1 - a) * th$honest + a * th$max, B, DEFF_W), 0),
             threshold = vapply(taus, function(a) hull_distance((1 - a) * th$honest + a * th$threshold, B, DEFF_W), 0))
}

# The share tau at which d (increasing in tau) first reaches b, by linear interpolation; NA when
# even tau = 1 falls short
equivalent_tau <- function(d, taus, b) {
  i <- which(d >= b)[1]
  if (is.na(i)) return(NA_real_)
  if (i == 1) return(taus[1])
  taus[i - 1] + (b - d[i - 1]) / (d[i] - d[i - 1]) * (taus[i] - taus[i - 1])
}

write_dw_tex <- function(DW, EQ, file) {
  num  <- function(v) ifelse(is.finite(v), sprintf("%.3f", v), "--")
  body <- paste0(format(DW$tau), " & ", num(1e3 * DW$max), " & ", num(1e3 * DW$threshold), " \\\\")
  if (!is.null(EQ)) {
    shr  <- function(v) ifelse(is.finite(v), sprintf("%.2f", v), "$>1$")
    body <- c(body, "\\midrule",
              "\\multicolumn{3}{l}{\\textit{Share $\\tau$ at which $d_w$ equals the application's $\\hat B$}} \\\\",
              paste0(EQ$design, " ($\\hat B = ", sprintf("%.1f", 1e3 * EQ$bhat), "$) & ", shr(EQ$tau_max),
                     " & ", shr(EQ$tau_threshold), " \\\\"))
  }
  pit <- if (identical(PI_HAT$type, "gamma")) {
    sprintf("a Gamma with shape %.3f and rate %.3f", PI_HAT$shape, PI_HAT$rate)
  } else "the NPMLE of the effect distribution"
  notes <- paste0(
    "Entries: $d_w \\times 10^3$, the weighted distance between the $t$-curve produced when a share $\\tau$ of studies $p$-hack and the set of $t$-curves possible without selective reporting, computed exactly rather than simulated, with the projection test's coefficients ($J = ",
    CFG$num_coeffs - 1, "$, $J_0 = ", CFG$J0, "$, symmetrized and shifted by ", format(CFG$shift),
    ") and its grid of ", format(CFG$numgrid, big.mark = "{,}"), " points on $[-", format(CFG$L), ", ", format(CFG$L),
    "]$ plus 9999. Effects are drawn from ", pit, ", as in the simulations, with no rounding. Maximum of ",
    CFG$max_draws, ": the largest of ", CFG$max_draws, " $t$-statistics; Threshold: a second draw if $t \\le ",
    format(CFG$cv), "$. The breakdown statistic $\\hat B$ is a lower 95\\% confidence bound on the same distance for the reported $t$-curve, so a design whose $\\hat B$ exceeds $d_w$ at $\\tau$ is at least as far from honest reporting as that form of $p$-hacking by a share $\\tau$ of its studies would make it. The shares depend on the distribution of effects and the form of $p$-hacking; they convey magnitudes, not estimates of the extent of $p$-hacking.")
  writeLines(c("\\begin{table}[htbp]", "\\centering",
               "\\caption{How far simulated $p$-hacking moves the $t$-curve from the honest set}",
               sprintf("\\label{tab:%s_dw}", CFG$name), "\\footnotesize",
               "\\begin{tabular}{lcc}", "\\toprule",
               sprintf("$\\tau$ & Maximum of %d & Threshold \\\\", as.integer(CFG$max_draws)), "\\midrule",
               body, "\\bottomrule", "\\end{tabular}", "\\par\\smallskip",
               sprintf("\\begin{minipage}{\\linewidth}\\footnotesize\\textit{Notes:} %s\\end{minipage}", notes),
               "\\end{table}"), file)
}


## ---- 7. Run --------------------------------------------------------------------------------

DESIGN <- build_design()
N_ART  <- max(DESIGN$gid)
# Rounding allowance (Appendix D.1, step 3): C mean(kappa^6) over the t-statistics used, with
# C = 0.13 for J0 = 20 and 0.31 for J0 = 30 (Section 5.4, weights min(1, 2^(J0 - j)))
ALLOW_C <- if (!is.null(CFG$allowance_C)) CFG$allowance_C else
  if (CFG$p_weight == 2 && CFG$J0 == 20) 0.13 else if (CFG$p_weight == 2 && CFG$J0 == 30) 0.31 else NA_real_
if (!is.finite(ALLOW_C) && any(CFG$options %in% c("smin", "smin_noise"))) {
  stop("The paper gives the rounding allowance only for p = 2 and J0 = 20 (0.13) or 30 (0.31); ",
       "set CFG$allowance_C for J0 = ", CFG$J0, " and p = ", CFG$p_weight, ".")
}
ALLOWANCE <- ALLOW_C * mean(DESIGN$kappa^6)

DATA_KEY <- paste(normalizePath(DATA_FILE), file.info(DATA_FILE)$size, format(file.info(DATA_FILE)$mtime))
CAL_KEY  <- paste(DATA_KEY, paste(SUBSAMPLE, collapse = "+"), CFG$kappa_cap, CFG$se_rounded_up, CFG$nu,
                  CFG$fit_option, CFG$deff_option, CFG$corr_in, CFG$n_cal, CFG$rho_max, CFG$se_cell, CFG$npmle_iters,
                  CFG$pi_fit, CFG$fit_cap, CFG$compare_kudrin,
                  CFG$seed, CFG$J0, CFG$p_weight, CFG$num_coeffs, CFG$shift)
cal_file <- file.path(OUT_DIR, paste0(CFG$name, "_calibration.rds"))
CAL <- if (file.exists(cal_file)) readRDS(cal_file) else NULL
if (is.null(CAL) || !identical(CAL$key, CAL_KEY)) {
  set.seed(CFG$seed + 11)
  t_fit  <- data_t(CFG$fit_option, runif(nrow(DESIGN), -0.5, 0.5), runif(nrow(DESIGN), -0.5, 0.5))
  if (CFG$pi_fit %in% c("kudrin_rate", "kudrin_scale")) {
    message("Effects: Kudrin's Gamma(shape 1.4134, ", if (CFG$pi_fit == "kudrin_rate") "rate" else "scale",
            " 1.2980), used as given; nothing is estimated.")
  } else {
    message("Fitting Pi-hat (", CFG$pi_fit, ") to the |t| of the subsample ...")
  }
  PI_HAT <- switch(CFG$pi_fit,
                   gamma        = fit_gamma(abs(t_fit)),
                   npmle        = fit_pi_hat(abs(t_fit)),
                   kudrin_rate  = gamma_fixed(1.4134, 1.2980),
                   kudrin_scale = gamma_fixed(1.4134, 1 / 1.2980))
  if (isTRUE(CFG$compare_kudrin)) {
    message("For the report only: log-likelihood of Kudrin's two readings on these |t| (no fitting) ...")
    PI_HAT$kudrin_loglik <- gamma_logliks(abs(t_fit), list(rate  = c(1.4134, 1.2980),
                                                           scale = c(1.4134, 1 / 1.2980)))
  }
  set.seed(CFG$seed + 12)
  t_deff <- data_t(CFG$deff_option, runif(nrow(DESIGN), -0.5, 0.5), runif(nrow(DESIGN), -0.5, 0.5))
  target <- deff_trace(t_deff, DESIGN$gid)
  message(sprintf("Calibrating rho to the data's design effect %.3f ...", target))
  RHO_FIT <- calibrate_rho(target)
  CAL <- list(key = CAL_KEY, pi_hat = PI_HAT, rho = RHO_FIT, data_sig = mean(abs(t_fit) > 1.96),
              data_q = quantile(abs(t_fit), c(0.1, 0.25, 0.5, 0.75, 0.9)))
  saveRDS(CAL, cal_file)
}
PI_HAT <- CAL$pi_hat
RHO    <- CAL$rho$rho

# Checks of the construction (cheap)
local({
  chk <- function(ok, what) if (!isTRUE(all(ok))) stop("Check failed: ", what)
  r   <- DESIGN$ratio
  chk(abs(DESIGN$S[r] * DESIGN$u_se[r] - DESIGN$sd_p[r]) <= 1e-9 * DESIGN$sd_p[r],
      "the printed se is S units of u")
  chk(DESIGN$c_rel[r] >= 1 & DESIGN$a_rel[r] >= 1 & DESIGN$b_rel[r] >= 1, "units no finer than u")
  chk(DESIGN$kappa[r] <= CFG$kappa_cap + 1e-12, "kappa under the cap")
  set.seed(1)
  u1 <- runif(nrow(DESIGN), -0.5, 0.5)
  u2 <- runif(nrow(DESIGN), -0.5, 0.5)
  td <- list(m = DESIGN$mu_p / DESIGN$u_se, tp = DESIGN$t_p)   # the data's own tables
  for (o in CFG$options) {
    chk(abs(analyst_t(o, td, u1, u2) - data_t(o, u1, u2)) <= 1e-8 * pmax(1, abs(data_t(o, u1, u2))),
        paste("the analyst reproduces the application's", o))
  }
  w <- draw_world(RHO, 2)
  chk(w$x[r] >= DESIGN$smin_rel[r] - 1e-12, "true se >= s_min (Proposition 4)")
  chk(round(w$x[r]) == DESIGN$S[r] | CFG$se_rounded_up, "the true se prints as the data's se")
  message("Construction checks passed.")
})

BY_METHOD <- table(DESIGN$method)
report <- c(
  sprintf("Subsample (%s): %s t-statistics in %s articles (%.1f per article, at most %d); %.1f%% coefficient/SE, %.1f%% printed t.",
          SUB_LABEL, format(nrow(DESIGN), big.mark = ","), format(N_ART, big.mark = ","), nrow(DESIGN) / N_ART,
          max(tabulate(DESIGN$gid)), 100 * mean(DESIGN$ratio), 100 * mean(!DESIGN$ratio)),
  if (length(BY_METHOD) > 1) paste0("By method: ", paste(sprintf("%s %s", names(BY_METHOD),
                                                                 format(as.vector(BY_METHOD), big.mark = ",")), collapse = ", "), "."),
  sprintf("Rounding: kappa median %.3f, 90th percentile %.3f; se significands (units of u) below 10: %.1f%%; rounding allowance %s x mean(kappa^6) = %.2e, which B-hat must exceed under s_min.",
          median(DESIGN$kappa), quantile(DESIGN$kappa, 0.9), 100 * mean(DESIGN$S[DESIGN$ratio] < 10),
          format(ALLOW_C), ALLOWANCE),
  if (identical(PI_HAT$type, "gamma")) {
    paste0(sprintf("Pi-hat: Gamma(shape %.4f, rate %.4f) = scale %.4f; mean %.2f, sd %.2f; ",
                   PI_HAT$shape, PI_HAT$rate, 1 / PI_HAT$rate, PI_HAT$shape / PI_HAT$rate,
                   sqrt(PI_HAT$shape) / PI_HAT$rate),
           if (CFG$pi_fit == "gamma") {
             sprintf("fitted by maximum likelihood (log-likelihood %.1f; %d |t| above %s censored); ",
                     PI_HAT$loglik, as.integer(PI_HAT$n_censored), format(CFG$fit_cap))
           } else {
             sprintf("Kudrin's values with 1.2980 a %s, not estimated; ", if (CFG$pi_fit == "kudrin_rate") "rate" else "scale")
           },
           if (is.null(PI_HAT$kudrin_loglik)) "" else {
             sprintf("log-likelihood of Kudrin's 1.4134, 1.2980 on these |t|: %.1f as a rate, %.1f as a scale; ",
                     PI_HAT$kudrin_loglik[["rate"]], PI_HAT$kudrin_loglik[["scale"]])
           },
           sprintf("P(|t| > 1.96): data %.3f, implied %.3f.", CAL$data_sig, PI_HAT$sig_share))
  } else {
    sprintf("Pi-hat (npmle; %d |t| above %s censored): %d atoms, mean %.2f, sd %.2f (%d EM iterations); log-likelihood %.1f; P(|t| > 1.96): data %.3f, fitted %.3f.",
            as.integer(PI_HAT$n_censored), format(CFG$fit_cap), length(PI_HAT$h), sum(PI_HAT$w * PI_HAT$h), sqrt(sum(PI_HAT$w * PI_HAT$h^2) - sum(PI_HAT$w * PI_HAT$h)^2),
            PI_HAT$iters, PI_HAT$loglik, CAL$data_sig, PI_HAT$sig_share)
  },
  sprintf("Within-article dependence (in the %s): design effect of the data %.3f; simulated %.3f at rho = %.3f (rho = 0: %.3f, rho = %s: %.3f).",
          CFG$corr_in, CAL$rho$target, CAL$rho$deff, RHO, CAL$rho$deff_0, format(CFG$rho_max), CAL$rho$deff_max))
message(paste(report, collapse = "\n"))
writeLines(report, file.path(OUT_DIR, paste0(CFG$name, "_calibration.txt")))

# Exact distortions d_w of the simulated p-hacking, and the share tau each B-hat corresponds to
if (requireNamespace("nnls", quietly = TRUE)) {
  DW <- dw_distortions(sort(unique(c(0, CFG$dw_taus))))
  EQ <- if (length(CFG$bhat) > 0) {
    data.frame(design = names(CFG$bhat), bhat = as.numeric(CFG$bhat),
               tau_max       = vapply(as.numeric(CFG$bhat), function(b) equivalent_tau(DW$max, DW$tau, b), 0),
               tau_threshold = vapply(as.numeric(CFG$bhat), function(b) equivalent_tau(DW$threshold, DW$tau, b), 0),
               stringsAsFactors = FALSE)
  } else NULL
  message("Exact distortion d_w (x 1000) of the simulated p-hacking (effects from Pi-hat, no rounding):")
  print(data.frame(tau = DW$tau, max_of_K = round(1e3 * DW$max, 3), threshold = round(1e3 * DW$threshold, 3)),
        row.names = FALSE)
  if (!is.null(EQ)) {
    message("Share tau at which d_w equals each B-hat (NA: beyond tau = 1):")
    print(data.frame(design = EQ$design, bhat_x1000 = 1e3 * EQ$bhat, tau_max = round(EQ$tau_max, 2),
                     tau_threshold = round(EQ$tau_threshold, 2)), row.names = FALSE)
  }
  write.csv(DW, file.path(OUT_DIR, paste0(CFG$name, "_dw.csv")), row.names = FALSE)
  if (!is.null(EQ)) write.csv(EQ, file.path(OUT_DIR, paste0(CFG$name, "_dw_bhat.csv")), row.names = FALSE)
  write_dw_tex(DW, EQ, file.path(OUT_DIR, paste0(CFG$name, "_dw.tex")))
} else {
  message("Skipping the d_w table: it needs the nnls package (install.packages(\"nnls\")).")
}

proj_key <- paste(CFG$num_coeffs, CFG$L, CFG$numgrid, if (is.null(CFG$epsilon_U_pre)) "computed" else CFG$epsilon_U_pre)
proj_file <- file.path(OUT_DIR, "proj_constants.rds")
PROJ <- if (file.exists(proj_file)) readRDS(proj_file) else NULL
if (is.null(PROJ) || !identical(PROJ$key, proj_key)) {
  eps <- CFG$epsilon_U_pre
  if (is.null(eps)) {
    message("Computing epsilon_U once; this is slow ...")
    eps <- compute_epsilons(CFG$L, nx = CFG$numgrid)$epsilon_U
  }
  PROJ <- list(key = proj_key, epsilon_U_pre = eps,
               U = create_basis_stable(CFG$num_coeffs, L = CFG$L, numgrid = CFG$numgrid, sigma_Y = 1))
  saveRDS(PROJ, proj_file)
}

PI_KEY  <- if (identical(PI_HAT$type, "gamma")) paste(PI_HAT$shape, PI_HAT$rate) else sum(PI_HAT$h * PI_HAT$w)
RUN_KEY <- paste(CAL_KEY, RHO, PI_KEY, paste(CFG$taus, collapse = ","),
                 paste(CFG$hack_methods, collapse = ","), CFG$max_draws, CFG$cv, CFG$hack_level,
                 paste(CFG$options, collapse = ","), CFG$noise_draws, CFG$proj_mode, CFG$proj_boots,
                 CFG$L, CFG$numgrid, PROJ$epsilon_U_pre, CFG$run_kudrin, CFG$kudrin_boots, CFG$run_ekw,
                 CFG$run_stouffer, "tests v2",   # v2: every EKW test separately, Stouffer
                 if (CFG$noise_draws > 1) "B-hat averaged" else "")

# Run time: one projection test now, times the tests per repetition
t0 <- Sys.time()
invisible(proj_test(data_t("smin")))
secs  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
per   <- length(CFG$options) * (1 + length(CFG$hack_methods) * sum(CFG$taus > 0))
per   <- per + 2 * (CFG$noise_draws - 1) * (1 + length(CFG$hack_methods) * sum(CFG$taus > 0)) *
  sum(CFG$options %in% c("uniform", "smin_noise")) / 2
message(sprintf("One projection test takes %.1f s; %d per repetition; %d repetitions on %d cores: about %.1f hours (more with Kudrin or EKW).",
                secs, as.integer(per), as.integer(CFG$reps), as.integer(CFG$n_cores),
                secs * per * CFG$reps / CFG$n_cores / 3600))

cl <- NULL
if (CFG$n_cores > 1) {
  cl <- parallel::makeCluster(CFG$n_cores)
  parallel::clusterExport(cl, "CODE_FILES", envir = globalenv())
  parallel::clusterEvalQ(cl, { for (f in CODE_FILES) source(f); NULL })
  parallel::clusterExport(cl, c("CFG", "OUT_DIR", "DESIGN", "N_ART", "PI_HAT", "RHO", "PROJ", "RUN_KEY",
                                "run_ewk", "draw_noise", "draw_effects", "draw_world", "print_tables",
                                "analyst_t", "proj_test", "kudrin_p", "ekw_p", "run_tests", "one_rep",
                                "rep_file", "one_rep_cached", "EKW_TESTS", "stouffer_p"), envir = globalenv())
}
RAW <- tryCatch({
  acc   <- NULL
  batch <- 4 * max(1L, CFG$n_cores)
  for (b0 in seq(1, CFG$reps, by = batch)) {
    rr  <- b0:min(CFG$reps, b0 + batch - 1)
    t0  <- Sys.time()
    out <- if (is.null(cl)) lapply(rr, one_rep_cached) else parallel::parLapply(cl, rr, one_rep_cached)
    acc <- rbind(acc, do.call(rbind, out))
    message(sprintf("Repetitions %d-%d of %d done (%.1f min).", min(rr), max(rr), as.integer(CFG$reps),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  acc
}, finally = {
  if (!is.null(cl)) parallel::stopCluster(cl)
})

SM <- summarise_rates(RAW)
SM <- SM[order(match(SM$method, CFG$hack_methods), match(SM$option, CFG$options), SM$tau), ]
rownames(SM) <- NULL
write.csv(SM, file.path(OUT_DIR, paste0(CFG$name, "_rates.csv")), row.names = FALSE)
print(SM[, c("method", "option", "tau", "reps", FIG_TESTS)], digits = 3, row.names = FALSE)

# Two lines above every figure: the sample, then the effect distribution and the dependence
PI_FIG <- if (identical(PI_HAT$type, "gamma")) {
  sprintf("Gamma(shape %.3f, rate %.3f), i.e. scale %.3f and mean %.2f, %s",
          PI_HAT$shape, PI_HAT$rate, 1 / PI_HAT$rate, PI_HAT$shape / PI_HAT$rate,
          switch(CFG$pi_fit, gamma = "fitted by maximum likelihood to |t|",
                 kudrin_rate = "Kudrin's, with 1.2980 a rate", kudrin_scale = "Kudrin's, with 1.2980 a scale"))
} else {
  sprintf("NPMLE fitted to |t| (%d atoms, mean %.2f, sd %.2f)", length(PI_HAT$h), sum(PI_HAT$w * PI_HAT$h),
          sqrt(sum(PI_HAT$w * PI_HAT$h^2) - sum(PI_HAT$w * PI_HAT$h)^2))
}
note <- c(sprintf("Brodeur et al. (2020), %s, kappa <= %s: %s t-statistics in %s articles; %d repetitions; shaded: 95%% band around %s%%.",
                  SUB_LABEL, format(CFG$kappa_cap), format(nrow(DESIGN), big.mark = ","), format(N_ART, big.mark = ","),
                  as.integer(CFG$reps), format(100 * CFG$alpha)),
          sprintf("Effects ~ %s; within-article rho = %.2f (in the %s).", PI_FIG, RHO, CFG$corr_in))
RULE_FIG <- sprintf("Under s_min and s_min + noise the projection test rejects only if B-hat > rounding allowance %s x mean(kappa^6) = %.1e; otherwise if p <= %s.",
                    format(ALLOW_C), ALLOWANCE, format(CFG$alpha))
for (tn in FIG_TESTS) {
  plot_power(SM, tn, file.path(OUT_DIR, paste0(CFG$name, "_", tn)), if (tn == "proj") c(note, RULE_FIG) else note)
}
like_bch <- function(what) sprintf("%s on meta-samples like Brodeur et al.\\ (2020), %s", what, SUB_LABEL)
RULE_TEX <- sprintf(" Under $s_{\\min}$ and $+$noise the test rejects only when its breakdown statistic $\\hat B$ exceeds the rounding allowance %s times the mean of $\\kappa^6$, %s (Appendix D.1, step 3), so that a rejection cannot be attributed to rounding; under None and Unif.\\ it rejects when $p \\le %s$.",
                    format(ALLOW_C), formatC(ALLOWANCE, format = "e", digits = 2), format(CFG$alpha))
TABLES <- list(proj = list(tests = "proj", notes = paste0(NOTE_PROJ, RULE_TEX), caption = like_bch("Projection test")))
if (CFG$run_stouffer) {
  TABLES$stouffer <- list(tests = "stouffer", notes = NOTE_STOUF,
                          caption = like_bch("The $p$-curve test for left skew (Stouffer)"))
}
if (CFG$run_ekw) {
  TABLES$EKW <- list(tests = names(EKW_TESTS), notes = NOTE_EKW,
                     caption = like_bch("Each EKW test"))
}
if (CFG$run_kudrin) {
  TABLES$kudrin <- list(tests = c("T_inf", "T_H", "T_H_hybrid"), notes = NOTE_KUD,
                        caption = like_bch("Kudrin's (2025) tests"))
}
for (nm in names(TABLES)) {
  write_rates_tex(SM, TABLES[[nm]]$tests, file.path(OUT_DIR, paste0(CFG$name, "_", nm, ".tex")),
                  TABLES[[nm]]$caption, table_notes(TABLES[[nm]]$notes))
}
message("Done: figures, rates and tables in ", OUT_DIR)