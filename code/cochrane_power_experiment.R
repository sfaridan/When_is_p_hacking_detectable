################################################################################
#  cochrane_power_experiment.R
#
#  Power against outcome switching, built from real trials. A randomized trial in the Cochrane
#  Database of Systematic Reviews usually reports several efficacy outcomes, so a meta-sample in
#  which some trials report their most significant outcome can be built from the trials themselves:
#  the effects, the correlation between a trial's outcomes and the noise are all real, and only the
#  selection, which is known, is imposed. The data are the continuous outcomes that
#  application_cochrane.R keeps (randomized, efficacy, at least n_min participants, repeated copies
#  of the same outcome data removed); run it once with continuous_only <- TRUE (inspect_only <- TRUE
#  is enough): it writes the sample file this script reads.
#
#  The DGP, for each trial in a meta-sample of n_trials trials drawn without replacement from the
#  trials with at least min_outcomes outcomes (a trial's outcomes: its first max_outcomes in review
#  order, so "max" is a best of at most 8, as in simulations_brodeur_rct.R):
#    Unrounded   the means and SDs are printed rounded; each repetition draws them once, uniformly
#                within their rounding cells, which gives every outcome its unrounded z
#    Honest      the trial reports one of its outcomes chosen at random (its "pre-specified" one)
#    Switching   with probability tau the trial instead reports, judging by unrounded |z|:
#                  "max"        its most significant outcome
#                  "threshold"  the pre-specified outcome if |z| > 1.96, else the larger of it
#                               and one other outcome chosen at random
#    Analyst     with a the unit of the means and b that of the SDs (as application_cochrane.R reads
#                them), s_j,min = max(s_j - b/2, 0) the smallest SDs consistent with the printed ones,
#                and SE_min the standard error from them (the smallest consistent one):
#                  "none"        the printed z
#                  "uniform"     means and SDs drawn uniformly within their rounding cells
#                  "smin"        (m1 - m2) / SE_min
#                  "smin_noise"  (m1 + V1 - m2 - V2) / SE_min, V_j ~ U(-a/2, a/2), the same V as "uniform"
#                The mean difference of two means each rounded to a can be off by up to a, so its
#                rounding unit in z-units is kappa = 2a / SE_min (the paper's kappa = a / s_min is the
#                width of the rounding interval in t-units). Outcomes with kappa > kappa_cap are left
#                out under every option, as the application does; this selects on how results are
#                printed, never on estimates. Under "smin" and "smin_noise" the projection test
#                rejects only when B-hat exceeds the rounding allowance C mean(kappa^6) of the
#                reported outcomes (C = 0.13 for J0 = 20, 0.31 for J0 = 30; Propositions 3 and 4),
#                as in the paper's tables; under "none" and "uniform" when p <= alpha
#  One reported z per trial, so trials are independent units (each its own cluster). Common random
#  numbers: within a repetition every tau and method share the trials, the pre-specified and other
#  outcomes, the switching draws and the rounding draws; tau = 0 gives size.
#
#  Tests and outputs as in simulations_brodeur_rct.R, in OUT_DIR:
#    NAME_<test>.pdf/.png     power curves of each test, one panel per method, one line per option
#    NAME_compare_<opt>.pdf/.png  every test in one figure, one panel per method (as the figures of
#                             kudrin_figures_engine.R), for each option
#    NAME_main.tex / .csv     the main results: size and power at chosen taus, every test, by option
#    NAME_<test>.tex          the tables of simulations_brodeur_rct.R; NAME_rates.csv every rate
#    NAME_sample.txt          the trials, their outcomes, and how much switching changes significance
#    reps/rep_NNNNN.rds       one file per repetition, so an interrupted run resumes
#
#  Run once with QUICK <- TRUE to check that everything works.
################################################################################

rm(list = ls())

## ---- 0. Settings -------------------------------------------------------------------------

ROOT     <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
CODE_DIR <- file.path(ROOT, "code")
QUICK    <- TRUE

CFG <- list(
  name          = "cochrane_switching_power",
  sample_file   = file.path(ROOT, "results", "cochrane", "Cochrane_sample_n100_all_trial_kappa0.5_cont.rds"),
  
  # Meta-samples
  n_trials      = 1500,                  # trials per meta-sample (one reported outcome each)
  min_outcomes  = 2,                     # eligible trials: at least this many outcomes
  max_outcomes  = 8,                     # a trial chooses among its first max_outcomes outcomes
  
  # Selection
  taus          = seq(0, 1, by = 0.1),   # share of trials that switch; tau = 0 gives size
  hack_methods  = c("max", "threshold"),
  cv            = 1.96,                  # "threshold": switch only when |z| <= cv
  
  # The analyst's treatment of rounding: the application's four options
  options       = c("none", "uniform", "smin", "smin_noise"),
  kappa_cap     = 0.5,                   # leave out outcomes with kappa = 2a / SE_min above this (all options); the
  # application's sample already applies its own cap, so this only matters if lower
  allowance_C   = NULL,                  # rounding allowance constant; NULL: 0.13 (J0 = 20) or 0.31 (J0 = 30)
  
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
  
  # Other tests (each adds its own run time)
  run_kudrin    = TRUE,
  kudrin_boots  = 199,
  run_ekw       = TRUE,
  run_stouffer  = TRUE,
  table_taus    = NULL,                  # taus in the tables of each test (NULL: all)
  main_taus     = c( 0.5, 1),       # taus in the main results table, besides tau = 0
  
  # Computing
  reps          = 500,
  alpha         = 0.05,
  n_cores       = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  seed          = 20261031,
  resume        = TRUE
)
OUT_DIR <- file.path(ROOT, "simulations", CFG$name)
if (QUICK) {
  CFG$reps          <- 50
  CFG$taus          <- seq(0, 1, by = 0.1) #c(0,0.25, 0.5, 0.75, 1)
  CFG$proj_boots    <- 99
  CFG$kudrin_boots  <- 99
  CFG$epsilon_U_pre <- 0.00063
  OUT_DIR           <- paste0(OUT_DIR, "_quick")
}
dir.create(file.path(OUT_DIR, "reps"), recursive = TRUE, showWarnings = FALSE)
stopifnot(all(CFG$options %in% c("none", "uniform", "smin", "smin_noise")), all(CFG$hack_methods %in% c("max", "threshold")),
          CFG$min_outcomes >= 2, CFG$max_outcomes >= CFG$min_outcomes, CFG$proj_mode %in% c("full", "warp"),
          any(CFG$taus == 0), all(CFG$taus >= 0 & CFG$taus <= 1), all(CFG$main_taus %in% CFG$taus))


## ---- 1. Your code ------------------------------------------------------------------------

CODE_FILES <- file.path(CODE_DIR, c("functions_all_weighted.R",
                                    if (CFG$run_ekw) "Tests.R",
                                    if (CFG$run_kudrin) "kudrin_tests.R"))
missing_files <- CODE_FILES[!file.exists(CODE_FILES)]
if (length(missing_files) > 0) stop("Not found: ", paste(missing_files, collapse = ", "))
for (f in CODE_FILES) source(f)
stopifnot(exists("run_test_weight_fast"), exists("create_basis_stable"), exists("compute_epsilons"))



## ---- 2. The trials, as application_cochrane.R keeps them ----------------------------------

if (!file.exists(CFG$sample_file)) {
  stop("No ", basename(CFG$sample_file), ": run application_cochrane.R with continuous_only <- TRUE ",
       "(inspect_only <- TRUE is enough), which writes it")
}
SAMPLE <- readRDS(CFG$sample_file)
if (!isTRUE(SAMPLE$settings$continuous_only) || !identical(SAMPLE$settings$outcomes_per_trial, "all")) {
  stop("The sample file must come from application_cochrane.R with continuous_only <- TRUE and outcomes_per_trial <- \"all\"")
}
inf_last <- function(x) ifelse(is.finite(x), x, Inf)
D <- SAMPLE$sample
D <- D[D$source == "continuous" & is.finite(D$z) & is.finite(D$z_raw), ]
D <- D[order(D$unit, D$review, inf_last(D$comparison), inf_last(D$outcome), D$row), ]   #review order within a trial
D$se     <- sqrt(D$s1^2 / D$n1 + D$s2^2 / D$n2)                                     #from the printed SDs
D$se_min <- sqrt(pmax(D$s1 - D$b_sd / 2, 0)^2 / D$n1 + pmax(D$s2 - D$b_sd / 2, 0)^2 / D$n2)   #smallest consistent
D$kappa  <- 2 * D$a_mean / D$se_min                                                   #rounding unit of m1 - m2, z-units
N_BEFORE_CAP <- c(outcomes = nrow(D), trials = length(unique(D$unit)))
keep_k <- D$se_min > 0 & is.finite(D$kappa)                                          #s_min needs SE_min > 0
if (!is.null(CFG$kappa_cap)) keep_k <- keep_k & D$kappa <= CFG$kappa_cap
D <- D[keep_k, ]
N_AFTER_CAP <- nrow(D)
D$k <- ave(seq_len(nrow(D)), D$unit, FUN = seq_along)
D   <- D[D$k <= CFG$max_outcomes, ]                                                  #a trial's first max_outcomes outcomes
K_ALL  <- table(D$unit)
D      <- D[D$unit %in% names(K_ALL)[K_ALL >= CFG$min_outcomes], ]
rownames(D) <- NULL
TRIALS <- split(seq_len(nrow(D)), D$unit)                                            #each trial's rows of D
if (CFG$n_trials > length(TRIALS)) stop("n_trials = ", CFG$n_trials, " but only ", length(TRIALS), " trials are eligible")

# z of rows i with their means and SDs drawn within their rounding cells (U: length(i) x 4 uniforms on
# [0, 1)), as application_cochrane.R de-rounds: the change in the z from the means and SDs is added to z
z_at <- function(i, U) {
  a <- D$a_mean[i]; b <- D$b_sd[i]
  sd_draw <- function(s, u) { lo <- pmax(s - b / 2, 0); lo + (s + b / 2 - lo) * u }
  m1 <- D$m1[i] + a * (U[, 1] - 0.5); m2 <- D$m2[i] + a * (U[, 2] - 0.5)
  s1 <- sd_draw(D$s1[i], U[, 3]);     s2 <- sd_draw(D$s2[i], U[, 4])
  zc <- (m1 - m2) / sqrt(s1^2 / D$n1[i] + s2^2 / D$n2[i])
  ifelse(is.finite(zc), D$z[i] + (zc - D$z_raw[i]), D$z[i])
}


# The analyst's z under option o for rows i; U: the analyst's uniforms for those rows, shared by "uniform"
# (means and SDs) and "smin_noise" (means). The file's z is rescaled from the printed SE to SE_min, so it
# keeps its own definition, as the uniform de-rounding adds its change to it
z_option <- function(o, i, U) {
  switch(o,
         none       = D$z[i],
         uniform    = z_at(i, U),
         smin       = D$z[i] * D$se[i] / D$se_min[i],
         smin_noise = D$z[i] * D$se[i] / D$se_min[i] + D$a_mean[i] * (U[, 1] - U[, 2]) / D$se_min[i])
}

ALLOW_C <- if (!is.null(CFG$allowance_C)) CFG$allowance_C else
  if (CFG$p_weight == 2 && CFG$J0 == 20) 0.13 else if (CFG$p_weight == 2 && CFG$J0 == 30) 0.31 else NA_real_
if (!is.finite(ALLOW_C) && any(CFG$options %in% c("smin", "smin_noise"))) {
  stop("The paper gives the rounding allowance only for p = 2 and J0 = 20 (0.13) or 30 (0.31); set CFG$allowance_C")
}


## ---- 3. Tests: as simulations_brodeur_rct.R, clustered by trial ------------------------------

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

proj_test <- function(t, ids) {
  pr <- tryCatch(run_test_weight_fast(data.frame(t = t, title = ids, stringsAsFactors = FALSE),
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

kudrin_p <- function(t, ids) {
  nm <- c("T_inf", "T_H", "T_H_hybrid")
  kt <- tryCatch(kudrin_tests(t, cluster = ids, boots = CFG$kudrin_boots), error = function(e) NULL)
  if (is.null(kt)) return(setNames(rep(NA_real_, 3), nm))
  vapply(nm, function(k) if (k %in% names(kt$pval)) as.numeric(kt$pval[[k]])[1] else NA_real_, 0)
}

ekw_p <- function(t, ids) {
  r <- tryCatch(run_ewk(2 * pnorm(-abs(t)), ids = as.integer(factor(ids))), error = function(e) NULL)
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

run_tests <- function(t, ids) {
  out <- proj_test(t, ids)
  if (CFG$run_kudrin)   out <- c(out, kudrin_p(t, ids))
  if (CFG$run_ekw)      out <- c(out, ekw_p(t, ids))
  if (CFG$run_stouffer) out <- c(out, stouffer_p(t))
  out
}


## ---- 4. One repetition ---------------------------------------------------------------------

# Repetition r: one meta-sample of trials, switched at every tau by every method (tau = 0 once, as
# "honest"), each seen by the analyst under every option, all on common random numbers
one_rep <- function(r) {
  set.seed(CFG$seed + 1000 * r)
  tr    <- sample.int(length(TRIALS), CFG$n_trials)
  rows  <- TRIALS[tr]
  ids   <- names(TRIALS)[tr]
  nt    <- length(rows)
  pick  <- function(ix) ix[sample.int(length(ix), 1L)]
  pre   <- vapply(rows, pick, 0L)                                              #pre-specified outcome
  oth   <- vapply(seq_len(nt), function(j) pick(setdiff(rows[[j]], pre[j])), 0L) #another outcome
  all_i <- unlist(rows, use.names = FALSE)
  z_true <- rep(NA_real_, nrow(D))
  z_true[all_i] <- z_at(all_i, matrix(runif(4 * length(all_i)), ncol = 4))     #unrounded: trials select on it
  U_an  <- matrix(runif(4 * length(all_i)), ncol = 4)                           #the analyst's de-rounding draws
  z_opt <- lapply(setNames(CFG$options, CFG$options), function(o) {
    z <- rep(NA_real_, nrow(D)); z[all_i] <- z_option(o, all_i, U_an); z
  })
  best   <- vapply(rows, function(ix) ix[which.max(abs(z_true[ix]))], 0L)
  thr    <- ifelse(abs(z_true[pre]) > CFG$cv, pre, ifelse(abs(z_true[oth]) > abs(z_true[pre]), oth, pre))
  hacked <- list(max = best, threshold = thr)
  hu     <- runif(nt)
  out    <- list()
  for (tau in CFG$taus) {
    for (m in if (tau == 0) "honest" else CFG$hack_methods) {
      ch <- if (m == "honest") pre else ifelse(hu < tau, hacked[[m]], pre)
      allowance <- ALLOW_C * mean(D$kappa[ch]^6)                             #of the reported outcomes
      for (o in CFG$options) {
        out[[length(out) + 1]] <- data.frame(rep = r, method = m, tau = tau, option = o, as.list(run_tests(z_opt[[o]][ch], ids)),
                                             switched = mean(ch != pre), allowance = allowance, stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, out)
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


## ---- 5. Summaries, figures and tables ------------------------------------------------------

TESTS     <- c("proj", if (CFG$run_kudrin) c("T_inf", "T_H", "T_H_hybrid"),
               if (CFG$run_ekw) names(EKW_TESTS), if (CFG$run_stouffer) "stouffer")
FIG_TESTS <- intersect(c("proj", "stouffer", "EKW_min", "T_inf", "T_H", "T_H_hybrid"), TESTS)   # figures
EKW_EACH  <- setdiff(names(EKW_TESTS), "EKW_min")

# Rejection rates by method, option and tau, as simulations_brodeur_rct.R: under "smin" and "smin_noise"
# the projection test rejects only when B-hat exceeds the reported outcomes' rounding allowance, under the
# others when p <= alpha; the tau = 0 ("honest") results count for every method. EKW_max: the largest
# rejection rate of the six EKW tests, as in Kudrin's figures.
summarise_rates <- function(raw) {
  h   <- raw[raw$method == "honest", , drop = FALSE]
  raw <- rbind(raw[raw$method != "honest", , drop = FALSE],
               do.call(rbind, lapply(CFG$hack_methods, function(m) { z <- h; z$method <- m; z })))
  keys <- unique(raw[c("method", "tau", "option")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    d <- raw[raw$method == keys$method[i] & raw$tau == keys$tau[i] & raw$option == keys$option[i], , drop = FALSE]
    if (CFG$proj_mode == "warp") {
      ok        <- is.finite(d$resid) & is.finite(d$boot1)
      crit      <- quantile(d$boot1[ok], 1 - CFG$alpha, names = FALSE)
      proj_pval <- mean(d$resid[ok] > crit + d$epsU[ok])
      proj_bhat <- mean(d$resid[ok] > crit + d$epsU[ok] + sqrt(CFG$n_trials) * d$allowance[ok])
    } else {
      proj_pval <- mean(d$pval <= CFG$alpha, na.rm = TRUE)
      ok        <- is.finite(d$breakdown)
      proj_bhat <- mean(d$breakdown[ok] > d$allowance[ok])
    }
    smin <- keys$option[i] %in% c("smin", "smin_noise")      # the paper's rule under s_min
    s <- data.frame(method = keys$method[i], tau = keys$tau[i], option = keys$option[i], reps = nrow(d),
                    proj = if (smin) proj_bhat else proj_pval, proj_pval = proj_pval, proj_bhat = proj_bhat,
                    proj_failed = mean(!is.finite(d$resid)), switched = mean(d$switched),
                    allowance = mean(d$allowance), stringsAsFactors = FALSE)
    for (tn in setdiff(TESTS, "proj")) s[[tn]] <- mean(d[[tn]] <= CFG$alpha, na.rm = TRUE)
    if (CFG$run_ekw) s$EKW_max <- suppressWarnings(max(unlist(s[EKW_EACH]), na.rm = TRUE))
    if (CFG$run_stouffer) s$stouffer_k <- mean(d$stouffer_k, na.rm = TRUE)   # p-values below .05
    s
  }))
}

XLAB       <- expression("Share of trials that switch ("*tau*"); at "*tau == 0*": size")
XLAB_SHORT <- expression("Share of trials that switch ("*tau*")")
OPT_STYLE <- list(
  none       = list(col = "#8C8C8C", lty = "longdash", lwd = 1.3, pch = 4,  label = "No de-rounding"),
  uniform    = list(col = "#E69F00", lty = "dashed",   lwd = 1.3, pch = 2,  label = "Uniform de-rounding"),
  smin       = list(col = "#0072B2", lty = "solid",    lwd = 1.3, pch = 1,  label = "Smallest consistent SE"),
  smin_noise = list(col = "#000000", lty = "solid",    lwd = 1.8, pch = 16, label = "Smallest consistent SE + noise"))
OPT_TEX      <- c(none = "None", uniform = "Unif.", smin = "$s_{\\min}$", smin_noise = "$+$noise")
METHOD_TITLE <- c(max = "Outcome switching: most significant outcome",
                  threshold = sprintf("Threshold: another outcome if |z| <= %s", format(CFG$cv)))
METHOD_TEX   <- c(max = "Most significant outcome",
                  threshold = sprintf("Another outcome if $|z| \\le %s$", format(CFG$cv)))
TEST_TITLE   <- c(proj = "Projection test", T_inf = "Kudrin's T_inf", T_H = "Kudrin's T_H",
                  T_H_hybrid = "Kudrin's T_H hybrid", EKW_min = "Any EKW test (smallest p-value)",
                  EKW_LCM = "EKW: LCM", EKW_Fisher = "EKW: Fisher", EKW_Discontinuity = "EKW: Discontinuity",
                  EKW_Binomial = "EKW: Binomial", EKW_CS1 = "EKW: CS1", EKW_CS2B = "EKW: CS2B",
                  stouffer = "p-curve test for left skew (Stouffer)", EKW_max = "EKW tests (largest rejection rate)")
TEST_TEX     <- c(proj = "Projection test", T_inf = "$T_\\infty$", T_H = "$T_H$",
                  T_H_hybrid = "$T_H^{\\mathrm{hyb}}$", EKW_min = "Any EKW test (smallest $p$-value)",
                  EKW_LCM = "LCM", EKW_Fisher = "Fisher", EKW_Discontinuity = "Discontinuity",
                  EKW_Binomial = "Binomial", EKW_CS1 = "CS1", EKW_CS2B = "CS2B",
                  stouffer = "$p$-curve test for left skew (Stouffer)", EKW_max = "EKW (largest of six)")

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
      mtext(XLAB, side = 1, line = 1.9, cex = 0.85)
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

# The comparison figure, in the style of kudrin_figures_engine.R: every test in one panel per method
STYLE <- list(
  sig        = list(col = "#000000", lty = "dotted",   lwd = 0.8),
  EKW_max    = list(col = "#8C8C8C", lty = "longdash", lwd = 1.1),
  T_inf      = list(col = "#0072B2", lty = "solid",    lwd = 1.1, pch = 2,  cex = 0.62, plwd = 0.9),
  T_H        = list(col = "#E69F00", lty = "dashed",   lwd = 1.1, pch = 3,  cex = 0.62, plwd = 0.9),
  T_H_hybrid = list(col = "#CC79A7", lty = "dotdash",  lwd = 1.1, pch = 5,  cex = 0.62, plwd = 0.9),
  proj       = list(col = "#000000", lty = "solid",    lwd = 1.6, pch = 16, cex = 0.55, plwd = 1.0)
)
AXIS_COL  <- "grey25"
GRID_COL  <- "grey90"
LEG_LABEL <- list(sig        = "5% level",
                  EKW_max    = "Max power of EKW tests",
                  T_inf      = expression(italic(T)[infinity]),
                  T_H        = expression(italic(T)[italic(H)]),
                  T_H_hybrid = expression(italic(T)[italic(H)]^hybrid),
                  proj       = "Projection test")
LEG_ORDER <- c("sig", "EKW_max", "T_inf", "T_H", "T_H_hybrid", "proj")

lighten <- function(col, f = 0.65) {     # opaque tint (works on every device)
  m <- col2rgb(col) / 255
  rgb(t(m + (1 - m) * f))
}

draw_series <- function(x, y, st, faint = FALSE) {
  ok   <- is.finite(y)
  col  <- if (faint) lighten(st$col) else st$col
  lwd  <- if (faint) 0.8 * st$lwd else st$lwd
  lines(x[ok], y[ok], col = col, lty = st$lty, lwd = lwd)
  if (!is.null(st$pch)) {
    plwd <- if (faint) 0.8 * st$plwd else st$plwd
    points(x[ok], y[ok], col = col, pch = st$pch, cex = st$cex, lwd = plwd, xpd = TRUE)
  }
  invisible(NULL)
}

# One panel. r: results for one (S, row, panel); pub: that panel's published curves or NULL.
# The EKW line comes from r when the EKW tests were run, otherwise from pub when available.
draw_panel <- function(r, title, pub = NULL, overlay_published = FALSE, add_projection = TRUE) {
  r <- r[order(r$tau), , drop = FALSE]
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  abline(h = c(0.25, 0.5, 0.75), v = c(0.25, 0.5, 0.75), col = GRID_COL, lwd = 0.5)
  at   <- c(0, 0.25, 0.5, 0.75, 1)
  labs <- c("0", ".25", ".5", ".75", "1")
  for (side in 1:2) {
    axis(side, at = at, labels = labs, lwd = 0, lwd.ticks = 0.5, col.ticks = AXIS_COL,
         col.axis = AXIS_COL, cex.axis = 0.8)
  }
  box(lwd = 0.6, col = AXIS_COL)
  mtext(XLAB_SHORT, side = 1, line = 1.35, cex = 0.78)
  mtext("Rejection rate", side = 2, line = 1.9, cex = 0.78, las = 0)
  mtext(title, side = 3, line = 0.25, cex = 0.8)
  
  abline(h = CFG$alpha, col = STYLE$sig$col, lty = STYLE$sig$lty, lwd = STYLE$sig$lwd)
  if ("EKW_max" %in% names(r)) {
    draw_series(r$tau, r$EKW_max, STYLE$EKW_max)
  } else if (!is.null(pub) && !is.null(pub$EKW_max)) {
    draw_series(pub$tau_ekw, pub$EKW_max, STYLE$EKW_max)
  }
  tests <- c("T_inf", "T_H", "T_H_hybrid")
  if (overlay_published && !is.null(pub)) {
    for (tn in tests) draw_series(pub$tau, pub[[tn]], STYLE[[tn]], faint = TRUE)
  }
  for (tn in tests) {
    if (tn %in% names(r)) draw_series(r$tau, r[[tn]], STYLE[[tn]])
  }
  if (add_projection && ("proj" %in% names(r))) draw_series(r$tau, r$proj, STYLE$proj)
  invisible(NULL)
}

open_device <- function(file, ext, height, width = 6.6) {
  if (ext == "pdf") {
    pdf(file, width = width, height = height, pointsize = 8)
  } else if (isTRUE(capabilities("cairo"))) {
    png(file, width = width, height = height, units = "in", res = 300, pointsize = 8, type = "cairo")
  } else {
    png(file, width = width, height = height, units = "in", res = 300, pointsize = 8)
  }
}

# Legend in a strip across the bottom of the device (drawn last, over the whole page)
draw_legend_strip <- function(keys, ncol) {
  sty  <- STYLE[keys]
  labs <- do.call(c, lapply(keys, function(k) as.expression(LEG_LABEL[[k]])))
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend("bottom", inset = c(0, 0.004), legend = labs, ncol = ncol, bty = "n", xpd = NA,
         col    = vapply(sty, function(s) s$col, ""),
         lty    = vapply(sty, function(s) s$lty, ""),
         lwd    = vapply(sty, function(s) s$lwd, 0),
         pch    = vapply(sty, function(s) if (is.null(s$pch)) NA_real_ else s$pch, 0),
         pt.cex = vapply(sty, function(s) if (is.null(s$cex)) 1 else s$cex, 0),
         pt.lwd = vapply(sty, function(s) if (is.null(s$plwd)) 1 else s$plwd, 0),
         cex = 0.8, seg.len = 2.8, x.intersp = 0.6, y.intersp = 1.05)
}

plot_compare <- function(sm, option, stem, note) {
  ms   <- CFG$hack_methods
  keys <- intersect(LEG_ORDER, c("sig", if (CFG$run_ekw) "EKW_max", if (CFG$run_kudrin) c("T_inf", "T_H", "T_H_hybrid"), "proj"))
  ncol <- 3
  leg_lines <- 0.9 + 1.1 * ceiling(length(keys) / ncol)
  top_lines <- 1.1 * length(note) + 0.4
  line_in   <- 8 * 1.2 / 72
  height    <- 2.35 + (leg_lines + top_lines) * line_in
  width     <- max(3.3, 2.7 * length(ms))
  for (ext in c("pdf", "png")) {
    open_device(paste0(stem, ".", ext), ext, height = height, width = width)
    par(mfrow = c(1, length(ms)))
    par(cex = 1, mar = c(2.6, 3.0, 2.5, 0.8), oma = c(leg_lines, 0, top_lines, 0),
        mgp = c(1.5, 0.35, 0), tcl = -0.22, las = 1, family = "sans")
    for (m in ms) draw_panel(sm[sm$method == m & sm$option == option, , drop = FALSE], METHOD_TITLE[[m]])
    for (j in seq_along(note)) mtext(note[j], side = 3, outer = TRUE, line = 0.2 + 1.1 * (length(note) - j), cex = 0.7)
    draw_legend_strip(keys, ncol)
    dev.off()
  }
  invisible(paste0(stem, c(".pdf", ".png")))
}

# The main results: rows the tests, columns size and power at CFG$main_taus for each method; one block
# of rows per option
write_main_tex <- function(sm, file, caption, notes) {
  ms    <- CFG$hack_methods
  tests <- c("proj", if (CFG$run_kudrin) c("T_inf", "T_H", "T_H_hybrid"), if (CFG$run_ekw) "EKW_max",
             if (CFG$run_stouffer) "stouffer")
  nt    <- length(CFG$main_taus)
  get   <- function(m, o, tau, tn) {
    z <- sm[sm$method == m & sm$option == o & abs(sm$tau - tau) < 1e-9, tn]
    if (length(z) == 0 || !is.finite(z[1])) "--" else sprintf("%.3f", z[1])
  }
  head <- c(paste0(" & Size & ", paste(sprintf("\\multicolumn{%d}{c}{%s}", nt, METHOD_TEX[ms]), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 3 + (seq_along(ms) - 1) * nt, 2 + seq_along(ms) * nt), collapse = " "),
            paste0("Test & $\\tau = 0$ & ", paste(rep(sprintf("$\\tau = %s$", format(CFG$main_taus)), length(ms)), collapse = " & "), " \\\\"))
  body <- character(0)
  for (o in CFG$options) {
    if (length(body)) body <- c(body, "\\midrule")
    body <- c(body, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", 2 + nt * length(ms), OPT_STYLE[[o]]$label))
    body <- c(body, vapply(tests, function(tn) {
      v <- c(get(ms[1], o, 0, tn), unlist(lapply(ms, function(m) vapply(CFG$main_taus, function(tau) get(m, o, tau, tn), ""))))
      paste0(TEST_TEX[[tn]], " & ", paste(v, collapse = " & "), " \\\\")
    }, ""))
  }
  writeLines(c("\\begin{table}[htbp]", "\\centering", sprintf("\\caption{%s}", caption),
               sprintf("\\label{tab:%s}", tools::file_path_sans_ext(basename(file))),
               "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
               sprintf("\\begin{tabular}{l%s}", strrep("c", 1 + nt * length(ms))), "\\toprule", head,
               "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\par\\smallskip",
               sprintf("\\begin{minipage}{\\linewidth}\\footnotesize\\textit{Notes:} %s\\end{minipage}", notes),
               "\\end{table}"), file)
  main <- do.call(rbind, lapply(CFG$options, function(o) do.call(rbind, lapply(tests, function(tn) {
    data.frame(option = o, test = tn, size = sm[sm$method == ms[1] & sm$option == o & sm$tau == 0, tn][1],
               as.list(setNames(unlist(lapply(ms, function(m) vapply(CFG$main_taus, function(tau)
                 sm[sm$method == m & sm$option == o & abs(sm$tau - tau) < 1e-9, tn][1], 0))),
                 unlist(lapply(ms, function(m) paste0(m, "_tau_", CFG$main_taus))))), stringsAsFactors = FALSE)
  }))))
  write.csv(main, sub("\\.tex$", ".csv", file), row.names = FALSE)
  main
}


## ---- 6. Run --------------------------------------------------------------------------------

# The trials, and how much switching changes significance (printed z, no simulation): the share of
# trials with a significant reported outcome when it is the pre-specified one (chosen at random),
# the most significant one, or as under "threshold" (averaged over the choices of the other outcome)
Kt    <- lengths(TRIALS)
sig   <- abs(D$z) > CFG$cv
p_pre <- mean(vapply(TRIALS, function(ix) mean(sig[ix]), 0))
p_max <- mean(vapply(TRIALS, function(ix) any(sig[ix]), 0))
p_thr <- mean(vapply(TRIALS, function(ix) {   #pre-specified j at random; if insignificant, the other outcome
  s <- sig[ix]                                  #(at random among the rest) is reported when it is significant
  mean(vapply(seq_along(s), function(j) if (s[j]) 1 else mean(s[-j]), 0))
}, 0))
report <- c(
  sprintf("Trials: %s randomized trials with at least %d continuous efficacy outcomes of at least %d participants (their first %d in review order), %s outcomes; from %s.",
          format(length(TRIALS), big.mark = ","), CFG$min_outcomes, SAMPLE$settings$n_min, CFG$max_outcomes,
          format(nrow(D), big.mark = ","), basename(CFG$sample_file)),
  sprintf("Outcomes per trial: %s.", paste(sprintf("%s: %s", names(table(Kt)), format(as.vector(table(Kt)), big.mark = ",")), collapse = "; ")),
  sprintf("Meta-samples: %s trials each, drawn without replacement (%.0f%% of the eligible trials); one reported outcome per trial.",
          format(CFG$n_trials, big.mark = ","), 100 * CFG$n_trials / length(TRIALS)),
  sprintf("Share of trials whose reported outcome is significant (|z| > %s, printed z): pre-specified %.3f; most significant %.3f; threshold rule %.3f.",
          format(CFG$cv), p_pre, p_max, p_thr),
  sprintf("Exact zeros among the printed z: %s (%.1f%%), from equal rounded means.", format(sum(D$z == 0), big.mark = ","), 100 * mean(D$z == 0)),
  sprintf("Rounding: kappa = 2a / SE_min; the cap %s leaves out %s of %s outcomes (trials: %s of %s remain eligible). Among the eligible outcomes kappa has median %.3f, 90th percentile %.3f; rounding allowance %s x mean(kappa^6) = %.2e over them (each meta-sample uses its reported outcomes').",
          if (is.null(CFG$kappa_cap)) "(none)" else format(CFG$kappa_cap),
          format(N_BEFORE_CAP[["outcomes"]] - N_AFTER_CAP, big.mark = ","),
          format(N_BEFORE_CAP[["outcomes"]], big.mark = ","), format(length(TRIALS), big.mark = ","),
          format(N_BEFORE_CAP[["trials"]], big.mark = ","), median(D$kappa), quantile(D$kappa, 0.9),
          format(ALLOW_C), ALLOW_C * mean(D$kappa^6)))
message(paste(report, collapse = "\n"))
writeLines(report, file.path(OUT_DIR, paste0(CFG$name, "_sample.txt")))

proj_key  <- paste(CFG$num_coeffs, CFG$L, CFG$numgrid, if (is.null(CFG$epsilon_U_pre)) "computed" else CFG$epsilon_U_pre)
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

SAMPLE_KEY <- paste(normalizePath(CFG$sample_file), file.info(CFG$sample_file)$size, format(file.info(CFG$sample_file)$mtime))
RUN_KEY    <- paste(SAMPLE_KEY, CFG$n_trials, CFG$min_outcomes, CFG$max_outcomes, paste(CFG$taus, collapse = ","),
                    paste(CFG$hack_methods, collapse = ","), CFG$cv, paste(CFG$options, collapse = ","), CFG$proj_mode,
                    if (is.null(CFG$kappa_cap)) "no cap" else CFG$kappa_cap, ALLOW_C, "four options",
                    CFG$proj_boots, CFG$J0, CFG$p_weight, CFG$num_coeffs, CFG$shift, CFG$L, CFG$numgrid,
                    PROJ$epsilon_U_pre, CFG$run_kudrin, CFG$kudrin_boots, CFG$run_ekw, CFG$run_stouffer, CFG$seed)

# Run time: one projection test now, times the test batteries per repetition
t0 <- Sys.time()
invisible(proj_test(D$z[vapply(TRIALS[seq_len(CFG$n_trials)], function(ix) ix[1], 0L)], names(TRIALS)[seq_len(CFG$n_trials)]))
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
per  <- length(CFG$options) * (1 + length(CFG$hack_methods) * sum(CFG$taus > 0))
message(sprintf("One projection test takes %.1f s; %d per repetition; %d repetitions on %d cores: about %.1f hours (more with Kudrin or EKW).",
                secs, as.integer(per), as.integer(CFG$reps), as.integer(CFG$n_cores), secs * per * CFG$reps / CFG$n_cores / 3600))

cl <- NULL
if (CFG$n_cores > 1) {
  cl <- parallel::makeCluster(CFG$n_cores)
  parallel::clusterExport(cl, "CODE_FILES", envir = globalenv())
  parallel::clusterEvalQ(cl, { for (f in CODE_FILES) source(f); NULL })
  parallel::clusterExport(cl, c("CFG", "OUT_DIR", "D", "TRIALS", "PROJ", "RUN_KEY", "run_ewk", "EKW_TESTS", "z_at",
                                "z_option", "ALLOW_C",
                                "proj_test", "kudrin_p", "ekw_p", "stouffer_p", "run_tests", "one_rep", "rep_file",
                                "one_rep_cached"), envir = globalenv())
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
print(SM[, c("method", "option", "tau", "reps", "switched", intersect(c(FIG_TESTS, "EKW_max"), names(SM)))], digits = 3, row.names = FALSE)

# Two lines above every figure: the trials, then the selection
note <- c(sprintf("Cochrane: continuous efficacy outcomes of randomized trials (n >= %d); meta-samples of %s of %s trials with >= %d outcomes; %d repetitions; shaded: 95%% band around %s%%.",
                  SAMPLE$settings$n_min, format(CFG$n_trials, big.mark = ","), format(length(TRIALS), big.mark = ","),
                  CFG$min_outcomes, as.integer(CFG$reps), format(100 * CFG$alpha)),
          sprintf("A share tau of trials reports its most significant outcome (max) or switches when |z| <= %s (threshold); the rest report a random one.",
                  format(CFG$cv)))
RULE_FIG <- sprintf("Under s_min and s_min + noise the projection test rejects only if B-hat > rounding allowance %s x mean(kappa^6) of the reported outcomes; otherwise if p <= %s.",
                    format(ALLOW_C), format(CFG$alpha))
for (tn in FIG_TESTS) plot_power(SM, tn, file.path(OUT_DIR, paste0(CFG$name, "_", tn)), if (tn == "proj") c(note, RULE_FIG) else note)
for (o in CFG$options) {
  plot_compare(SM, o, file.path(OUT_DIR, paste0(CFG$name, "_compare_", o)),
               c(note[1], sprintf("%s; %s", OPT_STYLE[[o]]$label, note[2]),
                 if (o %in% c("smin", "smin_noise")) RULE_FIG))
}

big <- function(x) gsub(",", "{,}", format(x, big.mark = ","), fixed = TRUE)
NOTES_ALL <- paste0(sprintf("Rejection rates of %s\\%%-level tests over %d repetitions. Each meta-sample draws %s of the %s randomized trials in the Cochrane Database of Systematic Reviews (Schwab's data set, as in van Zwet, Schwab and Senn 2021) with at least %d continuous efficacy outcomes of at least %d participants, without replacement, and reports one outcome per trial. ",
                            format(100 * CFG$alpha), as.integer(CFG$reps), big(CFG$n_trials), big(length(TRIALS)),
                            CFG$min_outcomes, SAMPLE$settings$n_min),
                    sprintf("A trial chooses among its first %d outcomes in review order. With probability $1 - \\tau$ it reports one chosen at random (its pre-specified outcome); with probability $\\tau$ it reports its most significant outcome (max) or, when the pre-specified outcome has $|z| \\le %s$, the larger of it and one other outcome (threshold). ",
                            CFG$max_outcomes, format(CFG$cv)),
                    "Trials select on unrounded $z$: each repetition draws the means and standard deviations uniformly within their rounding cells. Columns: the analyst's treatment of rounding. None: the printed $z$. Unif.: means and standard deviations drawn uniformly within their rounding cells. $s_{\\min}$: the mean difference over $SE_{\\min}$, the standard error from the smallest standard deviations consistent with the printed ones. $+$noise: the same with uniform noise at each mean's last digit. ",
                    sprintf("Outcomes whose rounding unit in $z$-units, $\\kappa = 2a/SE_{\\min}$ ($a$ the unit of the means; a difference of two rounded means can be off by $a$), exceeds %s are left out under every option. Every test treats trials as independent units.",
                            if (is.null(CFG$kappa_cap)) "(no cap)" else format(CFG$kappa_cap)))
like_coch <- function(what) sprintf("%s against outcome switching in Cochrane trials", what)
TABLES <- list(proj = list(tests = "proj", caption = like_coch("Projection test"),
                           notes = sprintf("Projection test: run\\_test\\_weight\\_fast() with $J = 40$, $J_0 = 20$ and $p = 2$. Under $s_{\\min}$ and $+$noise it rejects only when its breakdown statistic $\\hat B$ exceeds the rounding allowance %s times the mean of $\\kappa^6$ over the reported outcomes (Propositions 3 and 4); under None and Unif.\\ when $p \\le %s$.",
                                           format(ALLOW_C), format(CFG$alpha))))
if (CFG$run_stouffer) TABLES$stouffer <- list(tests = "stouffer", caption = like_coch("The $p$-curve test for left skew (Stouffer)"),
                                              notes = "$p$-curve test for left skew, as in simulations\\_brodeur\\_rct.R.")
if (CFG$run_ekw) TABLES$EKW <- list(tests = names(EKW_TESTS), caption = like_coch("Each EKW test"),
                                    notes = "EKW tests as in the application, on the two-sided $p$-values below 0.15.")
if (CFG$run_kudrin) TABLES$kudrin <- list(tests = c("T_inf", "T_H", "T_H_hybrid"), caption = like_coch("Kudrin's (2025) tests"),
                                          notes = "Kudrin (2025) tests $T_\\infty$, $T_H$ and $T_H^{\\mathrm{hyb}}$.")
for (nm in names(TABLES)) {
  write_rates_tex(SM, TABLES[[nm]]$tests, file.path(OUT_DIR, paste0(CFG$name, "_", nm, ".tex")),
                  TABLES[[nm]]$caption, paste(NOTES_ALL, TABLES[[nm]]$notes))
}
MAIN <- write_main_tex(SM, file.path(OUT_DIR, paste0(CFG$name, "_main.tex")),
                       "Size and power against outcome switching in Cochrane trials",
                       paste(NOTES_ALL, "Size: rejection rate at $\\tau = 0$, where every trial reports its pre-specified outcome. EKW: the largest rejection rate of its six tests at each $\\tau$, as in Kudrin's (2025) figures."))
message("Main results (rejection rates):")
print(MAIN, digits = 3, row.names = FALSE)
message("Outputs in ", OUT_DIR)