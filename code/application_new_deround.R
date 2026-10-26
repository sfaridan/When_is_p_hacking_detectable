# 1. Clear all objects from the global environment
rm(list = ls())

# 2. Force R to release the freed memory to your OS
gc(reset = TRUE)

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_all_weighted.R"))

#Load in the data
data_location <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data"
setwd(data_location)
MM_raw <- haven::read_dta('MM data.dta')
#MM_raw <- haven::read_dta('MM data fromkudrin.dta')
MM_data <- MM_raw

#Kudrin's (2025) drop
drop          <- is.na(MM_raw$mu) | MM_raw$mu == 0
MM_dropped    <- MM_raw[!drop, ]
MM_data       <- MM_dropped

#### Prep data: de-round the t-statistics (Online Appendix D, Propositions D.1-D.2)
# Four de-rounding options, chosen with deround_option below. For a t-statistic computed from a coefficient
# m and a standard error s printed to multiples of a and b, with s_min = s - b/2 and s_max = s + b/2 the
# smallest and largest standard errors consistent with s:
#   1 "none":       t = m / s
#   2 "uniform":    t = (m + V) / W,  V ~ U(-a/2, a/2), W ~ U(s_min, s_max)  (the usual de-rounding)
#   3 "smin":       t = m / s_min                                              (Appendix D.1)
#   4 "smin_noise": t = (m + V) / s_min
# A reported t-statistic is used as printed under options 1 and 3, and gets noise uniform on +-1/2 of its
# last digit under options 2 and 4. Options 2 and 4 are run on n_draws independent draws, and p-values and
# breakdown statistics are averaged over them. Not covered: t-statistics converted from rounded p-values.
#
# The sample is the same under every option: which t-statistics are dropped depends only on the printed
# values and kappa_cap, never on deround_option. The noise is drawn once, with fixed seeds, so options 2
# and 4 use identical numerator noise and differ only in the denominator.
#
# Rounding units:
#   b, per entry: the place value of sd's last digit. A dropped trailing zero only makes b too large,
#     which keeps s_min a lower bound; a table-wide b can be too small for a column printed more coarsely.
#   a, per table: the finest place value among the table's coefficients, but no finer than b, since a
#     coefficient and its standard error are almost always printed to the same decimals. Both are fixed
#     by the table's format and the standard error, not by the estimate. The width of V need not match
#     the true unit: Proposition D.1 allows any independent noise.

deround_option <- 4     #1 = none, 2 = uniform (standard), 3 = smallest SE, 4 = smallest SE + noise in the numerator
kappa_cap      <- 1     #drop coef/SE t-statistics with kappa = a/s_min above this (the same under every option)
se_rounded_up  <- FALSE #TRUE if papers may have rounded standard errors up: then s_min = sd - b and s_max = sd
n_draws        <- 100    #independent de-roundings for options 2 and 4; try 2 first to check the run time

stopifnot(deround_option %in% 1:4)
deround_spec <- list(key   = c("none", "uniform", "smin", "smin_noise")[deround_option],
                     label = c("No de-rounding", "Uniform de-rounding", "Smallest consistent standard error",
                               "Smallest consistent standard error with uniform noise")[deround_option],
                     noise = c(FALSE, TRUE, FALSE, TRUE)[deround_option])
add_noise <- deround_spec$noise
suffix    <- paste0("_", deround_spec$key)                   #labels every results file and the table
cat("De-rounding option ", deround_option, ": ", deround_spec$label, "\n", sep = "")

n_dec <- function(x, max_dec = 12) {   #decimals as printed; 6 significant digits clear float storage noise
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

n_obs  <- nrow(MM_data)
dec_mu <- n_dec(MM_data$mu)
dec_sd <- n_dec(MM_data$sd)
dec_t  <- n_dec(MM_data$t)
mu_p   <- round(MM_data$mu, dec_mu)   #printed values without float storage noise, which could
sd_p   <- round(MM_data$sd, dec_sd)   #otherwise put s_min a hair above the true SE
t_p    <- round(MM_data$t,  dec_t)
tab_id <- paste(MM_data$title, MM_data$table, sep = " | ")   #a table = an article x table number
b_sd   <- 10^(-dec_sd)
a_tab  <- 10^(-ave(dec_mu, tab_id, FUN = function(d) if (all(is.na(d))) NA_integer_ else max(d, na.rm = TRUE)))
a_mu   <- pmax(a_tab, b_sd)                                  #NA unless both mu and sd are present
s_min  <- sd_p - if (se_rounded_up) b_sd else b_sd / 2       #smallest SE consistent with sd
s_max  <- sd_p + if (se_rounded_up) 0 else b_sd / 2          #largest SE consistent with sd
a_t    <- 10^(-dec_t)                                        #rounding unit of a directly reported t

has_ratio  <- !is.na(mu_p) & !is.na(sd_p)
drop_se0   <- has_ratio & sd_p <= 0                          #SE printed as zero: t not recoverable
drop_kappa <- has_ratio & !drop_se0 & a_mu / s_min > kappa_cap
use_ratio  <- has_ratio & !drop_se0 & !drop_kappa
direct     <- !has_ratio & !is.na(t_p)

cat("Decimals printed in the standard errors (many at 6+ would suggest computed rather than printed SEs):\n")
print(table(dec_sd[has_ratio]))
cat("\nt-statistics dropped (the same under every de-rounding option):\n")
print(rbind(se_zero         = tapply(drop_se0,   MM_data$method, sum),
            kappa_above_cap = tapply(drop_kappa, MM_data$method, sum)))

#### How much do we drop, compared with Kranz and Putz's (2022) omission rule?
# Runs before MM_data is subset, while the vectors above still line up with the full sample.
#
# KP omit every test whose standard error has a significand below 37: the SE's significant digits read
# as an integer, e.g. 0.012 -> 12. With the coefficient and SE printed to the same decimals, this keeps
# rounding from moving any z-statistic across 1.96 unless it lies within 0.04 of it (their Lemma 1).
# Our data have lost trailing zeros, so we use KP's own fallback for that case: append zeros to the SE
# until it has as many decimals as its coefficient (KP online appendix, note to Figure A2), which they
# report agrees with the exact rule in 96 percent of cases. If your file records trailing zeros, read
# the significand from those instead. Reported t-statistics (no coefficient and SE) are kept by both rules.

kp_threshold <- 37
kp_signif    <- round(sd_p * 10^pmax(dec_sd, dec_mu))                   #SE significand, KP's convention
drop_kp      <- has_ratio & kp_signif < kp_threshold                     #includes SEs printed as zero
drop_ours    <- drop_se0 | drop_kappa
z_is_2       <- has_ratio & sd_p > 0 & abs(abs(mu_p / sd_p) - 2) < 1e-9   #KP's heap at exactly z = 2

drop_summary <- function(g) {                                            #g: logical row selector
  c(t_stats     = sum(g),
    ours        = sum(drop_ours[g]),            kp     = sum(drop_kp[g]),
    ours_pct    = 100 * mean(drop_ours[g]),     kp_pct = 100 * mean(drop_kp[g]),
    both        = sum((drop_ours & drop_kp)[g]),
    only_ours   = sum((drop_ours & !drop_kp)[g]),
    only_kp     = sum((!drop_ours & drop_kp)[g]),
    z_eq_2      = sum(z_is_2[g]),
    ours_z2_pct = 100 * mean(drop_ours[g & z_is_2]),
    kp_z2_pct   = 100 * mean(drop_kp[g & z_is_2]))
}
methods_kp      <- c("RCT", "IV", "DID", "RDD")
rows_kp         <- c(setNames(lapply(methods_kp, function(m) MM_data$method %in% m), methods_kp),
                     list(Pooled = rep(TRUE, nrow(MM_data))))
drop_comparison <- t(sapply(rows_kp, drop_summary))

cat("\nt-statistics dropped by our procedure (SE = 0 or kappa > ", kappa_cap,
    ") and by Kranz and Putz's rule (SE significand < ", kp_threshold, "):\n", sep = "")
print(round(drop_comparison, 1))
cat(sprintf("Our procedure drops %.1f to %.1f percent of the t-statistics by method (%.1f pooled); Kranz and Putz's rule would drop %.1f to %.1f percent (%.1f pooled).\n",
            min(drop_comparison[methods_kp, "ours_pct"]), max(drop_comparison[methods_kp, "ours_pct"]), drop_comparison["Pooled", "ours_pct"],
            min(drop_comparison[methods_kp, "kp_pct"]),   max(drop_comparison[methods_kp, "kp_pct"]),   drop_comparison["Pooled", "kp_pct"]))
cat(sprintf("Of the %d t-statistics at exactly z = 2, our procedure drops %.1f percent and theirs %.1f percent. (KP report 37.9 percent of all observations and 87.3 percent of those at z = 2, on their version of the data.)\n",
            drop_comparison["Pooled", "z_eq_2"], drop_comparison["Pooled", "ours_z2_pct"], drop_comparison["Pooled", "kp_z2_pct"]))
saveRDS(drop_comparison, file = file.path(root, "results", "Drop_comparison_KP"))

#Keep the printed values each option needs, then drop (the drops above do not depend on the option)
MM_data$t_source   <- ifelse(use_ratio, "coef/SE", ifelse(direct, "reported t", NA))
MM_data$dr_mu      <- mu_p                                                       #printed coefficient
MM_data$dr_sd      <- sd_p                                                       #printed standard error
MM_data$dr_smin    <- s_min                                                      #smallest SE consistent with it
MM_data$dr_smax    <- s_max                                                      #largest SE consistent with it
MM_data$dr_t       <- t_p                                                        #printed t-statistic
MM_data$dr_a       <- ifelse(use_ratio, a_mu, ifelse(direct, a_t, NA))           #rounding unit of the numerator
MM_data$kappa      <- ifelse(use_ratio, a_mu / s_min, ifelse(direct, a_t, NA))   #rounding unit in t-units
MM_data$lambda_max <- ifelse(use_ratio, s_max / s_min, NA)                       #largest inflation from dividing by s_min
MM_data            <- MM_data[!(drop_se0 | drop_kappa), ]

#The de-rounding noise, drawn once with fixed seeds, so every option and every run uses the same draws
if (!add_noise) n_draws <- 1                                 #without noise every draw would be the same
set.seed(1)
U_num <- matrix(runif(nrow(MM_data) * n_draws, -0.5, 0.5), ncol = n_draws)   #numerator noise, one column per draw
set.seed(2)
U_den <- matrix(runif(nrow(MM_data) * n_draws, -0.5, 0.5), ncol = n_draws)   #denominator noise (option 2 only)

deround_draw <- function(r) {                                #t-statistics for draw r under the chosen option
  d     <- MM_data
  ratio <- d$t_source %in% "coef/SE"
  u_num <- if (add_noise) U_num[, r] else 0
  u_den <- if (add_noise) U_den[, r] else 0
  num   <- d$dr_mu + d$dr_a * u_num                          #coefficient, plus noise under options 2 and 4
  den   <- switch(deround_spec$key,
                  none       = d$dr_sd,                                                       #printed SE
                  uniform    = (d$dr_smin + d$dr_smax) / 2 + (d$dr_smax - d$dr_smin) * u_den,  #uniform on [s_min, s_max]
                  smin       = d$dr_smin,
                  smin_noise = d$dr_smin)
  ifelse(ratio, num / den, d$dr_t + d$dr_a * u_num)          #a reported t: as printed, plus noise under options 2 and 4
}

if (any(!is.finite(deround_draw(1))) | any(!is.finite(MM_data$mu))){
  stop("Tscores not clean")
}

cat("\nt-statistics by source:\n")
print(table(MM_data$method, MM_data$t_source, useNA = "ifany"))
if (deround_spec$key %in% c("smin", "smin_noise")) {         #these diagnostics concern dividing by s_min
  cat("\nLargest possible inflation of t, s_max/s_min (median, 99th percentile, max):\n")
  print(round(sapply(split(MM_data$lambda_max, MM_data$method), quantile, probs = c(.5, .99, 1), na.rm = TRUE), 3))
  kappa6 <- tapply(MM_data$kappa^6, MM_data$method, mean, na.rm = TRUE)
  cat("\nRounding allowance (Appendix D.1, step 3), to compare with B-hat:\n")
  print(signif(rbind(mean_kappa6 = kappa6, bound_J0_20 = 0.13 * kappa6, bound_J0_30 = 0.31 * kappa6), 3))
} else {
  cat("\nNo rounding allowance is reported: Appendix D's bound covers only options 3 and 4.\n")
}

#Averaging over draws
meth            <- c("RCT","IV","DID","RDD")
as_num1         <- function(v) as.numeric(v)[1]              #test results may come back as lists or list-columns
draw_range      <- function(x) diff(range(x))                #max - min across draws
summarize_draws <- function(d, cols, f = mean) {             #one row per method, in the order of meth;
  out <- d[match(meth, d$methods), setdiff(names(d), "draw")] #cols summarized over draws by f
  for (k in cols) out[[k]] <- as.numeric(tapply(d[[k]], d$methods, f)[meth])
  rownames(out) <- NULL
  out
}

### Run the tests

cv         <- 1.96
sigma_Y    <- 1 
shift      <- cv
L          <- 6.5 #grid support [-L,L] \cup 9999
numgrid    <- 3000 #total number of elements in grid (density)
boots      <- 1000 #bootstrap repetitions
J          <- 40
coeffs     <- J+1
p          <- 2
J0_list    <- c(20) #one (p, B-hat) column group per J0 in the table, e.g. c(30, 20)
epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U #pre-compute for speed

#Run the projection test with more coefficients, on each de-rounding

proj_all <- NULL
for (r in seq_len(n_draws)) {
  t_start   <- Sys.time()
  MM_data$t <- deround_draw(r)
  set.seed(1000 * r + 1) #bootstrap seeds differ across draws, so bootstrap noise averages out too
  for (J0 in J0_list) {
    for (m in meth) {
      res <- run_test_weight_fast(MM_data[MM_data$method==m,],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
      proj_all <- rbind(proj_all, data.frame(draw = r, methods = m, n = as_num1(res$n), articles = as_num1(res$articles),
                                             pval = as_num1(res$pval), breakdown = as_num1(res$breakdown),
                                             numcoeffs = as_num1(res$numcoeffs), J0 = J0))
    }
  }
  cat(sprintf("Projection test (%s): de-rounding %d of %d done in %.1f minutes\n", deround_spec$key, r, n_draws,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

#Save the averages over draws (the files make_application_table() reads) and every draw, labeled by option
setwd(paste0(root,"/results"))
saveRDS(proj_all, file = paste0("MM_weights_results_draws", suffix))
for (J0 in J0_list) {
  name <- paste0("MM_weights_results_",J0,"_J0", suffix)
  saveRDS(summarize_draws(proj_all[proj_all$J0 == J0, ], c("pval", "breakdown")), file = name)
  print(readRDS(name))
}

#Run EWK Tests, on each de-rounding

source(paste0(root,"/code/Tests.R"))

run_ewk <- function(ps,ids=1:length(ps),pmax=0.15){
  pmin <- min(ps[ps>0])
  out <- list()
  
  out$lcms_EWK     <- LCM(ps,pmin,pmax)
  out$Fisher_EWK   <- Fisher(ps,pmin,pmax)
  out$disconts_EWK <- Discontinuity_test(ps,.05)
  out$binomial_EWK <- Binomial(ps,.04, .05, "c")
  out$CS1_EWK_s    <- CoxShi(ps, ids, pmin, pmax, 15, 1, 0)
  out$CS2_EWK_s    <- CoxShi(ps, ids, pmin, pmax, 15, 2, 1)
  
  out$min <- min(out$lcms_EWK,out$Fisher_EWK,out$disconts_EWK, out$binomial_EWK,  out$CS1_EWK_s, out$CS2_EWK_s)
  
  return(out)
}

ewk_draws <- vector("list", n_draws)
for (r in seq_len(n_draws)) {
  MM_data$t <- deround_draw(r)
  set.seed(1000 * r + 2)
  ewk_draws[[r]] <- sapply(setNames(meth, paste0("ewk_", tolower(meth))), function(m) {   #tests x methods
    s <- MM_data$method == m
    vapply(run_ewk(2 * pnorm(-abs(MM_data$t[s])), ids = as.integer(factor(MM_data$title[s]))), as_num1, numeric(1))
  })
}
ewk_all       <- simplify2array(ewk_draws)                   #tests x methods x draws
empirical_ewk <- apply(ewk_all, c(1, 2), mean)
ewk_tests     <- setdiff(rownames(empirical_ewk), "min")
empirical_ewk["min", ] <- apply(empirical_ewk[ewk_tests, , drop = FALSE], 2, min)   #smallest averaged p-value
setwd(paste0(root,"/results"))
saveRDS(ewk_all, file = paste0("EWK_results_draws", suffix))
saveRDS(empirical_ewk, file = paste0("EWK_results", suffix))
print(round(empirical_ewk, 3))

#Run Kudrin (2025) Tests: T_inf, T_H, T_H_hybrid, on each de-rounding

source(paste0(root,"/code/kudrin_tests.R"))
#kudrin_tests_selfcheck() #optional one-time check of kudrin_tests.R

run_kudrin <- function(data){
  kt <- kudrin_tests(data$t, cluster = data$title, boots = boots) #cluster by article, as in the projection test
  data.frame(n = length(data$t), articles = length(unique(data$title)), as.list(kt$pval),
             min = min(kt$pval), restricted = kt$restricted) #restricted = FALSE means T_H_hybrid = T_H
}
run_kudrin_nocluster <- function(data){
  kt <- kudrin_tests(data$t, cluster = NULL, boots = boots) #no clustering, to diagnose the gap with Kudrin's (2025) Table 2
  data.frame(n = length(data$t), articles = length(unique(data$title)), as.list(kt$pval),
             min = min(kt$pval), restricted = kt$restricted) #restricted = FALSE means T_H_hybrid = T_H
}

kud_all <- NULL
for (r in seq_len(n_draws)) {
  t_start   <- Sys.time()
  MM_data$t <- deround_draw(r)
  set.seed(1000 * r + 3) #T_inf is bootstrapped
  kud_all <- rbind(kud_all, cbind(draw = r, methods = meth,
                                  do.call(rbind, lapply(meth, function(m) run_kudrin(MM_data[MM_data$method == m, ])))))
  cat(sprintf("Kudrin tests (%s): de-rounding %d of %d done in %.1f minutes\n", deround_spec$key, r, n_draws,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
kud_tests        <- setdiff(names(kud_all), c("draw", "methods", "n", "articles", "min", "restricted"))
empirical_kudrin <- summarize_draws(kud_all, kud_tests)
empirical_kudrin$min        <- apply(empirical_kudrin[kud_tests], 1, min)   #smallest averaged p-value
empirical_kudrin$restricted <- as.logical(tapply(kud_all$restricted, kud_all$methods, any)[meth])   #bound in any draw
setwd(paste0(root,"/results"))
saveRDS(kud_all, file = paste0("Kudrin_results_draws", suffix))
saveRDS(empirical_kudrin, file = paste0("Kudrin_results", suffix))
print(readRDS(paste0("Kudrin_results", suffix)))

#How much do the results move across de-roundings? (max - min over draws)
if (n_draws > 1) {
  cat("\nSpread across de-roundings (max - min):\n")
  for (J0 in J0_list) {
    cat("Projection test, J0 =", J0, "\n")
    print(summarize_draws(proj_all[proj_all$J0 == J0, ], c("pval", "breakdown"), draw_range)[, c("methods", "pval", "breakdown")])
  }
  cat("EKW tests\n")
  print(round(apply(ewk_all, c(1, 2), draw_range), 3))
  cat("Kudrin (2025) tests\n")
  print(summarize_draws(kud_all, kud_tests, draw_range)[, c("methods", kud_tests)])
}

#### Table: projection test, EKW and Kudrin (2025) tests (same layout as Table 3 in the paper)
# Reads the results files saved above for one de-rounding option (so after a full run you can tweak and
# re-run this block alone), prints the table and writes results/table_application<suffix>.tex.
# LaTeX needs booktabs, graphicx (shrink-to-fit wrapper) and hyperref or url (for \url in the notes).

make_application_table <- function(J0_values = c(20),            # one (p, B-hat) group per J0, in this order; each must be in J0_list above
                                   n_draws   = 1,                # number of de-roundings the results files average over
                                   spec      = deround_spec,     # the de-rounding option the results files were run with
                                   suffix    = paste0("_", spec$key),
                                   res_dir   = paste0(root, "/results"),
                                   out_file  = paste0("table_application", suffix, ".tex")) {
  meth <- c("RCT", "IV", "DID", "RDD")
  
  # Effective sample sizes from the hand-collected sample of 100 articles per method (typed in from the paper)
  med_nu   <- c(RCT = 185,    IV = 262,    DID = 293.5,  RDD = NA)
  Delta_nu <- c(RCT = 0.0014, IV = 0.0013, DID = 0.0014, RDD = NA)
  
  # Citations used in the header and notes; swap in "\\citet{key}" to get hyperlinked citations
  cite <- c(brodeur = "Brodeur et al. (2020)", ekw  = "Elliott et al. (2022)", kudrin = "Kudrin (2025)",
            coxshi  = "Cox and Shi (2023)",    osqp = "Stellato et al. (2024)", brodeur16 = "Brodeur et al. (2016)")
  
  num  <- function(x) vapply(x, function(v) as.numeric(v)[1], numeric(1), USE.NAMES = FALSE) #results files hold list-columns
  pick <- function(d, what) {                                    #rows of d in the order of meth
    i <- match(meth, d$methods)
    if (anyNA(i)) stop(what, " has no results for ", paste(meth[is.na(i)], collapse = ", "))
    d[i, ]
  }
  
  #--- Read the saved results for this option
  res_file <- function(stem) {
    f <- file.path(res_dir, paste0(stem, suffix))
    if (!file.exists(f)) stop("No results file ", basename(f), ": run the tests with this de-rounding option first")
    f
  }
  proj <- lapply(J0_values, function(j0) {
    f <- res_file(paste0("MM_weights_results_", j0, "_J0"))
    r <- pick(readRDS(f), basename(f))
    if (any(num(r$numcoeffs) != coeffs)) warning(basename(f), " was run with numcoeffs != ", coeffs, " (stale file?)")
    data.frame(n = num(r$n), articles = num(r$articles), p = num(r$pval), B = num(r$breakdown))
  })
  ewk <- readRDS(res_file("EWK_results"))                        #matrix: tests x methods
  ewk <- matrix(num(ewk), nrow(ewk), dimnames = dimnames(ewk))[, paste0("ewk_", tolower(meth))]
  kud <- pick(readRDS(res_file("Kudrin_results")), "Kudrin_results")
  
  n_all <- cbind(sapply(proj, function(d) d$n),        num(kud$n))   #same samples in every file?
  a_all <- cbind(sapply(proj, function(d) d$articles), num(kud$articles))
  if (any(n_all != n_all[, 1]) || any(a_all != a_all[, 1]))
    warning("n or article counts differ across results files: were they all run on the same MM_data?")
  
  #--- Numbers
  tab <- data.frame(n = proj[[1]]$n, articles = proj[[1]]$articles,
                    med_nu = med_nu[meth], Delta_nu = Delta_nu[meth], row.names = meth)
  for (g in seq_along(J0_values)) {
    tab[[paste0("p_J0_", J0_values[g])]] <- proj[[g]]$p
    tab[[paste0("B_J0_", J0_values[g])]] <- proj[[g]]$B
  }
  ekw_rows <- c(LCM = "lcms_EWK", Fisher = "Fisher_EWK", Disc = "disconts_EWK", Binomial = "binomial_EWK",
                CS1 = "CS1_EWK_s", CS2B = "CS2_EWK_s", EKW_min = "min")
  for (k in names(ekw_rows)) tab[[k]] <- unname(ewk[ekw_rows[[k]], ])
  kud_cols <- c(T_inf = "T_inf", T_H = "T_H", T_H_hybrid = "T_H_hybrid", Kudrin_min = "min")
  for (k in names(kud_cols)) tab[[k]] <- num(kud[[kud_cols[[k]]]])
  tab$hybrid_restricted <- as.logical(kud$restricted)            #FALSE means T_H_hybrid = T_H
  
  #--- Cells in the paper's format: .59, 1.00, .0034, 0.0014; missing = centered dot
  fmt <- function(x, digits, lead0 = FALSE) {
    s <- sprintf(paste0("%.", digits, "f"), x)
    if (!lead0) s <- sub("^0\\.", ".", s)
    ifelse(is.na(x), "$\\cdot$", s)
  }
  cells <- cbind(meth, sprintf("%.0f", tab$n), sprintf("%.0f", tab$articles),
                 ifelse(is.na(tab$med_nu), "$\\cdot$", as.character(tab$med_nu)), fmt(tab$Delta_nu, 4, lead0 = TRUE))
  for (j0 in J0_values) cells <- cbind(cells, fmt(tab[[paste0("p_J0_", j0)]], 2), fmt(tab[[paste0("B_J0_", j0)]], 4))
  for (k in c(names(ekw_rows), names(kud_cols))) cells <- cbind(cells, fmt(tab[[k]], 2))
  unr <- !tab$hybrid_restricted                                  #dagger: hybrid restriction did not bind
  if (any(unr)) cells[unr, ncol(cells) - 1] <- paste0(cells[unr, ncol(cells) - 1], "$^{\\dagger}$")
  
  #--- Notes (one sentence per line in the .tex; edit freely)
  deround_notes <- c(
    none       = "De-rounding: none. Each $t$-statistic is the ratio of the reported coefficient and standard error, and reported $t$-statistics are used as printed.",
    uniform    = paste0("De-rounding: uniform, following ", cite[["brodeur16"]], ". The coefficient and the standard error are each replaced by an independent uniform draw from its rounding interval, ",
                        "and a reported $t$-statistic gets uniform noise at its last digit."),
    smin       = "De-rounding: smallest consistent standard error (Appendix D.1). Each coefficient is divided by the smallest standard error consistent with the reported one, and reported $t$-statistics are used as printed.",
    smin_noise = paste0("De-rounding: smallest consistent standard error with uniform noise (Appendix D.1). Each coefficient, plus uniform noise at its last digit, is divided by the smallest standard error ",
                        "consistent with the reported one, and a reported $t$-statistic gets uniform noise at its last digit."))
  rdd_insig <- isTRUE(all(unlist(tab["RDD", paste0("p_J0_", J0_values)]) > 0.05))
  if (!rdd_insig) warning("RDD projection-test p-value is now <= .05, so the notes' reason for omitting its Delta_nu no longer holds")
  notes <- c(
    paste0("Results for the ", cite[["brodeur"]], " data. Following ", cite[["kudrin"]], ", we drop tests whose reported coefficient is missing or zero."),
    deround_notes[[spec$key]],
    sprintf("The sample is the same for every de-rounding choice: we also drop $t$-statistics whose standard error is reported as zero or whose coefficient's rounding unit exceeds %s times the smallest standard error consistent with the reported one.", kappa_cap),
    if (n_draws > 1) sprintf("All $p$-values and breakdown statistics are averages over %d independent de-roundings.", n_draws),
    paste0("Columns headed $J_0$ report $p$-values and breakdown statistics $\\widehat{B}$ of the projection test with $J=", J, "$",
           if (length(J0_values) > 1) "; several values of $J_0$ are included to show low sensitivity to the smoothing parameter", "."),
    "Med.\\ $\\nu$ is the median effective sample size as measured by the median number of clusters in a random sample of 100 articles for each causal inference method.",
    "In ambiguous cases, the smallest cluster count was always used.",
    paste0("$\\Delta_\\nu=\\frac{1}{100}\\sum_{k=1}^{100}\\sup_h\\|\\varphi(t-h)-t_{\\nu_k,h}(t)\\|_T$ is the average difference between the PDFs of the standard normal ",
           "and Student-$t$ over the sample of 100 effective sample sizes collected (computed numerically)."),
    if (rdd_insig) "Since the $p$-value for RDDs is insignificant, we did not compute $\\Delta_\\nu$ for them.",
    sprintf("We form $\\mathcal{X}$ using a grid of %d elements evenly spaced between $\\pm %s$ plus the scalar 9999.", numgrid, L),
    sprintf("Discretization incurs error $\\mathcal{E}(L,\\delta)=%.5f$, which is already accounted for in $\\widehat{B}$ and the $p$-value.", epsilon_U_pre),
    "The discretization was chosen to be made of round numbers and result in error relatively small compared to $\\Delta_\\nu$.",
    paste0("EKW $p$-values are for the tests of ", cite[["ekw"]], ": the least concave majorant test (LCM) and Fisher's test (Fisher) on $p$-values in $(0,0.15]$, ",
           "the density discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$ (Binomial), and the ", cite[["coxshi"]],
           " histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) and for 2-monotonicity with bounds (CS2B)."),
    paste0(cite[["kudrin"]], " $p$-values are for his KS-type test ($T_\\infty$), histogram-based test ($T_H$) and hybrid histogram-based test ($T_H^{\\mathrm{hybrid}}$), ",
           "using the sinc kernel with $\\kappa=3$ and, for $T_H$ and $T_H^{\\mathrm{hybrid}}$, the points $0.01,0.02,\\ldots,0.15$."),
    "All three cluster by article: $T_H$ and $T_H^{\\mathrm{hybrid}}$ use cluster-robust covariance matrices and $T_\\infty$ uses an article-level bootstrap.",
    if (any(unr)) "$^\\dagger$The deconvolution estimate already satisfies the hybrid test's shape restrictions, so $T_H^{\\mathrm{hybrid}}=T_H$.",
    if (n_draws > 1) "Min is the smallest of the averaged $p$-values in each group of tests." else "Min is the smallest $p$-value in each group of tests.",
    sprintf("The projection test and $T_\\infty$ use %d bootstrap repetitions.", boots),
    paste0("Computation used osqp package in R by ", cite[["osqp"]], ". Data: ", cite[["brodeur"]], "."),
    "Transparency R package available at: \\url{https://github.com/sfaridan/When_is_p_hacking_detectable}")
  
  #--- LaTeX
  width <- c(rep(2, length(J0_values)), 7, 4)                    #column groups after the 5 leading columns
  last  <- 5 + cumsum(width)
  grp   <- c(sprintf("$J_0=%s$", J0_values), "EKW $p$-values", paste(cite[["kudrin"]], "$p$-values"))
  head2 <- c("", "$n$", "Articles", "Med.\\ $\\nu$", "$\\Delta_\\nu$", rep(c("$p$", "$\\widehat{B}$"), length(J0_values)),
             "LCM", "Fisher", "Disc.", "Binomial", "CS1", "CS2B", "Min",
             "$T_\\infty$", "$T_H$", "$T_H^{\\mathrm{hybrid}}$", "Min")
  stopifnot(length(head2) == ncol(cells), last[length(last)] == ncol(cells))
  tex_row <- function(x) paste0(paste(x, collapse = " & "), " \\\\")
  
  tex <- c("\\begin{table}[htbp]",
           "\\centering",
           paste0("\\caption{Empirical Application Results: ", spec$label, "}"),
           paste0("\\label{tab:application", suffix, "}"),
           "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{% shrinks the table only if it is wider than the text",
           paste0("\\begin{tabular}{l", strrep("c", ncol(cells) - 1), "}"),
           tex_row(c(rep("", 5), paste0("\\multicolumn{", width, "}{c}{", grp, "}"))),
           paste0("\\cmidrule(lr){", last - width + 1, "-", last, "}", collapse = " "),
           tex_row(head2),
           "\\midrule",
           apply(cells, 1, tex_row),
           "\\bottomrule",
           "\\end{tabular}}",
           "",
           "\\medskip",
           "\\begin{minipage}{\\linewidth}",
           notes,
           "\\end{minipage}",
           "\\end{table}")
  writeLines(tex, file.path(res_dir, out_file))
  saveRDS(tab, file.path(res_dir, paste0("Table_application_numbers", suffix)))
  
  #--- Console version
  txt <- gsub("$^{\\dagger}$", "+", gsub("$\\cdot$", "-", cells[, -1], fixed = TRUE), fixed = TRUE)
  dimnames(txt) <- list(meth, c("n", "Articles", "Med.nu", "Delta_nu", paste0(c("p_J0=", "B_J0="), rep(J0_values, each = 2)),
                                "LCM", "Fisher", "Disc.", "Binomial", "CS1", "CS2B", "EKW_min", "T_inf", "T_H", "T_H_hyb", "Kud_min"))
  cat("\nDe-rounding: ", spec$label, "\n", sep = "")
  print(noquote(txt), right = TRUE)
  if (any(unr)) cat("+ hybrid restriction did not bind, so T_H_hybrid = T_H\n")
  cat("LaTeX table written to", file.path(res_dir, out_file), "\n")
  invisible(list(numbers = tab, latex = tex))
}

table_application <- make_application_table(J0_values = J0_list, n_draws = n_draws)   #e.g. J0_list <- c(30, 20) mirrors the paper's two columns