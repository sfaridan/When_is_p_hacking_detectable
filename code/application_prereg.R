# Application of the projection test, the EKW tests and Kudrin's (2025) tests to the data of Brodeur, Cook,
# Hartley and Heyes (2024), "Do Preregistration and Preanalysis Plans Reduce p-Hacking and Publication Bias?
# Evidence from 15,992 Test Statistics and Suggestions for Improvement," Journal of Political Economy
# Microeconomics 2(3), 527-561 (BCHH below).
#
# The same as the application script for the Brodeur et al. (2020) data, with two changes:
#   1. Every result is broken down by pre-registration status (BCHH's prereg) instead of by method.
#   2. The sample is exactly BCHH's: the drops in their replication do file (do_file_figures.do) and no
#      others. Kudrin's (2025) drop of missing or zero coefficients and the drops of standard errors printed
#      as zero or with kappa above a cap are gone; the rounding diagnostics only report what they would drop.
# Every results file starts with file_tag, so the Brodeur et al. (2020) results are not overwritten.

# 1. Clear all objects from the global environment
rm(list = ls())

# 2. Force R to release the freed memory to your OS
gc(reset = TRUE)

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_all_weighted.R"))

#Load in the data: merged.dta from BCHH's replication package
data_location <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data/dataverse_files"   #the folder holding merged.dta
setwd(data_location)
PAP_raw <- haven::read_dta("merged.dta")

#### BCHH's sample, exactly as in their replication do file (do_file_figures.do)
# Their figures and tests load merged.dta, check it with datasignature (16,390 test statistics and 194
# variables) and run
#     drop if journal=="Experimental Economics"
#     drop if experiment==0
# before splitting by levelsof prereg, so a test statistic without a pre-registration status is in neither
# group, and one without a zstat enters no histogram, density or test. Stata's "drop if experiment==0"
# keeps a missing experiment, and so does %in% below. Their histograms show z < 10 only, but that is a
# display choice, not a drop. prereg = 1 if the study was registered before its trial end date. zstat is
# BCHH's z-statistic: |coefficient/SE|, or converted from a reported t-statistic or p-value.
# Every comparison by prereg uses these drops except figure_a12, which writes the input files for EKW's
# tests: it skips both and instead drops z < 0 (p > 1). bchh_sample = "ekw_input" reproduces that sample.
# (figure_5 and figure_a8 also skip them, but they compare kinds of pre-registration, not prereg itself.)

bchh_sample <- "main"   #"main" = the sample of BCHH's figures and tests; "ekw_input" = the sample of their figure_a12

stopifnot(bchh_sample %in% c("main", "ekw_input"))
if (!identical(dim(PAP_raw), c(16390L, 194L)))
  stop("merged.dta has ", nrow(PAP_raw), " rows and ", ncol(PAP_raw), " variables, but BCHH's do file asserts 16,390 and 194")
need <- c("zstat", "prereg", "journal", "experiment", "title", "coeff_string", "stderror_string")
if (!all(need %in% names(PAP_raw))) stop("merged.dta has no ", paste(setdiff(need, names(PAP_raw)), collapse = ", "))

as_chr <- function(x) as.character(haven::zap_labels(x))   #plain vectors, without Stata's value labels
as_dbl <- function(x) as.numeric(haven::zap_labels(x))

z_all     <- as_dbl(PAP_raw$zstat)
prereg    <- as_dbl(PAP_raw$prereg)
title_all <- as_chr(PAP_raw$title)
is_ee     <- as_chr(PAP_raw$journal) %in% "Experimental Economics"
is_nonexp <- as_dbl(PAP_raw$experiment) %in% 0
in_group  <- prereg %in% c(0, 1) & !is.na(z_all)
main_rows <- !is_ee & !is_nonexp & in_group                  #BCHH's sample
ekw_rows  <- in_group & z_all >= 0                           #the sample of their figure_a12

cat(sprintf(paste0("merged.dta holds %d test statistics. BCHH drop %d in Experimental Economics and %d more from studies that ",
                   "are not experiments; of the rest, %d have no pre-registration status and %d more no zstat.\n"),
            nrow(PAP_raw), sum(is_ee), sum(!is_ee & is_nonexp),
            sum(!is_ee & !is_nonexp & !(prereg %in% c(0, 1))),
            sum(!is_ee & !is_nonexp & prereg %in% c(0, 1) & is.na(z_all))))
sample_counts <- function(rows) {                            #tests, articles and p < .15 (EKW's window) by group
  sel <- list(not_prereg = rows & prereg %in% 0, prereg = rows & prereg %in% 1, total = rows)
  rbind(tests       = sapply(sel, sum),
        articles    = sapply(sel, function(s) length(unique(title_all[s]))),
        p_below_015 = sapply(sel, function(s) sum(2 * (1 - pnorm(z_all[s])) < 0.15)))   #their ptop1 < 0.15
}
cat("\nBCHH's sample (bchh_sample = \"main\"):\n")
print(sample_counts(main_rows))
cat("z < 0 in this sample:", sum(main_rows & z_all < 0), "\n")
cat("\nThe sample of their figure_a12 (bchh_sample = \"ekw_input\"):\n")
print(sample_counts(ekw_rows))
cat("For comparison: the title of BCHH (2024) counts 15,992 test statistics, and the EKW figure of their August 2022 draft",
    "reports 5,616 (not pre-registered) and 2,421 (pre-registered) observations.\n")

PAP_data       <- PAP_raw[if (bchh_sample == "main") main_rows else ekw_rows, ]
PAP_data$group <- ifelse(as_dbl(PAP_data$prereg) == 1, "PreReg", "NoPreReg")
PAP_data$title <- as_chr(PAP_data$title)                     #the article: the cluster in every test (their egen id = group(title))
grp            <- c("NoPreReg", "PreReg")                    #the groups, in the order of every table
grp_label      <- c(NoPreReg = "Not pre-registered", PreReg = "Pre-registered")
file_tag       <- if (bchh_sample == "main") "PAP" else "PAP_ekwinput"   #starts the name of every results file

#### Prep data: de-round the t-statistics (Online Appendix D, Propositions D.1-D.2)
# Four de-rounding options, chosen with deround_option below. For a t-statistic computed from a coefficient
# m and a standard error s printed to multiples of a and b, with s_min = s - b/2 and s_max = s + b/2 the
# smallest and largest standard errors consistent with s:
#   1 "none":       t = BCHH's zstat, exactly as stored (= |m / s| up to float storage)
#   2 "uniform":    t = (m + V) / W,  V ~ U(-a/2, a/2), W ~ U(s_min, s_max)  (the usual de-rounding)
#   3 "smin":       t = m / s_min                                              (Appendix D.1)
#   4 "smin_noise": t = (m + V) / s_min
# BCHH record each coefficient and standard error as printed (coeff_string, stderror_string), and zstat is
# |coefficient/SE| for tests reported that way and converted from the reported t-statistic or p-value
# otherwise. Under options 2 to 4 a t-statistic is de-rounded as m/s when both strings are plain numbers, the
# SE is positive and |m/s| reproduces zstat. Every other t-statistic is BCHH's zstat, used as stored under
# option 3 and with noise uniform on +-1/2 of its last digit under options 2 and 4. Option 1 uses zstat as
# stored for every t-statistic, with no rounding at all: a p-value converted from a reported round number
# such as 0.05 sits exactly on a threshold or bin edge of the EKW tests, so moving it by even 1e-7 changes
# their results. Options 2 and 4 are run on n_draws independent draws, and p-values and breakdown
# statistics are averaged over them. Not covered: t-statistics converted from rounded p-values.
#
# Nothing is dropped here, so the sample is BCHH's under every option. The noise is drawn once, with fixed
# seeds, so options 2 and 4 use identical numerator noise and differ only in the denominator.
#
# Rounding units:
#   b, per entry: the place value of the last digit of the SE as BCHH recorded it. A string can keep
#     trailing zeros (the diagnostics below check that it does), so b is exact wherever the coders typed
#     them; a dropped trailing zero only makes b too large, which keeps s_min a lower bound.
#   a, per table: the finest place value among the table's printed coefficients, but no finer than b, since
#     a coefficient and its standard error are almost always printed to the same decimals. Both are fixed
#     by the table's format and the standard error, not by the estimate. If merged.dta has no table
#     variable, the article stands in for the table. The width of V need not match the true unit:
#     Proposition D.1 allows any independent noise.

deround_option <- 2     #1 = none, 2 = uniform (standard), 3 = smallest SE, 4 = smallest SE + noise in the numerator
se_rounded_up  <- FALSE #TRUE if papers may have rounded standard errors up: then s_min = sd - b and s_max = sd
n_draws        <- 100   #independent de-roundings for options 2 and 4; try 2 first to check the run time

stopifnot(deround_option %in% 1:4)
deround_spec <- list(key   = c("none", "uniform", "smin", "smin_noise")[deround_option],
                     label = c("No de-rounding", "Uniform de-rounding", "Smallest consistent standard error",
                               "Smallest consistent standard error with uniform noise")[deround_option],
                     noise = c(FALSE, TRUE, FALSE, TRUE)[deround_option])
add_noise <- deround_spec$noise
suffix    <- paste0("_", deround_spec$key)                   #labels every results file and the table
cat("\nDe-rounding option ", deround_option, ": ", deround_spec$label, "\n", sep = "")

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

parse_printed <- function(x) {         #a number as BCHH recorded it: value, decimals and digits; NA unless a plain decimal
  s <- as_chr(x)
  s[is.na(s)] <- ""
  s <- gsub("[\\s*()\\[\\]]", "", s, perl = TRUE)             #spaces, stars, parentheses and brackets
  s <- gsub("\u2212", "-", s, fixed = TRUE)                   #Unicode minus, as text copied from a PDF has it
  s <- gsub("\u2013", "-", s, fixed = TRUE)                   #en dash used as a minus
  ok  <- grepl("^[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)$", s)
  dec <- ifelse(grepl(".", s, fixed = TRUE), nchar(sub("^[^.]*\\.", "", s)), 0L)
  list(str   = s,
       value = ifelse(ok, suppressWarnings(as.numeric(s)), NA_real_),
       dec   = ifelse(ok, dec, NA_integer_),
       sig   = ifelse(ok, suppressWarnings(as.numeric(gsub("[^0-9]", "", s))), NA_real_))   #the digits read as an integer
}

n_obs <- nrow(PAP_data)
cf    <- parse_printed(PAP_data$coeff_string)                #coefficient as printed
se    <- parse_printed(PAP_data$stderror_string)             #standard error as printed
z_p   <- as_dbl(PAP_data$zstat)                              #BCHH's z-statistic, never rounded (see option 1 above)
dec_z <- n_dec(z_p)
a_t   <- 10^(-dec_z)                                         #its rounding unit, if it is a printed t-statistic (noise width only)
mu_p  <- cf$value
sd_p  <- se$value

tab_var <- intersect(c("table", "table_number", "table_no", "tablenumber"), names(PAP_data))[1]
tab_id  <- if (is.na(tab_var)) PAP_data$title else paste(PAP_data$title, as_chr(PAP_data[[tab_var]]), sep = " | ")
b_sd    <- 10^(-se$dec)
a_tab   <- 10^(-ave(cf$dec, tab_id, FUN = function(d) if (all(is.na(d))) NA_integer_ else max(d, na.rm = TRUE)))
a_mu    <- pmax(a_tab, b_sd)                                 #NA unless both coefficient and SE are printed
s_min   <- sd_p - if (se_rounded_up) b_sd else b_sd / 2      #smallest SE consistent with sd
s_max   <- sd_p + if (se_rounded_up) 0 else b_sd / 2         #largest SE consistent with sd

has_cs    <- !is.na(mu_p) & !is.na(sd_p)                     #coefficient and SE both printed as plain numbers
agrees    <- has_cs & sd_p > 0 &                             #and their ratio is zstat, up to float storage or zstat's rounding
  abs(abs(mu_p / sd_p) - abs(z_p)) <= pmax(1e-3 * pmax(1, abs(z_p)), 0.5001 * a_t)
use_ratio <- agrees & s_min > 0                              #de-rounded as coefficient/SE; every other t-statistic is zstat

t_why <- ifelse(use_ratio, "coef/SE",
                ifelse(!(nzchar(cf$str) & nzchar(se$str)), "no coefficient and SE",
                       ifelse(!has_cs, "not plain numbers",
                              ifelse(sd_p <= 0, "SE zero or negative",
                                     ifelse(!agrees, "coef/SE is not zstat", "s_min is zero")))))
cat("\nHow each t-statistic is formed (coef/SE is de-rounded as a ratio; every other one is BCHH's zstat):\n")
print(table(grp_label[PAP_data$group], t_why))
if ("report" %in% names(PAP_data)) {
  cat("\nBy BCHH's reporting code (report == \"s\": coefficient and SE):\n")
  print(table(report = as_chr(PAP_data$report), t_why, useNA = "ifany"))
}
bad <- nzchar(cf$str) & nzchar(se$str) & !has_cs
if (any(bad)) {
  cat("\nCoefficients or SEs that are not plain numbers (their t-statistics stay BCHH's zstat), for example:\n")
  print(head(unique(paste(as_chr(PAP_data$coeff_string)[bad], as_chr(PAP_data$stderror_string)[bad], sep = "  |  ")), 10))
}
cat("\nDecimals printed in the standard errors (many at 6+ would suggest computed rather than printed SEs):\n")
print(table(se$dec[use_ratio]))
cat(sprintf("Standard errors ending in a zero after the decimal point, as in 0.010: %.1f percent (roughly 10 if the strings keep trailing zeros, none if they lost them).\n",
            100 * mean(grepl("\\.[0-9]*0$", se$str[use_ratio]))))
cat("Coefficient rounding unit a set per", if (is.na(tab_var)) "article (merged.dta has no table variable)" else paste0("table (variable ", tab_var, ")"), "\n")

#### Rounding flags: what the dropping rules would drop in BCHH's sample (nothing is dropped)
# The script for the Brodeur et al. (2020) data dropped t-statistics whose SE is printed as zero or whose
# kappa = a/s_min exceeds 1 ("ours" below), and compared that with Kranz and Putz's (2022) rule, which omits
# every test whose standard error has a significand below 37: the SE's significant digits read as an
# integer, e.g. 0.012 -> 12. With the coefficient and SE printed to the same decimals, this keeps rounding
# from moving any z-statistic across 1.96 unless it lies within 0.04 of it (their Lemma 1). BCHH record the
# SE as printed, so its significand is read from its own digits (0.010 -> 10) rather than from KP's fallback
# for data that lost trailing zeros. Both rules concern tests with a printed coefficient and SE; t-statistics
# reported otherwise pass both.

kp_threshold <- 37
kappa_flag   <- 1
flag_kp      <- has_cs & se$sig < kp_threshold                                  #includes SEs printed as zero
flag_ours    <- (has_cs & sd_p <= 0) | (use_ratio & a_mu / s_min > kappa_flag)
z_is_2       <- has_cs & sd_p > 0 & abs(abs(mu_p / sd_p) - 2) < 1e-9             #KP's heap at exactly z = 2

flag_summary <- function(g) {                                                   #g: logical row selector
  c(t_stats     = sum(g),
    coef_se     = sum((has_cs & g)),
    ours        = sum(flag_ours[g]),            kp     = sum(flag_kp[g]),
    ours_pct    = 100 * mean(flag_ours[g]),     kp_pct = 100 * mean(flag_kp[g]),
    both        = sum((flag_ours & flag_kp)[g]),
    only_ours   = sum((flag_ours & !flag_kp)[g]),
    only_kp     = sum((!flag_ours & flag_kp)[g]),
    z_eq_2      = sum(z_is_2[g]),
    ours_z2_pct = 100 * mean(flag_ours[g & z_is_2]),
    kp_z2_pct   = 100 * mean(flag_kp[g & z_is_2]))
}
rows_fl         <- c(setNames(lapply(grp, function(k) PAP_data$group %in% k), grp_label[grp]),
                     list(Pooled = rep(TRUE, n_obs)))
flag_comparison <- t(sapply(rows_fl, flag_summary))

cat("\nt-statistics the dropping rules would drop, all kept here: SE printed as zero or kappa > ", kappa_flag,
    " (ours) and SE significand < ", kp_threshold, " (Kranz and Putz):\n", sep = "")
print(round(flag_comparison, 1))
cat(sprintf("Kept here: our rule would drop %.1f and %.1f percent of the t-statistics (not pre-registered, pre-registered; %.1f pooled), Kranz and Putz's %.1f and %.1f percent (%.1f pooled).\n",
            flag_comparison[1, "ours_pct"], flag_comparison[2, "ours_pct"], flag_comparison["Pooled", "ours_pct"],
            flag_comparison[1, "kp_pct"],   flag_comparison[2, "kp_pct"],   flag_comparison["Pooled", "kp_pct"]))
cat(sprintf("Of the %d t-statistics at exactly z = 2, our rule would drop %.1f percent and Kranz and Putz's %.1f percent.\n",
            flag_comparison["Pooled", "z_eq_2"], flag_comparison["Pooled", "ours_z2_pct"], flag_comparison["Pooled", "kp_z2_pct"]))
saveRDS(flag_comparison, file = file.path(root, "results", paste0(file_tag, "_Flag_comparison_KP")))

#Keep the printed values each option needs (nothing is dropped)
PAP_data$t_source   <- ifelse(use_ratio, "coef/SE", "BCHH zstat")
PAP_data$dr_mu      <- mu_p                                       #printed coefficient
PAP_data$dr_sd      <- sd_p                                       #printed standard error
PAP_data$dr_smin    <- s_min                                      #smallest SE consistent with it
PAP_data$dr_smax    <- s_max                                      #largest SE consistent with it
PAP_data$dr_t       <- z_p                                        #BCHH's zstat, exactly as stored
PAP_data$dr_a       <- ifelse(use_ratio, a_mu, a_t)               #rounding unit of the numerator
PAP_data$kappa      <- ifelse(use_ratio, a_mu / s_min, a_t)       #rounding unit in t-units
PAP_data$lambda_max <- ifelse(use_ratio, s_max / s_min, NA)       #largest inflation from dividing by s_min

#The de-rounding noise, drawn once with fixed seeds, so every option and every run uses the same draws
if (!add_noise) n_draws <- 1                                 #without noise every draw would be the same
set.seed(1)
U_num <- matrix(runif(nrow(PAP_data) * n_draws, -0.5, 0.5), ncol = n_draws)   #numerator noise, one column per draw
set.seed(2)
U_den <- matrix(runif(nrow(PAP_data) * n_draws, -0.5, 0.5), ncol = n_draws)   #denominator noise (option 2 only)

deround_draw <- function(r) {                                #t-statistics for draw r under the chosen option
  d     <- PAP_data
  if (deround_spec$key == "none") return(d$dr_t)             #option 1: BCHH's zstat for every t-statistic
  ratio <- d$t_source %in% "coef/SE"
  u_num <- if (add_noise) U_num[, r] else 0
  u_den <- if (add_noise) U_den[, r] else 0
  num   <- d$dr_mu + d$dr_a * u_num                          #coefficient, plus noise under options 2 and 4
  den   <- switch(deround_spec$key,
                  uniform    = (d$dr_smin + d$dr_smax) / 2 + (d$dr_smax - d$dr_smin) * u_den,  #uniform on [s_min, s_max]
                  smin       = d$dr_smin,
                  smin_noise = d$dr_smin)
  ifelse(ratio, num / den, d$dr_t + d$dr_a * u_num)          #zstat: as recorded, plus noise under options 2 and 4
}

if (any(!is.finite(deround_draw(1)))) {
  stop("Tscores not clean")
}

if (deround_spec$key %in% c("smin", "smin_noise")) {         #these diagnostics concern dividing by s_min
  cat("\nLargest possible inflation of t, s_max/s_min (median, 99th percentile, max):\n")
  print(round(sapply(split(PAP_data$lambda_max, grp_label[PAP_data$group]), quantile, probs = c(.5, .99, 1), na.rm = TRUE), 3))
  kappa6 <- tapply(PAP_data$kappa^6, grp_label[PAP_data$group], mean, na.rm = TRUE)
  cat("\nRounding allowance (Appendix D.1, step 3), to compare with B-hat. kappa is not capped here, so a few coarsely",
      "printed coefficients can dominate the mean (the flags above count kappa > 1):\n")
  print(signif(rbind(mean_kappa6 = kappa6, bound_J0_20 = 0.13 * kappa6, bound_J0_30 = 0.31 * kappa6), 3))
} else {
  cat("\nNo rounding allowance is reported: Appendix D's bound covers only options 3 and 4.\n")
}

#Averaging over draws
as_num1         <- function(v) as.numeric(v)[1]              #test results may come back as lists or list-columns
draw_range      <- function(x) diff(range(x))                #max - min across draws
summarize_draws <- function(d, cols, f = mean) {             #one row per group, in the order of grp;
  out <- d[match(grp, d$group), setdiff(names(d), "draw")]   #cols summarized over draws by f
  for (k in cols) out[[k]] <- as.numeric(tapply(d[[k]], d$group, f)[grp])
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
  t_start    <- Sys.time()
  PAP_data$t <- deround_draw(r)
  set.seed(1000 * r + 1) #bootstrap seeds differ across draws, so bootstrap noise averages out too
  for (J0 in J0_list) {
    for (g in grp) {
      res <- run_test_weight_fast(PAP_data[PAP_data$group==g,],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
      proj_all <- rbind(proj_all, data.frame(draw = r, group = g, n = as_num1(res$n), articles = as_num1(res$articles),
                                             pval = as_num1(res$pval), breakdown = as_num1(res$breakdown),
                                             numcoeffs = as_num1(res$numcoeffs), J0 = J0))
    }
  }
  cat(sprintf("Projection test (%s): de-rounding %d of %d done in %.1f minutes\n", deround_spec$key, r, n_draws,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

#Save the averages over draws (the files make_application_table() reads) and every draw, labeled by option
setwd(paste0(root,"/results"))
saveRDS(proj_all, file = paste0(file_tag, "_weights_results_draws", suffix))
for (J0 in J0_list) {
  name <- paste0(file_tag, "_weights_results_", J0, "_J0", suffix)
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

#p-values built as in BCHH's EKW input (figure_a12 of their do file): gen ptop1 = 2*(1 - normal(zstat)), a float.
#This differs from 2*pnorm(-abs(t)) where it matters to EKW's tests: beyond |t| of about 8.3 it is exactly 0, outside
#the tests' window [pmin, 0.15], where 2*pnorm(-abs(t)) is tiny but positive and lands in the first bin; and the float
#decides on which side of a bin edge or threshold a p-value converted from a round number like 0.05 falls.
p_bchh <- function(t) {
  p <- 2 * (1 - pnorm(abs(t)))
  readBin(writeBin(p, raw(), size = 4), "double", n = length(p), size = 4)   #stored as a Stata float, gen's default type
}

#### Optional check against BCHH's own EKW input files
# Their do file writes figure_a12_input_to_elliot_prereg_0.csv and _1.csv (columns ptop1, ptop2, id), in the order of
# merged.dta. Set bchh_ekw_dir to the folder holding them and the script checks, row by row, that it rebuilds their
# p-values from merged.dta, then runs the EKW tests on their files. If those tests reproduce their published numbers,
# any remaining gap comes from the data; if not, from the test code or its settings. The files pass through export
# delimited, which writes floats to about 8 digits, so differences of order 1e-9 are formatting, not data.
bchh_ekw_dir <- NULL   #e.g. "C:/Users/stefa/Downloads/BCHH replication"; NULL skips the check

if (!is.null(bchh_ekw_dir)) {
  for (lv in 0:1) {
    f <- file.path(bchh_ekw_dir, paste0("figure_a12_input_to_elliot_prereg_", lv, ".csv"))
    if (!file.exists(f)) { cat("\nNo file", f, "\n"); next }
    theirs <- read.csv(f)
    ours   <- p_bchh(z_all[ekw_rows & prereg %in% lv])        #their EKW sample, in merged.dta order
    cat(sprintf("\nprereg = %d: their file has %d p-values, merged.dta gives %d", lv, nrow(theirs), length(ours)))
    if (nrow(theirs) == length(ours)) {
      gap <- abs(theirs$ptop1 - ours)
      cat(sprintf("; %d differ by more than 1e-12 (largest gap %.2g, %d by more than 1e-6)",
                  sum(gap > 1e-12), max(gap), sum(gap > 1e-6)))
    }
    cat("\nEKW tests on their file:\n")
    set.seed(1000 * 1 + 2)                                     #the seed of draw 1 below, so the two runs are comparable
    print(round(vapply(run_ewk(theirs$ptop1, ids = theirs$id), as_num1, numeric(1)), 3))
  }
}

ewk_draws <- vector("list", n_draws)
for (r in seq_len(n_draws)) {
  PAP_data$t <- deround_draw(r)
  set.seed(1000 * r + 2)
  ewk_draws[[r]] <- sapply(setNames(grp, paste0("ewk_", tolower(grp))), function(g) {   #tests x groups
    s <- PAP_data$group == g
    vapply(run_ewk(p_bchh(PAP_data$t[s]), ids = as.integer(factor(PAP_data$title[s]))), as_num1, numeric(1))
  })
}
ewk_all       <- simplify2array(ewk_draws)                   #tests x groups x draws
empirical_ewk <- apply(ewk_all, c(1, 2), mean)
ewk_tests     <- setdiff(rownames(empirical_ewk), "min")
empirical_ewk["min", ] <- apply(empirical_ewk[ewk_tests, , drop = FALSE], 2, min)   #smallest averaged p-value
setwd(paste0(root,"/results"))
saveRDS(ewk_all, file = paste0(file_tag, "_EWK_results_draws", suffix))
saveRDS(empirical_ewk, file = paste0(file_tag, "_EWK_results", suffix))
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
  kt <- kudrin_tests(data$t, cluster = NULL, boots = boots) #no clustering, for comparison
  data.frame(n = length(data$t), articles = length(unique(data$title)), as.list(kt$pval),
             min = min(kt$pval), restricted = kt$restricted) #restricted = FALSE means T_H_hybrid = T_H
}

kud_all <- NULL
for (r in seq_len(n_draws)) {
  t_start    <- Sys.time()
  PAP_data$t <- deround_draw(r)
  set.seed(1000 * r + 3) #T_inf is bootstrapped
  kud_all <- rbind(kud_all, cbind(draw = r, group = grp,
                                  do.call(rbind, lapply(grp, function(g) run_kudrin(PAP_data[PAP_data$group == g, ])))))
  cat(sprintf("Kudrin tests (%s): de-rounding %d of %d done in %.1f minutes\n", deround_spec$key, r, n_draws,
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
kud_tests        <- setdiff(names(kud_all), c("draw", "group", "n", "articles", "min", "restricted"))
empirical_kudrin <- summarize_draws(kud_all, kud_tests)
empirical_kudrin$min        <- apply(empirical_kudrin[kud_tests], 1, min)   #smallest averaged p-value
empirical_kudrin$restricted <- as.logical(tapply(kud_all$restricted, kud_all$group, any)[grp])   #bound in any draw
setwd(paste0(root,"/results"))
saveRDS(kud_all, file = paste0(file_tag, "_Kudrin_results_draws", suffix))
saveRDS(empirical_kudrin, file = paste0(file_tag, "_Kudrin_results", suffix))
print(readRDS(paste0(file_tag, "_Kudrin_results", suffix)))

#How much do the results move across de-roundings? (max - min over draws)
if (n_draws > 1) {
  cat("\nSpread across de-roundings (max - min):\n")
  for (J0 in J0_list) {
    cat("Projection test, J0 =", J0, "\n")
    print(summarize_draws(proj_all[proj_all$J0 == J0, ], c("pval", "breakdown"), draw_range)[, c("group", "pval", "breakdown")])
  }
  cat("EKW tests\n")
  print(round(apply(ewk_all, c(1, 2), draw_range), 3))
  cat("Kudrin (2025) tests\n")
  print(summarize_draws(kud_all, kud_tests, draw_range)[, c("group", kud_tests)])
}

#### Table: projection test, EKW and Kudrin (2025) tests by pre-registration status (layout of Table 3 in the paper)
# Reads the results files saved above for one de-rounding option and one sample (so after a full run you can
# tweak and re-run this block alone), prints the table and writes results/table_application_<file_tag><suffix>.tex.
# LaTeX needs booktabs, graphicx (shrink-to-fit wrapper) and hyperref or url (for \url in the notes).

make_application_table <- function(J0_values = c(20),            # one (p, B-hat) group per J0, in this order; each must be in J0_list above
                                   n_draws   = 1,                # number of de-roundings the results files average over
                                   spec      = deround_spec,     # the de-rounding option the results files were run with
                                   suffix    = paste0("_", spec$key),
                                   tag       = file_tag,         # the sample the results files were run on
                                   sample    = bchh_sample,
                                   res_dir   = paste0(root, "/results"),
                                   out_file  = paste0("table_application_", tag, suffix, ".tex")) {
  grp <- c("NoPreReg", "PreReg")
  lab <- c(NoPreReg = "Not pre-registered", PreReg = "Pre-registered")
  
  # Effective sample sizes: not collected for this sample, so both columns show a centered dot.
  # Type them in here (as for the hand-collected samples of Brodeur et al. (2020)) to fill them.
  med_nu   <- c(NoPreReg = NA_real_, PreReg = NA_real_)
  Delta_nu <- c(NoPreReg = NA_real_, PreReg = NA_real_)
  
  # Citations used in the header and notes; swap in "\\citet{key}" to get hyperlinked citations
  cite <- c(bchh    = "Brodeur et al. (2024)", ekw  = "Elliott et al. (2022)",  kudrin    = "Kudrin (2025)",
            coxshi  = "Cox and Shi (2023)",    osqp = "Stellato et al. (2024)", brodeur16 = "Brodeur et al. (2016)")
  
  num  <- function(x) vapply(x, function(v) as.numeric(v)[1], numeric(1), USE.NAMES = FALSE) #results files hold list-columns
  pick <- function(d, what) {                                    #rows of d in the order of grp
    i <- match(grp, d$group)
    if (anyNA(i)) stop(what, " has no results for ", paste(grp[is.na(i)], collapse = ", "))
    d[i, ]
  }
  
  #--- Read the saved results for this option and sample
  res_file <- function(stem) {
    f <- file.path(res_dir, paste0(tag, "_", stem, suffix))
    if (!file.exists(f)) stop("No results file ", basename(f), ": run the tests with this de-rounding option and sample first")
    f
  }
  proj <- lapply(J0_values, function(j0) {
    f <- res_file(paste0("weights_results_", j0, "_J0"))
    r <- pick(readRDS(f), basename(f))
    if (any(num(r$numcoeffs) != coeffs)) warning(basename(f), " was run with numcoeffs != ", coeffs, " (stale file?)")
    data.frame(n = num(r$n), articles = num(r$articles), p = num(r$pval), B = num(r$breakdown))
  })
  ewk <- readRDS(res_file("EWK_results"))                        #matrix: tests x groups
  ewk <- matrix(num(ewk), nrow(ewk), dimnames = dimnames(ewk))[, paste0("ewk_", tolower(grp))]
  kud <- pick(readRDS(res_file("Kudrin_results")), "Kudrin_results")
  
  n_all <- cbind(sapply(proj, function(d) d$n),        num(kud$n))   #same samples in every file?
  a_all <- cbind(sapply(proj, function(d) d$articles), num(kud$articles))
  if (any(n_all != n_all[, 1]) || any(a_all != a_all[, 1]))
    warning("n or article counts differ across results files: were they all run on the same PAP_data?")
  
  #--- Numbers
  tab <- data.frame(n = proj[[1]]$n, articles = proj[[1]]$articles,
                    med_nu = med_nu[grp], Delta_nu = Delta_nu[grp], row.names = grp)
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
  cells <- cbind(unname(lab[grp]), sprintf("%.0f", tab$n), sprintf("%.0f", tab$articles),
                 ifelse(is.na(tab$med_nu), "$\\cdot$", as.character(tab$med_nu)), fmt(tab$Delta_nu, 4, lead0 = TRUE))
  for (j0 in J0_values) cells <- cbind(cells, fmt(tab[[paste0("p_J0_", j0)]], 2), fmt(tab[[paste0("B_J0_", j0)]], 4))
  for (k in c(names(ekw_rows), names(kud_cols))) cells <- cbind(cells, fmt(tab[[k]], 2))
  unr <- !tab$hybrid_restricted                                  #dagger: hybrid restriction did not bind
  if (any(unr)) cells[unr, ncol(cells) - 1] <- paste0(cells[unr, ncol(cells) - 1], "$^{\\dagger}$")
  
  #--- Notes (one sentence per line in the .tex; edit freely)
  sample_note <- c(
    main      = paste0("Results for the ", cite[["bchh"]], " data by pre-registration status, on their sample: following their replication files, we drop articles in ",
                       "\\textit{Experimental Economics} and studies that are not experiments, and keep every test statistic with a recorded pre-registration status and $z$-statistic."),
    ekw_input = paste0("Results for the ", cite[["bchh"]], " data by pre-registration status, on the sample their replication files pass to the tests of ", cite[["ekw"]],
                       ": every test statistic with a recorded pre-registration status and a nonnegative $z$-statistic."))[[sample]]
  deround_notes <- c(
    none       = paste0("De-rounding: none. Each $t$-statistic is the $z$-statistic of ", cite[["bchh"]],
                        ": the ratio of the reported coefficient and standard error, or converted from a reported $t$-statistic or $p$-value."),
    uniform    = paste0("De-rounding: uniform, following ", cite[["brodeur16"]], ". Where the coefficient and standard error are reported, each is replaced by an independent uniform draw from its rounding interval, ",
                        "and any other $t$-statistic gets uniform noise at its last digit."),
    smin       = paste0("De-rounding: smallest consistent standard error (Appendix D.1). Where the coefficient and standard error are reported, the coefficient is divided by the smallest standard error consistent with the reported one, ",
                        "and any other $t$-statistic is used as reported."),
    smin_noise = paste0("De-rounding: smallest consistent standard error with uniform noise (Appendix D.1). Where the coefficient and standard error are reported, the coefficient, plus uniform noise at its last digit, ",
                        "is divided by the smallest standard error consistent with the reported one, and any other $t$-statistic gets uniform noise at its last digit."))
  nu_missing <- is.na(tab$med_nu) | is.na(tab$Delta_nu)
  nu_notes <- if (all(nu_missing)) {
    "We did not collect effective sample sizes for this sample, so Med.\\ $\\nu$ and $\\Delta_\\nu$ are not reported."
  } else c(
    "Med.\\ $\\nu$ is the median effective sample size as measured by the median number of clusters in a random sample of 100 articles for each group.",
    "In ambiguous cases, the smallest cluster count was always used.",
    paste0("$\\Delta_\\nu=\\frac{1}{100}\\sum_{k=1}^{100}\\sup_h\\|\\varphi(t-h)-t_{\\nu_k,h}(t)\\|_T$ is the average difference between the PDFs of the standard normal ",
           "and Student-$t$ over the sample of 100 effective sample sizes collected (computed numerically)."),
    if (any(nu_missing)) "A centered dot marks a group for which they were not computed.")
  notes <- c(
    sample_note,
    "A study counts as pre-registered if it was registered before its trial end date.",
    deround_notes[[spec$key]],
    "No test statistic is dropped for rounding, so the sample is the same for every de-rounding choice.",
    if (n_draws > 1) sprintf("All $p$-values and breakdown statistics are averages over %d independent de-roundings.", n_draws),
    paste0("Columns headed $J_0$ report $p$-values and breakdown statistics $\\widehat{B}$ of the projection test with $J=", J, "$",
           if (length(J0_values) > 1) "; several values of $J_0$ are included to show low sensitivity to the smoothing parameter", "."),
    nu_notes,
    sprintf("We form $\\mathcal{X}$ using a grid of %d elements evenly spaced between $\\pm %s$ plus the scalar 9999.", numgrid, L),
    sprintf("Discretization incurs error $\\mathcal{E}(L,\\delta)=%.5f$, which is already accounted for in $\\widehat{B}$ and the $p$-value.", epsilon_U_pre),
    if (all(nu_missing)) "The discretization was chosen to be made of round numbers and result in a small error."
    else "The discretization was chosen to be made of round numbers and result in error relatively small compared to $\\Delta_\\nu$.",
    paste0("EKW $p$-values are for the tests of ", cite[["ekw"]], ": the least concave majorant test (LCM) and Fisher's test (Fisher) on $p$-values in $(0,0.15]$, ",
           "the density discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$ (Binomial), and the ", cite[["coxshi"]],
           " histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) and for 2-monotonicity with bounds (CS2B)."),
    paste0("As in the input files of ", cite[["bchh"]], ", these $p$-values are $2(1-\\Phi(|t|))$ in single precision."),
    paste0(cite[["kudrin"]], " $p$-values are for his KS-type test ($T_\\infty$), histogram-based test ($T_H$) and hybrid histogram-based test ($T_H^{\\mathrm{hybrid}}$), ",
           "using the sinc kernel with $\\kappa=3$ and, for $T_H$ and $T_H^{\\mathrm{hybrid}}$, the points $0.01,0.02,\\ldots,0.15$."),
    "All three cluster by article: $T_H$ and $T_H^{\\mathrm{hybrid}}$ use cluster-robust covariance matrices and $T_\\infty$ uses an article-level bootstrap.",
    if (any(unr)) "$^\\dagger$The deconvolution estimate already satisfies the hybrid test's shape restrictions, so $T_H^{\\mathrm{hybrid}}=T_H$.",
    if (n_draws > 1) "Min is the smallest of the averaged $p$-values in each group of tests." else "Min is the smallest $p$-value in each group of tests.",
    sprintf("The projection test and $T_\\infty$ use %d bootstrap repetitions.", boots),
    paste0("Computation used osqp package in R by ", cite[["osqp"]], ". Data: ", cite[["bchh"]], "."),
    "Transparency R package available at: \\url{https://github.com/sfaridan/When_is_p_hacking_detectable}")
  
  #--- LaTeX
  width   <- c(rep(2, length(J0_values)), 7, 4)                  #column groups after the 5 leading columns
  last    <- 5 + cumsum(width)
  col_grp <- c(sprintf("$J_0=%s$", J0_values), "EKW $p$-values", paste(cite[["kudrin"]], "$p$-values"))
  head2   <- c("", "$n$", "Articles", "Med.\\ $\\nu$", "$\\Delta_\\nu$", rep(c("$p$", "$\\widehat{B}$"), length(J0_values)),
               "LCM", "Fisher", "Disc.", "Binomial", "CS1", "CS2B", "Min",
               "$T_\\infty$", "$T_H$", "$T_H^{\\mathrm{hybrid}}$", "Min")
  stopifnot(length(head2) == ncol(cells), last[length(last)] == ncol(cells))
  tex_row <- function(x) paste0(paste(x, collapse = " & "), " \\\\")
  
  tex <- c("\\begin{table}[htbp]",
           "\\centering",
           paste0("\\caption{Empirical Application Results by Pre-Registration Status: ", spec$label, "}"),
           paste0("\\label{tab:application_", tag, suffix, "}"),
           "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{% shrinks the table only if it is wider than the text",
           paste0("\\begin{tabular}{l", strrep("c", ncol(cells) - 1), "}"),
           tex_row(c(rep("", 5), paste0("\\multicolumn{", width, "}{c}{", col_grp, "}"))),
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
  saveRDS(tab, file.path(res_dir, paste0(tag, "_Table_application_numbers", suffix)))
  
  #--- Console version
  txt <- gsub("$^{\\dagger}$", "+", gsub("$\\cdot$", "-", cells[, -1], fixed = TRUE), fixed = TRUE)
  dimnames(txt) <- list(unname(lab[grp]), c("n", "Articles", "Med.nu", "Delta_nu", paste0(c("p_J0=", "B_J0="), rep(J0_values, each = 2)),
                                            "LCM", "Fisher", "Disc.", "Binomial", "CS1", "CS2B", "EKW_min", "T_inf", "T_H", "T_H_hyb", "Kud_min"))
  cat("\nBCHH data by pre-registration status (", sample, " sample). De-rounding: ", spec$label, "\n", sep = "")
  print(noquote(txt), right = TRUE)
  if (any(unr)) cat("+ hybrid restriction did not bind, so T_H_hybrid = T_H\n")
  cat("LaTeX table written to", file.path(res_dir, out_file), "\n")
  invisible(list(numbers = tab, latex = tex))
}

table_application <- make_application_table(J0_values = J0_list, n_draws = n_draws)   #e.g. J0_list <- c(30, 20) mirrors the paper's two columns

## effective sample size

hf <- function(x, J) {                      # psi_j(x) * phi(x), j = 0..J
  G <- matrix(0, length(x), J + 1)
  G[, 1] <- dnorm(x); G[, 2] <- x * G[, 1]
  for (j in 1:(J - 1)) G[, j + 2] <- (x * G[, j + 1] - sqrt(j) * G[, j]) / sqrt(j + 1)
  G
}
xi <- function(t, J = 40) 0.5 * (hf(t - 1.96, J) + hf(-t - 1.96, J))
w  <- pmin(1, 2^(20 - (0:40)))              # J0 = 20 weights

deff <- function(d) {
  X <- sweep(xi(d$t), 2, w, `*`)
  R <- sweep(X, 2, colMeans(X))
  S <- rowsum(R, d$title)                   # paper-level sums
  c(trace = sum(S^2) / sum(R^2), j0 = sum(S[, 1]^2) / sum(R[, 1]^2))
}
print(sapply(split(PAP_data, grp_label[PAP_data$group]), deff)) # by pre-registration status