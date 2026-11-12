#### =====================================================================================================
#### application_MM_analysis.R: every computation for the application to Brodeur, Cook and Heyes (2020)
#### Under each de-rounding option it runs the projection test (J0 = 20 and 30, and without the 1.96 shift),
#### the EKW tests, Kudrin's (2025) tests and Stouffer's p-curve test; then it computes the design effects,
#### the data behind the two figures and the run times. It saves all of them, and every setting of the run,
#### in <root>/results. It writes no table or figure: application_MM_output.R reads these files, the settings
#### included, and makes all of them, so their notes always describe the run that produced the numbers.
#### Every file is tagged with run_tag (see drop_zero_mu, require_se, drop_duplicates and run_label): runs on different samples never
#### overwrite each other. Settings_application<run_tag> records the run. Also written, as application_cochrane.R does:
####   Sample_report_application<run_tag>.txt  every step from the file to the final sample, by design, reconciled with
####                                           Kudrin's (2025) Table 2 and Kranz and Putz's (2022) rule
####   Sample_rows_application<run_tag>.csv    every row of the file: kept or not, and why (merge on row)
####   Stability_application<run_tag>.csv      each test's p-value (and B-hat) across de-rounding draws: mean, standard
####                                           deviation, range, and the share of draws that reject
#### =====================================================================================================
rm(list = ls())
gc(reset = TRUE)

#### 0. Settings: every choice of this run, saved with its results -------------------------------------
root          <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
data_location <- "C:/Users/sfaridani6/Documents/Data/120246-V1"
res_dir       <- file.path(root, "results")

# Sample
drop_zero_mu   <- FALSE  #Kudrin's (2025) drop of coefficients printed as zero; FALSE keeps them (missing ones are
#dropped either way). Proposition 4 covers selection on how results are printed, not on
#the estimates, so keeping them is the sample the validity argument covers
kappa_cap      <- 0.5    #drop t-statistics with kappa = a/s_min above this (the same under every option)
require_se     <- TRUE   #keep a test only if its coefficient and its standard error are both printed, so every t is
#coefficient / SE and de-rounded as in Section 5.2. FALSE (a robustness check; the run tag
#then ends in "_keept") also keeps t-statistics printed without a standard error, under the
#same cap with kappa = a_t; t's the file computed from printed p-values are dropped either way
se_rounded_up  <- FALSE  #TRUE if papers may have rounded standard errors up: then s_min = sd - b and s_max = sd
drop_duplicates <- TRUE #TRUE: keep only the first copy of a test the file lists more than once (same design, article,
#printed coefficient, standard error and t), as application_cochrane.R keeps the first copy of each
#trial's outcome data; the run tag then gains "_dedup". FALSE keeps them, as Kudrin (2025) does
run_label      <- "_drop_se_mi"     #anything else that should tell this run's files apart, e.g. "_check"

# De-rounding
options_to_run <- 1:4    #1 = none, 2 = uniform (standard), 3 = smallest SE, 4 = smallest SE + noise in the numerator
n_draws_noise  <- 200    #independent de-roundings for options 2 and 4; try 2 first to check the run time

# Tests
run_tests      <- TRUE   #FALSE: keep the test results an earlier run with these settings saved, and redo the rest
cv             <- 1.96
sigma_Y        <- 1
shift          <- cv
L              <- 6.5    #grid support [-L, L] plus 9999
numgrid        <- 3000   #number of grid points
boots          <- 1000   #bootstrap repetitions (projection test and Kudrin's T_inf)
boot_seed_by_draw <- FALSE #FALSE: the same bootstrap seeds in every de-rounding draw, as application_cochrane.R, so a
#test's p-value moves across draws only through the de-rounding; TRUE: new seeds in every
#draw (the earlier runs)
J              <- 40
p              <- 2
J0_list        <- c(20, 30) #J0 = 20 is the main specification; J0 = 30 is the robustness check (R5 3iii)
C_round        <- c("20" = 0.13, "30" = 0.31)   #C of Propositions 3 and 4, by J0 (for p = 2)
run_stouffer   <- TRUE   #p-curve test for left skew by Stouffer's method (Simonsohn, point 5)
ekw_failed_p   <- NA_real_ #p-value an EKW test gets in a draw where it fails (an error, or Cox-Shi's codes 888/999;
#see ekw_safe.R): NA leaves that draw out of the test's averages, and EWK_failures<suffix>
#counts the failures; 1 counts it as a non-rejection, as kudrin_figures_engine.R does
run_noshift    <- TRUE   #projection test without the 1.96 shift, J0 = 20 (R1 17a)
noshift_draws  <- 200    #de-roundings for the unshifted test under options 2 and 4 (at most n_draws_noise)
time_tests     <- TRUE   #time one call of each test per design, for a run-time sentence (R5 3iv)

# Figures (R1 17b): their data are computed here, and drawn by application_MM_output.R
fig_resid_option <- 3    #the residual figure's de-rounding option (3 = smin, deterministic)
fig_resid_boots  <- 200  #article-level bootstrap draws for its 90% bands (0: no bands)
fig_shape_option <- 4    #the reported-versus-honest figure's option (4 = smin + noise, its first draw: the noise smooths
#the rounding lattice of smin t-statistics, which a histogram would otherwise show)
fig_binwidth     <- 0.1  #its bins of |t|
fig_max_t        <- 5    #its right end
fig_bandwidth    <- 0.25 #Gaussian kernel of its curvature check (0.35 as a robustness check)
fig_window       <- c(1.6, 2.6)   #where the curvature check's simultaneous band applies
fig_curv_boots   <- 500  #article-level bootstrap draws for the curvature check
fig_level        <- 0.95
fig_seed         <- 20261030      #seeds of the two figures' bootstraps
fig_shape_seed   <- 20261102

run_tag <- paste0(if (drop_zero_mu) "" else "_keepmu0", if (require_se) "" else "_keept",
                  if (drop_duplicates) "_dedup" else "", run_label)   #in every file name of this run

#### 1. Code ---------------------------------------------------------------------------------------------
source(file.path(root, "code", "functions_all_weighted.R"))
if (!requireNamespace("nnls", quietly = TRUE)) stop("The figure data need the nnls package: install.packages(\"nnls\")")
coeffs        <- J + 1
epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U   #E(L, delta), included in every p-value and B-hat
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

#### 2. Data: Brodeur et al. (2020), with or without Kudrin's (2025) drop of zero coefficients -----------
MM_raw    <- haven::read_dta(file.path(data_location, "MM data.dta"))
n_zero_mu <- sum(MM_raw$mu %in% 0)
drop      <- is.na(MM_raw$mu) | (drop_zero_mu & MM_raw$mu %in% 0)
cat(sprintf("Coefficients printed as zero: %d, %s.\n", n_zero_mu,
            if (drop_zero_mu) "dropped, following Kudrin (2025)" else "kept"))
MM_data   <- MM_raw[!drop, ]

#### 3. De-round the t-statistics (Section 5.2, Propositions 3 and 4)
# Four de-rounding options. For a t-statistic computed from a coefficient m and a standard error s printed to
# multiples of a and b, with s_min = s - b/2 and s_max = s + b/2 the smallest and largest standard errors
# consistent with s:
#   1 "none":       t = m / s
#   2 "uniform":    t = (m + V) / W,  V ~ U(-a/2, a/2), W ~ U(s_min, s_max)  (the usual de-rounding)
#   3 "smin":       t = m / s_min                                              (Section 5.2)
#   4 "smin_noise": t = (m + V) / s_min
# With require_se, every test in the sample is a coefficient and a standard error. Some papers print a t-statistic
# or a p-value in their place; with require_se = FALSE, a printed t-statistic is used as printed under options 1
# and 3, and gets noise uniform on +-1/2 of its rounding unit a_t under options 2 and 4. Never covered:
# t-statistics converted from rounded p-values, which are dropped either way. Options 2 and 4 are run on
# n_draws_noise independent draws, and p-values and breakdown statistics are averaged over them.
#
# The sample is the same under every option: which t-statistics are dropped depends only on the printed
# values and kappa_cap, never on the option. The noise is drawn once, with fixed seeds, so options 2
# and 4 use identical numerator noise and differ only in the denominator.
#
# Rounding units:
#   b, per entry: the place value of sd's last digit. A dropped trailing zero only makes b too large,
#     which keeps s_min a lower bound; a table-wide b can be too small for a column printed more coarsely.
#   a, per table: the finest place value among the table's coefficients, but no finer than b, since a
#     coefficient and its standard error are almost always printed to the same decimals. Both are fixed
#     by the table's format and the standard error, not by the estimate. The width of V need not match
#     the true unit: Proposition 4 allows any independent noise.
#   Used only with require_se = FALSE:
#   a_t, per article, for a t printed without a standard error: the finest place value among the article's
#     printed t-statistics. A t has no units, so its decimals are the article's style, and the file has lost
#     trailing zeros: a t printed 2.00 is stored as 2, and its own last digit would give kappa = 1, dropping
#     it for its value. Taking the unit from the article keeps it fixed by the format, as a is.
#   A t with 4 or more decimals was computed by the file, not printed: converted back, it gives a p-value with
#     few decimals (the sample report shows them), so it came from a printed p-value. It never sets an article's
#     a_t, and is dropped: rounding a p-value is not rounding t on a lattice (Propositions 3 and 4).

stopifnot(all(options_to_run %in% 1:4))
deround_keys   <- c("none", "uniform", "smin", "smin_noise")
deround_labels <- c("No de-rounding", "Uniform de-rounding", "Smallest consistent standard error",
                    "Smallest consistent standard error with uniform noise")
deround_noise  <- c(FALSE, TRUE, FALSE, TRUE)

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

has_ratio  <- !is.na(mu_p) & !is.na(sd_p)
direct     <- !has_ratio & !is.na(t_p)                       #a t-statistic in the file, and no standard error
pvalue_t   <- direct & dec_t >= 4                            #4+ decimals: computed by the file, from a printed p-value
printed_t  <- direct & !pvalue_t
a_art      <- 10^(-ave(ifelse(printed_t, dec_t, NA), MM_data$title,   #the finest printed t in the article: a t
                       FUN = function(d) if (all(is.na(d))) NA_real_ else max(d, na.rm = TRUE)))   #printed 2.00 is stored as 2
a_t        <- ifelse(printed_t, a_art, 10^(-dec_t))          #rounding unit of a t without a standard error
p_decimals <- n_dec(round(2 * pnorm(-abs(MM_data$t)), 6))   #decimals of the p-value each t converts back to: few
pval_evidence <- table(design = MM_data$method[pvalue_t], p_decimals = p_decimals[pvalue_t])   #mean a printed p
no_se_kinds   <- table(design = MM_data$method[direct],                                         #what papers printed
                       printed = ifelse(pvalue_t, "p-value", "t-statistic")[direct])           #in place of the SE

drop_nose  <- !has_ratio & require_se                        #no standard error printed (a t or a p-value instead)
drop_se0   <- has_ratio & sd_p <= 0                          #SE printed as zero: t not recoverable
drop_kappa <- (has_ratio & !drop_se0 & a_mu / s_min > kappa_cap) | (printed_t & !require_se & a_t > kappa_cap)
drop_pval  <- pvalue_t & !require_se                         #t computed from a printed p-value (robustness run)
use_ratio  <- has_ratio & !drop_se0 & !drop_kappa

# Exact duplicates: a later row repeating an earlier row's design, article, printed coefficient, standard error and t
# (the key the sample report counts as exact_duplicates). The rule applies after every other rule, among the rows
# they keep, so a copy that survives them is never dropped in favor of one that does not: a coefficient's rounding
# unit is set per table, so two copies of a test can differ in kappa. The first copy in the file is kept, as
# application_cochrane.R keeps the first copy of each trial's outcome data.
keep_other <- !(drop_nose | drop_se0 | drop_kappa | drop_pval)   #the sample before duplicates are considered
dup_later  <- rep(FALSE, n_obs)
dup_later[keep_other] <- duplicated(paste(MM_data$method, MM_data$title, mu_p, sd_p, t_p, sep = "|")[keep_other])
drop_dup   <- drop_duplicates & dup_later

cat("Decimals printed in the standard errors (many at 6+ would suggest computed rather than printed SEs):\n")
print(table(dec_sd[has_ratio]))
cat("\nt-statistics dropped (the same under every de-rounding option):\n")
drop_counts <- rbind(se_missing      = tapply(drop_nose,  MM_data$method, sum),
                     se_zero         = tapply(drop_se0,   MM_data$method, sum),
                     kappa_above_cap = tapply(drop_kappa, MM_data$method, sum),
                     t_from_pvalue   = tapply(drop_pval,  MM_data$method, sum))
if (drop_duplicates) drop_counts <- rbind(drop_counts, exact_duplicate = tapply(drop_dup, MM_data$method, sum))
print(drop_counts)

#### How much do we drop, compared with Kranz and Putz's (2022) omission rule?
# Runs before MM_data is subset, while the vectors above still line up with the full sample.
#
# KP omit every test whose standard error has a significand below 37: the SE's significant digits read
# as an integer, e.g. 0.012 -> 12. With the coefficient and SE printed to the same decimals, this keeps
# rounding from moving any z-statistic across 1.96 unless it lies within 0.04 of it (their Lemma 1).
# Our data have lost trailing zeros, so we use KP's own fallback for that case: append zeros to the SE
# until it has as many decimals as its coefficient (KP online appendix, note to Figure A2), which they
# report agrees with the exact rule in 96 percent of cases. If your file records trailing zeros, read
# the significand from those instead. Both rules act on tests printed as a coefficient and a standard error,
# so the comparison is on those rows; tests printed without a standard error are left out of it.

kp_threshold <- 37
kp_signif    <- round(sd_p * 10^pmax(dec_sd, dec_mu))                   #SE significand, KP's convention
drop_kp      <- has_ratio & kp_signif < kp_threshold                     #includes SEs printed as zero
drop_ours    <- has_ratio & (drop_se0 | drop_kappa)                   #rounding drops only: duplicates are not part of it
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
rows_kp         <- c(setNames(lapply(methods_kp, function(m) has_ratio & MM_data$method %in% m), methods_kp),
                     list(Pooled = has_ratio))                         #rows printed with a standard error
drop_comparison <- t(sapply(rows_kp, drop_summary))

cat("\nTests printed with a standard error, dropped by our procedure (SE = 0 or kappa > ", kappa_cap,
    ") and by Kranz and Putz's rule (SE significand < ", kp_threshold, "):\n", sep = "")
print(round(drop_comparison, 1))
cat(sprintf("Of the tests printed with a standard error, our procedure drops %.1f to %.1f percent by method (%.1f pooled); Kranz and Putz's rule would drop %.1f to %.1f percent (%.1f pooled).\n",
            min(drop_comparison[methods_kp, "ours_pct"]), max(drop_comparison[methods_kp, "ours_pct"]), drop_comparison["Pooled", "ours_pct"],
            min(drop_comparison[methods_kp, "kp_pct"]),   max(drop_comparison[methods_kp, "kp_pct"]),   drop_comparison["Pooled", "kp_pct"]))
cat(sprintf("Of the %d t-statistics at exactly z = 2, our procedure drops %.1f percent and theirs %.1f percent. (KP report 37.9 percent of all observations and 87.3 percent of those at z = 2, on their version of the data.)\n",
            drop_comparison["Pooled", "z_eq_2"], drop_comparison["Pooled", "ours_z2_pct"], drop_comparison["Pooled", "kp_z2_pct"]))
saveRDS(drop_comparison, file = file.path(res_dir, paste0("Drop_comparison_KP", run_tag)))

#The row-level record of the sample: every row of the file, kept or not, and why (for the sample report)
strip       <- function(x) as.vector(haven::zap_labels(x))
kept_rows   <- which(!drop)                                     #rows of MM_raw that are rows of MM_data here
sample_rows <- data.frame(row = seq_len(nrow(MM_raw)), article = strip(MM_raw$title), table = strip(MM_raw$table),
                          design = strip(MM_raw$method), coefficient = strip(MM_raw$mu), se = strip(MM_raw$sd),
                          t = strip(MM_raw$t), stringsAsFactors = FALSE)
reason <- rep("", nrow(MM_raw))
reason[is.na(MM_raw$mu)] <- "coefficient missing"
if (drop_zero_mu) reason[!is.na(MM_raw$mu) & MM_raw$mu %in% 0] <- "coefficient printed as zero"
reason[kept_rows[drop_nose]]  <- "standard error missing"
reason[kept_rows[drop_se0]]   <- "standard error printed as zero"
reason[kept_rows[drop_kappa]] <- "kappa above cap"
reason[kept_rows[drop_pval]]  <- "t computed from a p-value"
reason[kept_rows[drop_dup]]   <- "exact duplicate"
sample_rows$kept     <- reason == ""
sample_rows$reason   <- reason
sample_rows$t_source <- NA_character_
sample_rows$t_source[kept_rows] <- ifelse(use_ratio, "coef/SE", ifelse(pvalue_t, "t from a p-value", ifelse(direct, "reported t", NA)))
sample_rows$kappa    <- NA_real_
sample_rows$kappa[kept_rows]    <- ifelse(has_ratio & !drop_se0, a_mu / s_min, ifelse(direct, a_t, NA))

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
MM_data            <- MM_data[keep_other, ]

#The de-rounding noise, drawn once with fixed seeds, so every option and every run uses the same draws. It is drawn
#before duplicates are dropped, so each remaining test gets the same noise whether or not drop_duplicates is set
set.seed(1)
U_num <- matrix(runif(nrow(MM_data) * n_draws_noise, -0.5, 0.5), ncol = n_draws_noise)   #numerator noise, one column per draw
set.seed(2)
U_den <- matrix(runif(nrow(MM_data) * n_draws_noise, -0.5, 0.5), ncol = n_draws_noise)   #denominator noise (option 2 only)
if (drop_duplicates) {
  first_copy <- !dup_later[keep_other]
  MM_data    <- MM_data[first_copy, ]
  U_num      <- U_num[first_copy, , drop = FALSE]
  U_den      <- U_den[first_copy, , drop = FALSE]
}
stopifnot(nrow(MM_data) == sum(sample_rows$kept))   #the sample report and the tests use the same rows

deround_draw <- function(r) {                                #t-statistics for draw r under the current option
  d     <- MM_data                                           #(reads deround_spec and add_noise, set in the loop below)
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

cat("\nt-statistics by source:\n")
print(table(MM_data$method, MM_data$t_source, useNA = "ifany"))
cat("\nLargest possible inflation of t under s_min, s_max/s_min (median, 99th percentile, max):\n")
print(round(sapply(split(MM_data$lambda_max, MM_data$method), quantile, probs = c(.5, .99, 1), na.rm = TRUE), 3))
kappa6 <- tapply(MM_data$kappa^6, MM_data$method, mean, na.rm = TRUE)
cat("\nRounding allowance under options 3 and 4 (Section 5.2.1, step 3), to compare with B-hat:\n")
print(signif(rbind(mean_kappa6 = kappa6, bound_J0_20 = C_round[["20"]] * kappa6, bound_J0_30 = C_round[["30"]] * kappa6), 3))

#### Sample report: every step from the file to the final sample, by design, reconciled with Kudrin (2025) ----
meth    <- c("RCT","IV","DID","RDD")
designs <- c(meth, setdiff(sort(unique(sample_rows$design[!is.na(sample_rows$design)])), meth))
count_by <- function(v) {                                     #rows, by design and in total
  x <- table(factor(sample_rows$design[v], levels = designs))
  c(setNames(as.vector(x), designs), Total = sum(v))
}
papers_by <- function(v) c(vapply(designs, function(m) length(unique(sample_rows$article[v & sample_rows$design %in% m])), 0),
                           Total = length(unique(sample_rows$article[v])))
is_zero  <- !is.na(sample_rows$coefficient) & sample_rows$coefficient %in% 0
steps_mm <- rbind("In the file"                         = count_by(rep(TRUE, nrow(sample_rows))),
                  "Dropped: coefficient missing"        = count_by(reason == "coefficient missing"),
                  "Dropped: coefficient printed as 0"   = count_by(reason == "coefficient printed as zero"),
                  "Dropped: standard error missing"     = count_by(reason == "standard error missing"),
                  "Dropped: standard error printed as 0" = count_by(reason == "standard error printed as zero"),
                  "Dropped: kappa above the cap"        = count_by(reason == "kappa above cap"),
                  "Dropped: t computed from a p-value"  = count_by(reason == "t computed from a p-value"),
                  "Dropped: exact duplicate"            = count_by(reason == "exact duplicate"),
                  "Final sample"                        = count_by(sample_rows$kept),
                  "Articles in the final sample"        = papers_by(sample_rows$kept))
if (!drop_zero_mu) steps_mm <- steps_mm[rownames(steps_mm) != "Dropped: coefficient printed as 0", ]
if (!drop_duplicates) steps_mm <- steps_mm[rownames(steps_mm) != "Dropped: exact duplicate", ]
steps_mm <- steps_mm[rownames(steps_mm) != (if (require_se) "Dropped: t computed from a p-value" else "Dropped: standard error missing"), ]

# Kudrin (2025), Table 2: tests (papers) after excluding observations with missing information. Our reading of his
# rule: a coefficient is present and not printed as zero. Our sample differs from it by exactly
# (zero coefficients we keep) - (standard errors missing) - (standard errors printed as zero) - (kappa above the cap),
# or, with require_se = FALSE, - (t's from p-values) instead of - (standard errors missing); with drop_duplicates,
# also - (exact duplicates)
kud_n      <- c(RCT = 7101, IV = 5158, DID = 5780, RDD = 3117)
kud_papers <- c(RCT = 141, IV = 281, DID = 241, RDD = 85)
kud_rule   <- !is.na(sample_rows$coefficient) & !is_zero
recon <- rbind("Kudrin (2025), Table 2: tests"              = c(kud_n[meth], Total = sum(kud_n)),
               "Kudrin (2025), Table 2: papers"             = c(kud_papers[meth], Total = sum(kud_papers)),
               "His rule on this file: tests"               = count_by(kud_rule)[c(meth, "Total")],
               "His rule on this file: papers"              = c(papers_by(kud_rule)[meth], Total = sum(papers_by(kud_rule)[meth])),
               "+ zero coefficients we keep"                = if (drop_zero_mu) 0 * count_by(is_zero)[c(meth, "Total")] else count_by(is_zero)[c(meth, "Total")],
               "- standard errors missing"                  = count_by(reason == "standard error missing")[c(meth, "Total")],
               "- standard errors printed as zero"          = count_by(reason == "standard error printed as zero")[c(meth, "Total")],
               "- kappa above the cap"                      = count_by(reason == "kappa above cap")[c(meth, "Total")],
               "- t computed from a p-value"                = count_by(reason == "t computed from a p-value")[c(meth, "Total")],
               "- exact duplicates"                         = count_by(reason == "exact duplicate")[c(meth, "Total")],
               "= our final sample"                         = count_by(sample_rows$kept)[c(meth, "Total")])
recon <- recon[rownames(recon) != (if (require_se) "- t computed from a p-value" else "- standard errors missing"), ]
if (!drop_duplicates) recon <- recon[rownames(recon) != "- exact duplicates", ]

t_print <- ifelse(MM_data$t_source %in% "coef/SE", MM_data$dr_mu / MM_data$dr_sd, MM_data$dr_t)
dup_key <- paste(MM_data$method, MM_data$title, MM_data$dr_mu, MM_data$dr_sd, MM_data$dr_t, sep = "|")
final_mm <- t(vapply(meth, function(m) {
  s <- MM_data$method %in% m
  k <- MM_data$kappa[s]
  c(tests = sum(s), articles = length(unique(MM_data$title[s])),
    coef_over_SE = sum(MM_data$t_source[s] %in% "coef/SE"), reported_t = sum(MM_data$t_source[s] %in% "reported t"),
    share_sig = mean(abs(t_print[s]) > cv, na.rm = TRUE), zero_coef_kept = sum(MM_data$dr_mu[s] %in% 0),
    exact_duplicates = sum(duplicated(dup_key[s])),
    kappa_median = median(k, na.rm = TRUE), kappa_p90 = unname(quantile(k, 0.9, na.rm = TRUE)), kappa_max = max(k, na.rm = TRUE),
    allowance_J0_20_x1000 = 1e3 * C_round[["20"]] * kappa6[[m]], allowance_J0_30_x1000 = 1e3 * C_round[["30"]] * kappa6[[m]])
}, numeric(12)))

show_lines <- function(x, digits = NULL) {                     #a printed table as report lines
  op <- options(width = 250); on.exit(options(op))
  out <- capture.output(print(if (is.null(digits)) x else round(x, digits)))
  paste0("    ", out[nzchar(trimws(out))])
}
report_mm <- c(
  "SAMPLE REPORT: the application to Brodeur, Cook and Heyes (2020)",
  sprintf("Written %s by application_MM_analysis.R; run tag '%s'. Data: %s (%s rows, %s articles).",
          format(Sys.time(), "%Y-%m-%d %H:%M"), run_tag, file.path(data_location, "MM data.dta"),
          format(nrow(MM_raw), big.mark = ","), format(length(unique(sample_rows$article)), big.mark = ",")),
  "One row of the file is one test: a coefficient and its standard error, or a reported t-statistic.",
  "",
  "1. SETTINGS THAT DEFINE THE SAMPLE",
  sprintf("  drop_zero_mu = %s (coefficients printed as zero %s); kappa_cap = %s; se_rounded_up = %s; require_se = %s;",
          drop_zero_mu, if (drop_zero_mu) "dropped, as Kudrin (2025)" else "kept", format(kappa_cap), se_rounded_up, require_se),
  sprintf("  drop_duplicates = %s (%s)", drop_duplicates,
          if (drop_duplicates) "a later copy of a test the file lists more than once is dropped, after every other rule" else
            "tests the file lists more than once are kept; section 3 counts them"),
  "  kappa = a / s_min: the coefficient's rounding unit a (the finest place value among its table's coefficients, but no",
  "  finer than its standard error's) over the smallest standard error consistent with the printed one, s - b/2.",
  "  Some tests have a coefficient but no standard error: the paper printed a t-statistic or a p-value in its place, and",
  "  the file stores a t. One computed from a p-value has 4 or more decimals; converted back, it gives a p-value with the",
  "  number of decimals in the second table (few decimals mean a printed p-value).",
  show_lines(no_se_kinds),
  show_lines(pval_evidence),
  if (require_se) "  require_se = TRUE drops all of them, so every t-statistic in the sample is coefficient / standard error." else
    c("  require_se = FALSE keeps the printed t-statistics, with kappa = a_t, the finest place value among the article's",
      "  printed t-statistics (the file has lost trailing zeros, so a t printed 2.00 is stored as 2), under the same cap.",
      "  It drops those computed from p-values: rounding a p-value is not rounding t on a lattice (Propositions 3 and 4)."),
  "  Which rows are dropped never depends on the de-rounding option.",
  "",
  "2. EVERY STEP, IN ORDER, BY DESIGN (rows dropped by each step; the first and last rows count rows kept)",
  show_lines(steps_mm),
  sprintf("  Rows of the final sample with no usable t-statistic: %d (must be 0).", sum(sample_rows$kept & is.na(sample_rows$t_source))),
  "",
  "3. THE FINAL SAMPLE, BY DESIGN",
  show_lines(final_mm, 4),
  "  share_sig: share with printed |t| > 1.96. exact_duplicates: rows repeating an earlier row's article, coefficient,",
  sprintf("  standard error and t (%s). Allowances: C mean(kappa^6) x 1000, C = 0.13 (J0 = 20), 0.31 (J0 = 30).",
          if (drop_duplicates) "0 by construction: drop_duplicates dropped them, as section 2 shows" else "kept, as Kudrin keeps them"),
  "",
  "4. RECONCILIATION WITH KUDRIN (2025)",
  "  His Table 2 reports the tests left 'after excluding observations with missing information'. We read his rule as:",
  "  the coefficient is present and not printed as zero. Our final sample differs from that rule's sample by exactly the",
  "  rows below it:",
  show_lines(recon),
  if (any(recon["His rule on this file: tests", meth] != kud_n[meth]))
    "  NOTE: our reading of his rule does not reproduce his counts for every design; the gap is in the rows above." else
      "  Our reading of his rule reproduces his counts for every design.",
  "  Other differences in the analysis, not the sample: every test here clusters by article (Kudrin's T_inf resamples",
  "  p-values independently), and every test is run under four de-rounding options (Section 5.2).",
  "",
  "5. COMPARISON WITH KRANZ AND PUTZ'S (2022) OMISSION RULE (standard-error significand below 37), ON THE TESTS PRINTED",
  "   WITH A STANDARD ERROR (the only ones either rule acts on)",
  show_lines(drop_comparison, 1),
  "",
  sprintf("6. ROW-LEVEL FILE: Sample_rows_application%s.csv lists every row of the file (row = its position in MM data.dta)", run_tag),
  "  with its article, table, design, coefficient, standard error and t as stored, whether it is kept, the reason if not,",
  "  how its t-statistic is formed, and its kappa.")
writeLines(report_mm, file.path(res_dir, paste0("Sample_report_application", run_tag, ".txt")))
write.csv(sample_rows, file.path(res_dir, paste0("Sample_rows_application", run_tag, ".csv")), row.names = FALSE)
cat("\n", paste(report_mm, collapse = "\n"), "\n\nSample report written to ",
    file.path(res_dir, paste0("Sample_report_application", run_tag, ".txt")), "\n", sep = "")

#Averaging over draws
as_num1         <- function(v) as.numeric(v)[1]              #test results may come back as lists or list-columns
draw_range      <- function(x) { x <- x[is.finite(x)]; if (length(x)) diff(range(x)) else NA_real_ }   #max - min across draws
summarize_draws <- function(d, cols, f = mean) {             #one row per method, in the order of meth;
  out <- d[match(meth, d$methods), setdiff(names(d), "draw")] #cols summarized over draws by f
  for (k in cols) out[[k]] <- as.numeric(tapply(d[[k]], d$methods, f)[meth])
  rownames(out) <- NULL
  out
}


#### 4. The run's settings, saved now (status "running") and again at the end (status "complete") -------
settings_now <- function(status) list(
  status = status, saved = format(Sys.time()), R = R.version.string, script = "application_MM_analysis.R",
  tests_from = if (exists("tests_source")) tests_source else "this run",
  run_tag = run_tag, run_label = run_label, data_file = file.path(data_location, "MM data.dta"),
  drop_zero_mu = drop_zero_mu, n_zero_mu = n_zero_mu, kappa_cap = kappa_cap, se_rounded_up = se_rounded_up,
  require_se = require_se, drop_duplicates = drop_duplicates, no_se_kinds = no_se_kinds, pval_evidence = pval_evidence,
  drop_counts = drop_counts,
  n_by_design        = setNames(as.vector(table(factor(MM_data$method, levels = meth))), meth),
  articles_by_design = vapply(meth, function(m) length(unique(MM_data$title[MM_data$method == m])), 0),
  kappa6    = setNames(as.vector(kappa6[meth]), meth), C_round = C_round,
  allowance = lapply(C_round, function(C) setNames(C * as.vector(kappa6[meth]), meth)),
  options_run = options_to_run, keys_run = deround_keys[options_to_run], deround_keys = deround_keys,
  deround_labels = deround_labels, deround_noise = deround_noise, n_draws_noise = n_draws_noise,
  n_draws = setNames(ifelse(deround_noise[options_to_run], n_draws_noise, 1), deround_keys[options_to_run]),
  noshift_n_draws = setNames(ifelse(deround_noise[options_to_run], min(noshift_draws, n_draws_noise), 1),
                             deround_keys[options_to_run]),
  cv = cv, sigma_Y = sigma_Y, shift = shift, L = L, numgrid = numgrid, boots = boots, boot_seed_by_draw = boot_seed_by_draw,
  J = J, coeffs = coeffs, p = p,
  J0_list = J0_list, epsilon_U_pre = epsilon_U_pre, run_stouffer = run_stouffer, run_noshift = run_noshift,
  ekw_failed_p = ekw_failed_p,
  noshift_draws = noshift_draws, time_tests = time_tests,
  fig_resid_option = fig_resid_option, fig_resid_boots = fig_resid_boots, fig_shape_option = fig_shape_option,
  fig_binwidth = fig_binwidth, fig_max_t = fig_max_t, fig_bandwidth = fig_bandwidth, fig_window = fig_window,
  fig_curv_boots = fig_curv_boots, fig_level = fig_level, fig_seed = fig_seed, fig_shape_seed = fig_shape_seed,
  deff = if (exists("deff_by_design")) deff_by_design else NULL)
settings_file <- file.path(res_dir, paste0("Settings_application", run_tag))
save_settings <- function(status) saveRDS(settings_now(status), settings_file)

tests_source <- "this run"
if (!run_tests) {                       #reuse saved test results only if they were run with these settings
  now <- settings_now("check")
  if (file.exists(settings_file)) {
    old  <- readRDS(settings_file)
    if (is.null(old$drop_duplicates)) old$drop_duplicates <- FALSE   #saved before the setting existed: duplicates were kept
    same <- c("drop_zero_mu", "drop_duplicates", "kappa_cap", "se_rounded_up", "require_se", "options_run", "n_draws_noise", "cv", "shift", "L",
              "numgrid", "boots", "boot_seed_by_draw", "J", "p", "J0_list", "run_stouffer", "run_noshift", "noshift_draws",
              "n_by_design")
    diff <- same[!vapply(same, function(k) isTRUE(all.equal(old[[k]], now[[k]])), TRUE)]
    if (length(diff)) stop("run_tests = FALSE, but these settings differ from the saved run's: ", paste(diff, collapse = ", "))
    if (!identical(old$status, "complete")) stop("The saved run with run_tag '", run_tag, "' did not finish")
    tests_source <- paste("the run saved", old$saved)
  } else {                              #results saved before settings files existed: adopt them if their sample matches
    f <- file.path(res_dir, paste0("MM_weights_results_20_J0_", deround_keys[options_to_run], run_tag))
    if (!all(file.exists(f))) stop("run_tests = FALSE, but there are no saved results with run_tag '", run_tag, "'")
    for (x in f) {
      r <- readRDS(x)
      n_file <- setNames(vapply(r$n, function(v) as.numeric(v)[1], 0), r$methods)[meth]
      if (!isTRUE(all(n_file == now$n_by_design)))
        stop(basename(x), " was run on a different sample: n = ", paste(n_file, collapse = "/"), ", not ",
             paste(now$n_by_design, collapse = "/"))
    }
    warning("No settings file for run_tag '", run_tag, "': adopting its saved results, which match this sample, and ",
            "assuming they were run with the test settings above", call. = FALSE)
    tests_source <- "an earlier run without a settings file (adopted; same sample)"
  }
  cat("Keeping the test results of ", tests_source, "\n", sep = "")
} else {
  save_settings("running")
}

#### 5. The tests ---------------------------------------------------------------------------------------
### EKW and Kudrin (2025) test wrappers

source(paste0(root,"/code/Tests.R"))
source(paste0(root,"/code/ekw_safe.R"))   #ekw_six(), mean_finite(), min_finite(): a failed EKW test is never averaged in

run_ewk <- function(ps, ids = 1:length(ps), pmax = 0.15) {   #a test that fails gets ekw_failed_p (ekw_safe.R)
  e   <- ekw_six(ps, ids, pmax, failed_p = ekw_failed_p)
  out <- list(lcms_EWK = e$p[["LCM"]], Fisher_EWK = e$p[["Fisher"]], disconts_EWK = e$p[["Disc"]],
              binomial_EWK = e$p[["Binomial"]], CS1_EWK_s = e$p[["CS1"]], CS2_EWK_s = e$p[["CS2B"]])
  out$min <- min_finite(e$p)                                  #smallest p-value among the tests that ran
  attr(out, "status") <- setNames(e$status, names(out)[1:6])  #0 if the test ran; else why not (ekw_safe.R)
  out
}

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

#p-curve test for left skew (Simonsohn, Nelson and Simmons 2014; Simonsohn's report, point 5): every p < .05
#is multiplied by 20, transformed by the normal quantile function and combined by Stouffer's method; the
#one-sided p-value is for left skew. Like p-curve it treats the p-values as independent. qnorm(20 p) is
#computed on the log scale, so huge |t| do not underflow to -Inf.
stouffer_p <- function(t) {
  a   <- abs(t)
  sig <- is.finite(a) & 2 * pnorm(-a) < 0.05
  if (!any(sig)) return(NA_real_)
  z   <- qnorm(log(40) + pnorm(-a[sig], log.p = TRUE), log.p = TRUE)
  pnorm(sum(z) / sqrt(sum(sig)), lower.tail = FALSE)
}


#### Run every test under every de-rounding option
for (deround_option in if (run_tests) options_to_run else integer(0)) {
  deround_spec <- list(key   = deround_keys[deround_option],
                       label = deround_labels[deround_option],
                       noise = deround_noise[deround_option])
  add_noise <- deround_spec$noise
  suffix    <- paste0("_", deround_spec$key, run_tag)        #labels every results file and the table (run_tag: see drop_zero_mu)
  n_draws   <- if (add_noise) n_draws_noise else 1           #without noise every draw would be the same
  cat("\n========== De-rounding option ", deround_option, ": ", deround_spec$label, " ==========\n", sep = "")
  
  if (any(!is.finite(deround_draw(1))) | any(!is.finite(MM_data$mu))){
    stop("Tscores not clean")
  }
  
  #Projection test, on each de-rounding
  proj_all <- NULL
  for (r in seq_len(n_draws)) {
    t_start   <- Sys.time()
    MM_data$t <- deround_draw(r)
    set.seed(if (boot_seed_by_draw) 1000 * r + 1 else 1)   #see boot_seed_by_draw
    for (J0 in J0_list) {
      for (m in meth) {
        res <- run_test_weight_fast(MM_data[MM_data$method==m,],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
        proj_all <- rbind(proj_all, data.frame(draw = r, methods = m, n = as_num1(res$n), articles = as_num1(res$articles),
                                               pval = as_num1(res$pval), breakdown = as_num1(res$breakdown),
                                               numcoeffs = as_num1(res$numcoeffs), J0 = J0, shift = shift))
      }
    }
    if (run_noshift && r <= noshift_draws) {        #the same test without the 1.96 shift (R1 17a)
      for (m in meth) {
        res <- run_test_weight_fast(MM_data[MM_data$method==m,],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=0,L=L,numgrid=numgrid,boots=boots,J0=20,p=p)
        proj_all <- rbind(proj_all, data.frame(draw = r, methods = m, n = as_num1(res$n), articles = as_num1(res$articles),
                                               pval = as_num1(res$pval), breakdown = as_num1(res$breakdown),
                                               numcoeffs = as_num1(res$numcoeffs), J0 = 20, shift = 0))
      }
    }
    cat(sprintf("Projection test (%s): de-rounding %d of %d done in %.1f minutes\n", deround_spec$key, r, n_draws,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
  
  #Save the averages over draws (the files application_MM_output.R reads) and every draw, labeled by option
  saveRDS(proj_all, file = file.path(res_dir, paste0("MM_weights_results_draws", suffix)))
  for (J0 in J0_list) {
    name <- file.path(res_dir, paste0("MM_weights_results_",J0,"_J0", suffix))
    saveRDS(summarize_draws(proj_all[proj_all$J0 == J0 & proj_all$shift == shift, ], c("pval", "breakdown")), file = name)
    print(readRDS(name))
  }
  if (run_noshift) {                                 #the unshifted test, for Table 3 (R1 17a)
    name <- file.path(res_dir, paste0("MM_weights_results_20_J0_noshift", suffix))
    saveRDS(summarize_draws(proj_all[proj_all$J0 == 20 & proj_all$shift == 0, ], c("pval", "breakdown")), file = name)
    print(readRDS(name))
  }
  
  #EKW tests, on each de-rounding. A test that fails in a draw (an error, or Cox-Shi's codes 888/999) gets
  #ekw_failed_p there: with NA, its averages leave that draw out. EWK_failures<suffix> counts the failures
  ewk_draws <- vector("list", n_draws)
  ewk_codes <- vector("list", n_draws)
  for (r in seq_len(n_draws)) {
    MM_data$t <- deround_draw(r)
    set.seed(if (boot_seed_by_draw) 1000 * r + 2 else 2)
    ewk_r <- lapply(setNames(meth, paste0("ewk_", tolower(meth))), function(m) {
      s <- MM_data$method == m
      run_ewk(2 * pnorm(-abs(MM_data$t[s])), ids = as.integer(factor(MM_data$title[s])))
    })
    ewk_draws[[r]] <- sapply(ewk_r, function(x) vapply(x, as_num1, numeric(1)))   #tests x methods
    ewk_codes[[r]] <- sapply(ewk_r, function(x) attr(x, "status"))               #the six tests x methods; 0 = ran
  }
  ewk_all       <- simplify2array(ewk_draws)                 #tests x methods x draws
  ewk_code_all  <- simplify2array(ewk_codes)                 #the six tests x methods x draws
  empirical_ewk <- apply(ewk_all, c(1, 2), mean_finite)      #each test's mean over the draws in which it ran
  ewk_tests     <- setdiff(rownames(empirical_ewk), "min")
  empirical_ewk["min", ] <- apply(empirical_ewk[ewk_tests, , drop = FALSE], 2, min_finite)   #smallest averaged p-value
  ewk_failed    <- apply(ewk_code_all != 0L, c(1, 2), sum)   #draws in which each test failed, tests x methods
  saveRDS(ewk_all, file = file.path(res_dir, paste0("EWK_results_draws", suffix)))
  saveRDS(empirical_ewk, file = file.path(res_dir, paste0("EWK_results", suffix)))
  saveRDS(list(failed = ewk_failed, codes = ewk_code_all, failed_p = ekw_failed_p, draws = n_draws),
          file = file.path(res_dir, paste0("EWK_failures", suffix)))
  print(round(empirical_ewk, 3))
  if (any(ewk_failed > 0)) {
    cat("EKW tests that failed, in how many of the ", n_draws, " draws (their p-value is ", format(ekw_failed_p),
        " there):\n", sep = "")
    print(ewk_failed)
  }
  
  #Stouffer's p-curve test for left skew, on each de-rounding (Simonsohn's report, point 5)
  if (run_stouffer) {
    st_all <- t(vapply(seq_len(n_draws), function(r) {
      tt <- deround_draw(r)
      vapply(meth, function(m) stouffer_p(tt[MM_data$method == m]), numeric(1))
    }, numeric(length(meth))))                                 #draws x methods
    st_res <- data.frame(methods = meth, pval = colMeans(st_all), stringsAsFactors = FALSE)
    saveRDS(st_all, file = file.path(res_dir, paste0("Stouffer_results_draws", suffix)))
    saveRDS(st_res, file = file.path(res_dir, paste0("Stouffer_results", suffix)))
    print(st_res)
  }
  
  #Kudrin (2025) tests: T_inf, T_H, T_H_hybrid, on each de-rounding
  kud_all <- NULL
  for (r in seq_len(n_draws)) {
    t_start   <- Sys.time()
    MM_data$t <- deround_draw(r)
    set.seed(if (boot_seed_by_draw) 1000 * r + 3 else 3)   #T_inf is bootstrapped
    kud_all <- rbind(kud_all, cbind(draw = r, methods = meth,
                                    do.call(rbind, lapply(meth, function(m) run_kudrin(MM_data[MM_data$method == m, ])))))
    cat(sprintf("Kudrin tests (%s): de-rounding %d of %d done in %.1f minutes\n", deround_spec$key, r, n_draws,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
  }
  kud_tests        <- setdiff(names(kud_all), c("draw", "methods", "n", "articles", "min", "restricted"))
  empirical_kudrin <- summarize_draws(kud_all, kud_tests)
  empirical_kudrin$min        <- apply(empirical_kudrin[kud_tests], 1, min)   #smallest averaged p-value
  empirical_kudrin$restricted <- as.logical(tapply(kud_all$restricted, kud_all$methods, any)[meth])   #bound in any draw
  saveRDS(kud_all, file = file.path(res_dir, paste0("Kudrin_results_draws", suffix)))
  saveRDS(empirical_kudrin, file = file.path(res_dir, paste0("Kudrin_results", suffix)))
  print(empirical_kudrin)
  
  #How much do the results move across de-roundings? (max - min over draws)
  if (n_draws > 1) {
    cat("\nSpread across de-roundings (max - min):\n")
    for (J0 in J0_list) {
      cat("Projection test, J0 =", J0, "\n")
      print(summarize_draws(proj_all[proj_all$J0 == J0 & proj_all$shift == shift, ], c("pval", "breakdown"), draw_range)[, c("methods", "pval", "breakdown")])
    }
    cat("EKW tests\n")
    print(round(apply(ewk_all, c(1, 2), draw_range), 3))
    cat("Kudrin (2025) tests\n")
    print(summarize_draws(kud_all, kud_tests, draw_range)[, c("methods", kud_tests)])
  }
}

#### Stability across de-rounding draws, as application_cochrane.R records it --------------------------------
# For each option, design and test: the p-value's mean, standard deviation and range across draws, and the share of
# draws that reject at 5%. For B-hat: the share of draws with B-hat > 0 (the projection test rejecting) or, under
# smin and smin_noise, with B-hat above the design's rounding allowance (a certified rejection). Read from the saved
# draws, so it also works when run_tests = FALSE.
stab <- list()
add_stab <- function(option, design, statistic, v, allow = NULL) {
  v <- as.numeric(v); v <- v[is.finite(v)]
  stab[[length(stab) + 1]] <<- data.frame(option = option, design = design, statistic = statistic, draws = length(v),
                                          mean = if (length(v)) mean(v) else NA_real_, sd = if (length(v) > 1) stats::sd(v) else NA_real_,
                                          min = if (length(v)) min(v) else NA_real_, max = if (length(v)) max(v) else NA_real_,
                                          share_rejecting_5pct = if (!length(v)) NA_real_ else if (is.null(allow)) mean(v < 0.05) else mean(v > allow),
                                          stringsAsFactors = FALSE)
}
for (o in options_to_run) {
  key  <- deround_keys[o]
  sfx  <- paste0("_", key, run_tag)
  smin <- key %in% c("smin", "smin_noise")
  f <- file.path(res_dir, paste0("MM_weights_results_draws", sfx))
  if (file.exists(f)) {
    pa <- readRDS(f)
    for (m in meth) {
      for (J0 in J0_list) {
        x <- pa[pa$methods == m & pa$J0 == J0 & pa$shift == shift, ]
        add_stab(key, m, paste0("proj_p_J0_", J0), x$pval)
        add_stab(key, m, paste0("proj_B_J0_", J0), x$breakdown, allow = if (smin) C_round[[as.character(J0)]] * kappa6[[m]] else 0)
      }
      x <- pa[pa$methods == m & pa$J0 == 20 & pa$shift == 0, ]
      if (nrow(x)) { add_stab(key, m, "proj_p_noshift", x$pval); add_stab(key, m, "proj_B_noshift", x$breakdown, allow = 0) }
    }
  }
  f <- file.path(res_dir, paste0("EWK_results_draws", sfx))
  if (file.exists(f)) {
    ea <- readRDS(f)                                           #tests x methods x draws
    if (length(dim(ea)) == 2) ea <- array(ea, c(dim(ea), 1), dimnames = c(dimnames(ea), list(NULL)))
    for (m in meth) for (tn in setdiff(dimnames(ea)[[1]], "min")) add_stab(key, m, tn, ea[tn, paste0("ewk_", tolower(m)), ])
  }
  f <- file.path(res_dir, paste0("Stouffer_results_draws", sfx))
  if (file.exists(f)) {
    sa <- readRDS(f)                                           #draws x methods
    for (m in meth) add_stab(key, m, "Stouffer", sa[, m])
  }
  f <- file.path(res_dir, paste0("Kudrin_results_draws", sfx))
  if (file.exists(f)) {
    ka <- readRDS(f)
    for (m in meth) for (tn in intersect(c("T_inf", "T_H", "T_H_hybrid"), names(ka))) add_stab(key, m, tn, ka[[tn]][ka$methods == m])
  }
}
stability_mm <- do.call(rbind, stab)
if (!is.null(stability_mm)) {
  write.csv(stability_mm, file.path(res_dir, paste0("Stability_application", run_tag, ".csv")), row.names = FALSE)
  u <- stability_mm[stability_mm$draws > 1 & !grepl("^proj_B", stability_mm$statistic), ]
  if (nrow(u)) {
    cat("\nStandard deviation of each test's p-value across the de-rounding draws (options with noise):\n")
    print(round(tapply(u$sd, list(paste(u$option, u$design), u$statistic), identity), 3))
  }
}

#### 6. Design effects of the projection test's features, by design, under smin (deterministic) ---------
hf <- function(x, J) {                      # psi_j(x) * phi(x), j = 0..J
  G <- matrix(0, length(x), J + 1)
  G[, 1] <- dnorm(x); G[, 2] <- x * G[, 1]
  for (j in 1:(J - 1)) G[, j + 2] <- (x * G[, j + 1] - sqrt(j) * G[, j]) / sqrt(j + 1)
  G
}
xi <- function(t, J = 40) 0.5 * (hf(t - 1.96, J) + hf(-t - 1.96, J))   #the test's features (J0 = 20, shift 1.96)
w  <- pmin(1, 2^(20 - (0:40)))              # J0 = 20 weights
t_under <- function(option, r = 1) {        #t-statistics under an option (draw r), setting what deround_draw() reads
  deround_spec <<- list(key = deround_keys[option], label = deround_labels[option], noise = deround_noise[option])
  add_noise    <<- deround_noise[option]
  deround_draw(r)
}
deff <- function(d) {
  X <- sweep(xi(d$t), 2, w, `*`)
  R <- sweep(X, 2, colMeans(X))
  S <- rowsum(R, d$title)                   # paper-level sums
  c(trace = sum(S^2) / sum(R^2), j0 = sum(S[, 1]^2) / sum(R[, 1]^2))
}
MM_data$t      <- t_under(3)
deff_by_design <- sapply(split(MM_data, MM_data$method), deff)[, meth]
cat("\nDesign effects (smin):\n"); print(round(deff_by_design, 2))
saveRDS(deff_by_design, file.path(res_dir, paste0("Deff", run_tag)))

#### 7. Data behind the figures (R1 17b); application_MM_output.R draws them -----------------------------
# Both use the projection test's features above and its honest set: point masses on its grid plus 9999
# (Lemma 3, eq. 34), fitted by nnls with a heavy sum-to-one row. The closest honest t-curve is that mixture.
bx_coeffs      <- function(x, J = 40) t(hf(x / sqrt(2), J) / sqrt(2)) * 2^(-(0:J) / 2)   #c_j of a point mass at x
x_grid         <- seq(-L, L, length.out = numgrid)
B_grid         <- cbind(bx_coeffs(x_grid), 0)                                  #the test's grid, plus 9999
honest_weights <- function(theta) nnls::nnls(rbind(w * B_grid, 1e3), c(w * theta, 1e3))$x

# 7a. Residual figure: R(|t|) = sum_j r_j psi_j(|t| - 1.96) sqrt(phi(|t| - 1.96)), r = w * (theta-hat - c-hat),
#     the part of the reported t-curve no honest t-curve explains, as the test weights it; its square
#     integrates to the squared distance. 90% pointwise bands from an article-level bootstrap.
resid_key <- deround_keys[fig_resid_option]
tt        <- t_under(fig_resid_option)
u_grid    <- seq(0, 5, by = 0.01)
W_u       <- hf(u_grid - 1.96, 40) / sqrt(dnorm(u_grid - 1.96))              #psi_j(s) sqrt(phi(s)), s = |t| - 1.96
resid_curve <- function(theta) {
  r <- w * (theta - drop(B_grid %*% honest_weights(theta)))
  list(R = drop(W_u %*% r), d = sqrt(sum(r^2)))
}
set.seed(fig_seed)
fig_res <- lapply(setNames(meth, meth), function(m) {
  s    <- MM_data$method == m
  X    <- xi(tt[s])                                                  #features of each t-statistic
  fit  <- resid_curve(colMeans(X))
  band <- NULL
  if (fig_resid_boots > 0) {
    rows <- split(seq_len(sum(s)), MM_data$title[s])                 #article-level bootstrap
    Rb   <- vapply(seq_len(fig_resid_boots), function(b) {
      idx <- unlist(rows[sample.int(length(rows), length(rows), replace = TRUE)], use.names = FALSE)
      resid_curve(colMeans(X[idx, , drop = FALSE]))$R
    }, numeric(length(u_grid)))
    band <- apply(Rb, 1, quantile, probs = c(0.05, 0.95), names = FALSE)   #2 x length(u_grid)
  }
  list(R = fit$R, d = fit$d, band = band, n = sum(s), articles = length(unique(MM_data$title[s])))
})
saveRDS(list(option = fig_resid_option, key = resid_key, u = u_grid, boots = fig_resid_boots, results = fig_res),
        file.path(res_dir, paste0("Figdata_residual_tcurve_", resid_key, run_tag)))

# 7b. Reported-versus-honest figure: histogram of reported |t|, the closest honest t-curve, and a curvature
#     check. Every honest t-curve is a mixture of unit-variance normals, so f'' + f >= 0; a Gaussian kernel
#     estimate (bandwidth h) estimates an honest curve with noise variance 1 + h^2, so f'' + f / (1 + h^2) >= 0.
#     Stretches where an article-level bootstrap upper bound on it is below zero are too peaked for honest
#     reporting: simultaneously over fig_window ("sim") and pointwise ("pw").
shape_key <- deround_keys[fig_shape_option]
tt        <- t_under(fig_shape_option)
breaks    <- seq(0, fig_max_t, by = fig_binwidth)
u_c       <- seq(0, fig_max_t, by = 0.02)
u_h       <- seq(0, fig_max_t, by = 0.01)
win       <- u_c >= fig_window[1] & u_c <= fig_window[2]
contrib <- function(a) {                                             #per t-statistic: K'' + K / (1 + h^2), reflected at 0
  z1 <- outer(a, u_c, function(a, u) (u - a) / fig_bandwidth)
  z2 <- outer(a, u_c, function(a, u) (u + a) / fig_bandwidth)
  ((z1^2 - 1) * dnorm(z1) + (z2^2 - 1) * dnorm(z2)) / fig_bandwidth^3 +
    (dnorm(z1) + dnorm(z2)) / (fig_bandwidth * (1 + fig_bandwidth^2))
}
stretches <- function(mask) {                                        #runs of TRUE as c(from, to) rows
  if (!any(mask)) return(matrix(numeric(0), 0, 2))
  k <- which(mask); brk <- c(0, which(diff(k) > 1), length(k))
  t(vapply(seq_len(length(brk) - 1), function(i) c(u_c[k[brk[i] + 1]], u_c[k[brk[i + 1]]]), numeric(2)))
}
honest_density <- function(td) {                                     #density of |t| under the closest honest t-curve
  a <- honest_weights(colMeans(xi(td)))[seq_along(x_grid)]           #(the point 9999 adds no mass here)
  drop((dnorm(outer(u_h - 1.96, x_grid, "-")) + dnorm(outer(-u_h - 1.96, x_grid, "-"))) %*% a)
}
set.seed(fig_shape_seed)
fig_shape <- lapply(setNames(meth, meth), function(m) {
  s      <- MM_data$method == m
  td     <- tt[s]
  counts <- tabulate(findInterval(abs(td), breaks), nbins = length(breaks))[seq_len(length(breaks) - 1)]
  S  <- rowsum(contrib(abs(td)), MM_data$title[s])                   #articles x grid
  N  <- drop(rowsum(rep(1, length(td)), MM_data$title[s]))
  D  <- colSums(S) / sum(N)
  G  <- nrow(S)
  C  <- t(vapply(seq_len(fig_curv_boots), function(b) tabulate(sample.int(G, G, replace = TRUE), nbins = G), numeric(G)))
  Ds <- (C %*% S) / drop(C %*% N)                                    #boots x grid
  sdv <- pmax(apply(Ds, 2, sd), 1e-12)
  dev <- sweep(-sweep(Ds, 2, D), 2, sdv, "/")                         #(D - D*) / sd
  c_w <- quantile(apply(dev[, win, drop = FALSE], 1, max), fig_level, names = FALSE)
  list(density = counts / (length(td) * fig_binwidth), D = D,
       sim = stretches(win & D + c_w * sdv < 0), pw = stretches(apply(Ds, 2, quantile, probs = fig_level, names = FALSE) < 0),
       honest = honest_density(td), n = sum(s), articles = length(unique(MM_data$title[s])))
})
saveRDS(list(option = fig_shape_option, key = shape_key, breaks = breaks, u_c = u_c, u_h = u_h,
             bandwidth = fig_bandwidth, window = fig_window, level = fig_level, boots = fig_curv_boots, results = fig_shape),
        file.path(res_dir, paste0("Figdata_tcurve_shape_", shape_key, run_tag)))

#### 8. Run times (R5 3iv): one call of each test per design on one core, under smin ---------------------
if (time_tests) {
  MM_data$t <- t_under(3)
  tsec <- function(f) {
    t0 <- Sys.time()
    invisible(f())
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  set.seed(1)
  run_times <- t(vapply(meth, function(m) {
    d <- MM_data[MM_data$method == m, ]
    c(n          = nrow(d),
      projection = tsec(function() run_test_weight_fast(d, numcoeffs = coeffs, epsilon_U_pre = epsilon_U_pre, sigma_Y = 1,
                                                        shift = shift, L = L, numgrid = numgrid, boots = boots, J0 = 20, p = p)),
      kudrin     = tsec(function() run_kudrin(d)),
      ekw        = tsec(function() run_ewk(2 * pnorm(-abs(d$t)), ids = as.integer(factor(d$title)))),
      stouffer   = tsec(function() stouffer_p(d$t)))
  }, numeric(5)))
  print(round(run_times, 2))
  saveRDS(run_times, file.path(res_dir, paste0("Run_times", run_tag)))
}

#### 9. Done: the settings again, now with everything the run computed -----------------------------------
save_settings("complete")
cat("\nAll results and settings saved in ", res_dir, " with run_tag '", run_tag, "'. ",
    "Now run application_MM_output.R with run_tag <- \"", run_tag, "\".\n", sep = "")