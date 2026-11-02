# 1. Clear all objects from the global environment
rm(list = ls())

# 2. Force R to release the freed memory to your OS
gc(reset = TRUE)

#load in functions
root               <-  "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_all_weighted.R"))

#Load in the data
data_location <- "C:/Users/sfaridani6/Documents/Data/120246-V1"
setwd(data_location)
MM_raw <- haven::read_dta('MM data.dta')
#MM_raw <- haven::read_dta('MM data fromkudrin.dta')
MM_data <- MM_raw

#Kudrin's (2025) drop of coefficients printed as zero, and a switch for it.
#  TRUE:  as before; every saved file keeps its name.
#  FALSE: keep them (a coefficient printed as zero gives t = 0 under smin, and s_min-scaled noise under smin + noise),
#         and every file this script and make_MM_application_rounding_table_withstouffer.R save gets the tag
#         "_keepmu0", so no earlier result is overwritten. Coefficients that are missing are dropped either way.
drop_zero_mu  <- FALSE
run_tag       <- if (drop_zero_mu) "" else "_keepmu0"
drop          <- is.na(MM_raw$mu) | (drop_zero_mu & MM_raw$mu %in% 0)
cat(sprintf("Coefficients printed as zero: %d, %s.\n", sum(MM_raw$mu %in% 0),
            if (drop_zero_mu) "dropped, following Kudrin (2025)" else "kept; saved files are tagged _keepmu0"))
MM_dropped    <- MM_raw[!drop, ]
MM_data       <- MM_dropped

#### Prep data: de-round the t-statistics (Online Appendix D, Propositions D.1-D.2)
# Four de-rounding options; the script now runs every option in options_to_run, one after another. For a
# t-statistic computed from a coefficient m and a standard error s printed to multiples of a and b, with
# s_min = s - b/2 and s_max = s + b/2 the smallest and largest standard errors consistent with s:
#   1 "none":       t = m / s
#   2 "uniform":    t = (m + V) / W,  V ~ U(-a/2, a/2), W ~ U(s_min, s_max)  (the usual de-rounding)
#   3 "smin":       t = m / s_min                                              (Appendix D.1)
#   4 "smin_noise": t = (m + V) / s_min
# A reported t-statistic is used as printed under options 1 and 3, and gets noise uniform on +-1/2 of its
# last digit under options 2 and 4. Options 2 and 4 are run on n_draws_noise independent draws, and p-values
# and breakdown statistics are averaged over them. Not covered: t-statistics converted from rounded p-values.
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
#     the true unit: Proposition D.1 allows any independent noise.

options_to_run <- 1:4   #1 = none, 2 = uniform (standard), 3 = smallest SE, 4 = smallest SE + noise in the numerator
kappa_cap      <- 0.5   #drop coef/SE t-statistics with kappa = a/s_min above this (the same under every option)
se_rounded_up  <- FALSE #TRUE if papers may have rounded standard errors up: then s_min = sd - b and s_max = sd
n_draws_noise  <- 200   #independent de-roundings for options 2 and 4; try 2 first to check the run time

#Additions for the referees (each can be switched off)
run_stouffer   <- TRUE  #p-curve test for left skew by Stouffer's method, for Table 3 (Simonsohn, point 5)
run_noshift    <- TRUE  #projection test without the 1.96 shift (J0 = 20), for Table 3 (R1 17a)
noshift_draws  <- 200   #de-roundings for the unshifted test under options 2 and 4 (at most n_draws_noise; fewer is faster)
fig_option     <- 3     #de-rounding option of the residual t-curve figure (R1 17b); 3 = smin, which is deterministic
fig_boots      <- 200   #article-level bootstrap draws for the figure's 90% bands (0: no bands)
time_tests     <- TRUE  #time one call of each test per design, for a run-time sentence (R5 3iv)

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
saveRDS(drop_comparison, file = file.path(root, "results", paste0("Drop_comparison_KP", run_tag)))

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
set.seed(1)
U_num <- matrix(runif(nrow(MM_data) * n_draws_noise, -0.5, 0.5), ncol = n_draws_noise)   #numerator noise, one column per draw
set.seed(2)
U_den <- matrix(runif(nrow(MM_data) * n_draws_noise, -0.5, 0.5), ncol = n_draws_noise)   #denominator noise (option 2 only)

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
cat("\nRounding allowance under options 3 and 4 (Appendix D.1, step 3), to compare with B-hat:\n")
print(signif(rbind(mean_kappa6 = kappa6, bound_J0_20 = 0.13 * kappa6, bound_J0_30 = 0.31 * kappa6), 3))

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

### Test settings

cv         <- 1.96
sigma_Y    <- 1 
shift      <- cv
L          <- 6.5 #grid support [-L,L] \cup 9999
numgrid    <- 3000 #total number of elements in grid (density)
boots      <- 1000 #bootstrap repetitions
J          <- 40
coeffs     <- J+1
p          <- 2
J0_list    <- c(20, 30) #one (p, B-hat) column group per J0; J0 = 30 is Table 3's robustness block (R5 3iii)
epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U #pre-compute for speed

### EKW and Kudrin (2025) test wrappers

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

#### Table: projection test, EKW and Kudrin (2025) tests (same layout as Table 3 in the paper)
# Reads the results files saved below for one de-rounding option (so after a full run you can tweak and
# re-run this block alone), prints the table and writes results/table_application<suffix>.tex.
# LaTeX needs booktabs, graphicx (shrink-to-fit wrapper) and hyperref or url (for \url in the notes).

make_application_table <- function(J0_values = c(20),            # one (p, B-hat) group per J0, in this order; each must be in J0_list above
                                   n_draws   = 1,                # number of de-roundings the results files average over
                                   spec      = deround_spec,     # the de-rounding option the results files were run with
                                   suffix    = paste0("_", spec$key, run_tag),
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
    paste0("Results for the ", cite[["brodeur"]], " data. ",
           if (drop_zero_mu) paste0("Following ", cite[["kudrin"]], ", we drop tests whose reported coefficient is missing or zero.") else
             paste0("We drop tests whose reported coefficient is missing; unlike ", cite[["kudrin"]], ", we keep those reported as zero.")),
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
           paste0("\\caption{Empirical Application Results: ", spec$label,
                  if (drop_zero_mu) "" else ", keeping coefficients reported as zero", "}"),
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

#### Run every test under every de-rounding option
tables <- list()
for (deround_option in options_to_run) {
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
    set.seed(1000 * r + 1) #bootstrap seeds differ across draws, so bootstrap noise averages out too
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
  
  #Save the averages over draws (the files make_application_table() reads) and every draw, labeled by option
  setwd(paste0(root,"/results"))
  saveRDS(proj_all, file = paste0("MM_weights_results_draws", suffix))
  for (J0 in J0_list) {
    name <- paste0("MM_weights_results_",J0,"_J0", suffix)
    saveRDS(summarize_draws(proj_all[proj_all$J0 == J0 & proj_all$shift == shift, ], c("pval", "breakdown")), file = name)
    print(readRDS(name))
  }
  if (run_noshift) {                                 #the unshifted test, for Table 3 (R1 17a)
    name <- paste0("MM_weights_results_20_J0_noshift", suffix)
    saveRDS(summarize_draws(proj_all[proj_all$J0 == 20 & proj_all$shift == 0, ], c("pval", "breakdown")), file = name)
    print(readRDS(name))
  }
  
  #EKW tests, on each de-rounding
  ewk_draws <- vector("list", n_draws)
  for (r in seq_len(n_draws)) {
    MM_data$t <- deround_draw(r)
    set.seed(1000 * r + 2)
    ewk_draws[[r]] <- sapply(setNames(meth, paste0("ewk_", tolower(meth))), function(m) {   #tests x methods
      s <- MM_data$method == m
      vapply(run_ewk(2 * pnorm(-abs(MM_data$t[s])), ids = as.integer(factor(MM_data$title[s]))), as_num1, numeric(1))
    })
  }
  ewk_all       <- simplify2array(ewk_draws)                 #tests x methods x draws
  empirical_ewk <- apply(ewk_all, c(1, 2), mean)
  ewk_tests     <- setdiff(rownames(empirical_ewk), "min")
  empirical_ewk["min", ] <- apply(empirical_ewk[ewk_tests, , drop = FALSE], 2, min)   #smallest averaged p-value
  setwd(paste0(root,"/results"))
  saveRDS(ewk_all, file = paste0("EWK_results_draws", suffix))
  saveRDS(empirical_ewk, file = paste0("EWK_results", suffix))
  print(round(empirical_ewk, 3))
  
  #Stouffer's p-curve test for left skew, on each de-rounding (Simonsohn's report, point 5)
  if (run_stouffer) {
    st_all <- t(vapply(seq_len(n_draws), function(r) {
      tt <- deround_draw(r)
      vapply(meth, function(m) stouffer_p(tt[MM_data$method == m]), numeric(1))
    }, numeric(length(meth))))                                 #draws x methods
    st_res <- data.frame(methods = meth, pval = colMeans(st_all), stringsAsFactors = FALSE)
    setwd(paste0(root,"/results"))
    saveRDS(st_all, file = paste0("Stouffer_results_draws", suffix))
    saveRDS(st_res, file = paste0("Stouffer_results", suffix))
    print(st_res)
  }
  
  #Kudrin (2025) tests: T_inf, T_H, T_H_hybrid, on each de-rounding
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
      print(summarize_draws(proj_all[proj_all$J0 == J0 & proj_all$shift == shift, ], c("pval", "breakdown"), draw_range)[, c("methods", "pval", "breakdown")])
    }
    cat("EKW tests\n")
    print(round(apply(ewk_all, c(1, 2), draw_range), 3))
    cat("Kudrin (2025) tests\n")
    print(summarize_draws(kud_all, kud_tests, draw_range)[, c("methods", kud_tests)])
  }
  
  tables[[deround_spec$key]] <- make_application_table(J0_values = J0_list, n_draws = n_draws, spec = deround_spec)
}

#### Table 4 with the EKW tests, once all four options have run (set kappa_cap in that script to the value above)
# source(paste0(root, "/code/table_rounding_ekw.R"))

## effective sample size (uses the t-statistics of the last option and draw run above)

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
sapply(split(MM_data, MM_data$method), deff) # use your design variable

#### Where each design's distortion lies (R1 17b): the part of the reported t-curve no honest t-curve explains
# Under de-rounding option fig_option, each design's t-statistics give the projection test's coefficients
# theta-hat (symmetrized, shifted by 1.96, J = 40). Their projection c-hat onto the honest set (point masses on
# the test's grid plus 9999, J0 = 20 weights; Lemma 3, eq. 34) is unique even though the mixture behind it is
# not (Section 4.1), so the figure plots the residual r = w * (theta-hat - c-hat), not a fitted density:
#   R(|t|) = sum_j r_j psi_j(|t| - 1.96) sqrt(phi(|t| - 1.96)),
# the part of the reported t-curve that no honest t-curve can explain, as the test weights it. Its square
# integrates to the squared distance (printed in each panel; B-hat is it minus the critical value and E(L,
# delta)). Above zero: more t-statistics than honest reporting can produce; below: fewer. In simulations,
# thresholding at 1.96 leaves a deficit just below 1.96, an excess just above and some ringing further out;
# choosing the largest of several specifications leaves a deficit at small |t| and an excess around 1.5-2.2.
# Shaded: 90% pointwise bands from an article-level bootstrap of the whole computation. Needs nnls.
if (!requireNamespace("nnls", quietly = TRUE)) stop("The residual figure needs the nnls package: install.packages(\"nnls\")")
fig_key      <- deround_keys[fig_option]
deround_spec <- list(key = fig_key, label = deround_labels[fig_option], noise = deround_noise[fig_option])
add_noise    <- deround_spec$noise
MM_data$t    <- deround_draw(1)                                     #deterministic under smin; the first draw otherwise
bx_coeffs    <- function(x, J = 40) t(hf(x / sqrt(2), J) / sqrt(2)) * 2^(-(0:J) / 2)   #c_j of a point mass at x (Lemma 3)
B_grid       <- cbind(bx_coeffs(seq(-L, L, length.out = numgrid)), 0)               #the test's grid, plus 9999
u_grid       <- seq(0, 5, by = 0.01)                                                #|t|
W_u          <- hf(u_grid - 1.96, 40) / sqrt(dnorm(u_grid - 1.96))                  #psi_j(s) sqrt(phi(s)), s = |t| - 1.96
resid_curve  <- function(theta) {                                   #projection by nnls with a sum-to-one row
  a <- nnls::nnls(rbind(w * B_grid, 1e3), c(w * theta, 1e3))$x
  r <- w * (theta - drop(B_grid %*% a))
  list(R = drop(W_u %*% r), d = sqrt(sum(r^2)))
}
runs_text <- function(mask) {                                       #"1.42-1.71, 2.03-2.39": stretches of |t| where mask holds
  if (!any(mask)) return("")
  k <- which(mask); brk <- c(0, which(diff(k) > 1), length(k))
  paste(vapply(seq_len(length(brk) - 1), function(i) sprintf("%.2f-%.2f", u_grid[k[brk[i] + 1]], u_grid[k[brk[i + 1]]]), ""),
        collapse = ", ")
}
set.seed(20261030)
fig_res <- lapply(setNames(meth, meth), function(m) {
  s    <- MM_data$method == m
  X    <- xi(MM_data$t[s])                                           #features of each t-statistic
  fit  <- resid_curve(colMeans(X))
  band <- NULL
  if (fig_boots > 0) {
    rows <- split(seq_len(sum(s)), MM_data$title[s])                 #article-level bootstrap
    Rb   <- vapply(seq_len(fig_boots), function(b) {
      idx <- unlist(rows[sample.int(length(rows), length(rows), replace = TRUE)], use.names = FALSE)
      resid_curve(colMeans(X[idx, , drop = FALSE]))$R
    }, numeric(length(u_grid)))
    band <- apply(Rb, 1, quantile, probs = c(0.05, 0.95), names = FALSE)   #2 x length(u_grid)
  }
  list(R = fit$R, d = fit$d, band = band, n = sum(s), articles = length(unique(MM_data$title[s])))
})
B_fig <- tryCatch({                                                 #B-hat from the main run, if it exists
  d <- readRDS(file.path(root, "results", paste0("MM_weights_results_20_J0_", fig_key, run_tag)))
  setNames(vapply(d$breakdown, as_num1, numeric(1)), d$methods)
}, error = function(e) NULL)
fig_summary <- do.call(rbind, lapply(meth, function(m) {
  z <- fig_res[[m]]
  data.frame(design = m, distance_x1000 = round(1e3 * z$d, 2),
             peak_excess_at = u_grid[which.max(z$R)], peak_deficit_at = u_grid[which.min(z$R)],
             excess_where  = if (is.null(z$band)) NA else runs_text(z$band[1, ] > 0),
             deficit_where = if (is.null(z$band)) NA else runs_text(z$band[2, ] < 0), stringsAsFactors = FALSE)
}))
cat("\nWhere each design's distortion lies (", deround_labels[fig_option], "); excess_where and deficit_where: ",
    "stretches of |t| whose 90% band excludes zero:\n", sep = "")
print(fig_summary, row.names = FALSE)
fig_file <- file.path(root, "results", paste0("Figure_residual_tcurve_", fig_key, run_tag))
ylim_fig <- 1e3 * range(0, unlist(lapply(fig_res, function(z) c(z$R, z$band))), na.rm = TRUE)
for (ext in c("pdf", "png")) {
  if (ext == "pdf") {
    grDevices::pdf(paste0(fig_file, ".pdf"), width = 6.5, height = 5.2, pointsize = 9)
  } else {
    grDevices::png(paste0(fig_file, ".png"), width = 6.5, height = 5.2, units = "in", res = 300, pointsize = 9)
  }
  par(mfrow = c(2, 2), mar = c(3.4, 3.6, 2.8, 0.8), mgp = c(2.1, 0.6, 0), las = 1)
  for (m in meth) {
    z <- fig_res[[m]]
    plot(u_grid, 1e3 * z$R, type = "n", ylim = ylim_fig, xlab = "|t|", ylab = expression("Unexplained part (" %*% 10^3 * ")"))
    if (!is.null(z$band)) polygon(c(u_grid, rev(u_grid)), 1e3 * c(z$band[1, ], rev(z$band[2, ])), col = "grey85", border = NA)
    abline(h = 0, col = "grey40")
    abline(v = c(1.645, 1.96, 2.576), lty = 3, col = "grey40")
    lines(u_grid, 1e3 * z$R, lwd = 1.6)
    title(main = sprintf("%s (n = %s, %d articles)", m, format(z$n, big.mark = ","), as.integer(z$articles)),
          line = 1.3, cex.main = 1)
    mtext(sprintf("distance %.1f x 10^-3%s", 1e3 * z$d,
                  if (!is.null(B_fig) && m %in% names(B_fig)) sprintf("; B-hat %.1f x 10^-3", 1e3 * B_fig[[m]]) else ""),
          side = 3, line = 0.3, cex = 0.75)
  }
  grDevices::dev.off()
}
saveRDS(list(u = u_grid, results = fig_res, summary = fig_summary, option = fig_key), paste0(fig_file, "_data"))
cat("Figure written to", paste0(fig_file, ".pdf/.png"), "\n")

#### Run times (R5 3iv): one call of each test per design on one core, on the t-statistics above
if (time_tests) {
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
  saveRDS(run_times, file.path(root, "results", paste0("Run_times", run_tag)))
  sec_range <- function(x) if (diff(range(x)) < 1) sprintf("%.1f", max(x)) else sprintf("%.0f to %.0f", min(x), max(x))
  cat(sprintf(paste0("\nRun-time sentence: On one core, one run of the projection test with %s bootstrap repetitions ",
                     "takes %s seconds per design (n = %s to %s), compared with %s seconds for the three tests of ",
                     "Kudrin (2025) with the same number of bootstrap repetitions and %s seconds for the six EKW tests.\n"),
              format(boots, big.mark = ","), sec_range(run_times[, "projection"]), format(min(run_times[, "n"]), big.mark = ","),
              format(max(run_times[, "n"]), big.mark = ","), sec_range(run_times[, "kudrin"]), sec_range(run_times[, "ekw"])))
}

source(paste0(root,"/code/make_MM_application_rounding_table_withstouffer.R"))
#### =================================================================================================
#### Figure for R1 17(b): the reported t-curve, shaded where its shape rules out honest reporting, by design
#### Paste this block at the end of application_MM_all_four_rounding_withstouffer.R and run it on its own:
#### it runs no tests (the curvature check below is a kernel estimate and a bootstrap, a few minutes). It
#### reads the projection test's p-value and B-hat from the full run's results files.
#### If this R session does not hold the script's de-rounded data (MM_data, deround_draw), it first runs
#### the script's data preparation, which stops before the tests (about a minute) and keeps the other
#### objects in your session. The optional closest honest t-curve needs the nnls package.
####
#### The closest honest t-curve is the mixture of normal t-curves, centered on the projection test's grid,
#### whose coefficients are closest to the data's in the test's metric (J = 40, J0 = 20, shift 1.96). Which
#### distribution of effects produces it is not identified, but the t-curve is: honest t-curves are smooth
#### and their coefficients decay like 2^(-j/2), so the 41 matched coefficients pin it down. The test weights
#### |t| near 1.96 most, so the curve need not follow the bars far from it.
#### =================================================================================================

# 1. The script's data preparation, if this session does not have it
fig_prepare_data <- function(app_file = NULL) {
  GE <- globalenv()
  ready <- exists("deround_draw", envir = GE, inherits = FALSE) && exists("MM_data", envir = GE, inherits = FALSE) &&
    all(c("dr_mu", "dr_smin", "dr_a", "t_source") %in% names(get("MM_data", envir = GE)))
  if (ready) return(invisible(FALSE))
  if (is.null(app_file)) {
    app_file <- if (exists("this_script", envir = GE)) get("this_script", envir = GE) else
      tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
  }
  if (!nzchar(app_file) || !file.exists(app_file))
    stop("Cannot find the application script: run fig_prepare_data(app_file = \"<its path>\")", call. = FALSE)
  message("Running the data preparation of ", basename(app_file), " (it stops before the tests) ...")
  for (e in parse(app_file, keep.source = FALSE)) {
    code <- paste(deparse(e), collapse = " ")
    if (grepl("^rm\\(list = ls\\(\\)\\)", code)) next             #keep the rest of this session
    tryCatch(eval(e, envir = GE), error = function(err) {
      if (grepl("rstudioapi", code)) return(NULL)                 #the script's own path, not needed here
      stop("While preparing the data: ", conditionMessage(err), call. = FALSE)
    })
    if (exists("deround_draw", envir = GE, inherits = FALSE)) return(invisible(TRUE))
  }
  stop(basename(app_file), " never defines deround_draw", call. = FALSE)
}
fig_prepare_data()

# 2. The figure: reported |t|, shaded where its shape rules out honest reporting
# Every honest t-curve is a mixture of unit-variance normals, so its density f satisfies f'' + f >= 0: it cannot
# bend down faster than its own height (phi''(z) = (z^2 - 1) phi(z)). A Gaussian kernel estimate with bandwidth
# h estimates an honest curve with noise variance 1 + h^2, which satisfies the tighter f'' + f / (1 + h^2) >= 0,
# so smoothing cannot create violations. Shaded: stretches of |t| where an article-level bootstrap upper bound
# on f'' + f / (1 + h^2) is below zero, i.e. where the reported t-curve is more sharply peaked than any
# distribution of true effects allows. Dark: 95% simultaneously over the window (prespecified around 1.96);
# light: 95% pointwise. In simulations clustered like these data (200 articles of 25, within-article correlation
# 0.8, h = 0.25), the dark shading flagged threshold p-hacking by every study 35% of the time and honest
# reporting never; the light shading flagged them 70% and 8% of the time. So no shading is not evidence of
# honest reporting: the projection test, which uses the whole t-curve, is the test.
#   option:     de-rounding option, as in the script (1 none, 2 uniform, 3 smin, 4 smin + noise); 4, its first
#               draw, because its noise smooths the rounding lattice of smin t-statistics
#   bandwidth:  of the Gaussian kernel (0.35 is a robustness check); window: the dark shading's range of |t|
#   boots:      article-level bootstrap draws; show_pointwise: light shading; show_honest: also draw the closest
#               honest t-curve (needs nnls), which invites reading its gaps as results added or removed
fig_tcurve_vs_honest <- function(option = 4, binwidth = 0.1, max_t = 5, bandwidth = 0.25, window = c(1.6, 2.6),
                                 boots = 500, level = 0.95, show_pointwise = FALSE, show_honest = TRUE,
                                 results_dir = file.path(root, "results"), seed = 20261102) {
  if (show_honest && !requireNamespace("nnls", quietly = TRUE)) stop("show_honest needs the nnls package: install.packages(\"nnls\")")
  GE    <- globalenv()
  tag   <- if (exists("run_tag", envir = GE)) get("run_tag", envir = GE) else ""   #"_keepmu0" without the zero drop
  keys  <- c("none", "uniform", "smin", "smin_noise")
  labs  <- c("No de-rounding", "Uniform de-rounding", "Smallest consistent standard error",
             "Smallest consistent standard error with uniform noise")
  key   <- keys[option]
  noisy <- option %in% c(2, 4)
  meth  <- c("RCT", "IV", "DID", "RDD")
  
  # t-statistics under this option (its first draw), from the script's deround_draw(); the two globals it
  # reads are put back on exit
  old <- mget(c("deround_spec", "add_noise"), envir = GE, ifnotfound = list(NULL))
  on.exit(for (nm in names(old)) {
    if (is.null(old[[nm]])) {
      if (exists(nm, envir = GE, inherits = FALSE)) rm(list = nm, envir = GE)
    } else assign(nm, old[[nm]], envir = GE)
  }, add = TRUE)
  assign("deround_spec", list(key = key, label = labs[option], noise = noisy), envir = GE)
  assign("add_noise", noisy, envir = GE)
  dat <- get("MM_data", envir = GE)
  tt  <- get("deround_draw", envir = GE)(1)
  
  # The curvature check: per t-statistic, its contribution to f'' + f / (1 + h^2) of the density of |t|
  # (Gaussian kernel, reflected at 0), summed by article for the bootstrap
  u_c <- seq(0, max_t, by = 0.02)
  win <- u_c >= window[1] & u_c <= window[2]
  contrib <- function(a) {
    z1 <- outer(a, u_c, function(a, u) (u - a) / bandwidth)
    z2 <- outer(a, u_c, function(a, u) (u + a) / bandwidth)
    ((z1^2 - 1) * dnorm(z1) + (z2^2 - 1) * dnorm(z2)) / bandwidth^3 +
      (dnorm(z1) + dnorm(z2)) / (bandwidth * (1 + bandwidth^2))
  }
  stretches <- function(mask) {                                      #runs of TRUE as c(from, to) rows
    if (!any(mask)) return(matrix(numeric(0), 0, 2))
    k <- which(mask); brk <- c(0, which(diff(k) > 1), length(k))
    t(vapply(seq_len(length(brk) - 1), function(i) c(u_c[k[brk[i] + 1]], u_c[k[brk[i + 1]]]), numeric(2)))
  }
  as_text <- function(m) if (nrow(m) == 0) "none" else paste(sprintf("%.2f-%.2f", m[, 1], m[, 2]), collapse = ", ")
  
  # The closest honest t-curve, only if asked for (Lemma 3, eq. 34; J = 40, J0 = 20, shift 1.96)
  if (show_honest) {
    Lg <- if (exists("L", envir = GE)) get("L", envir = GE) else 6.5
    ng <- if (exists("numgrid", envir = GE)) get("numgrid", envir = GE) else 3000
    hf <- function(x, J) {
      M <- matrix(0, length(x), J + 1)
      M[, 1] <- dnorm(x)
      M[, 2] <- x * M[, 1]
      for (j in 1:(J - 1)) M[, j + 2] <- (x * M[, j + 1] - sqrt(j) * M[, j]) / sqrt(j + 1)
      M
    }
    w      <- pmin(1, 2^(20 - (0:40)))
    x_grid <- seq(-Lg, Lg, length.out = ng)
    B_grid <- cbind(t(hf(x_grid / sqrt(2), 40) / sqrt(2)) * 2^(-(0:40) / 2), 0)
    u_h    <- seq(0, max_t, by = 0.01)
    honest_density <- function(td) {
      theta <- colMeans(0.5 * (hf(td - 1.96, 40) + hf(-td - 1.96, 40)))
      a     <- nnls::nnls(rbind(w * B_grid, 1e3), c(w * theta, 1e3))$x[seq_along(x_grid)]
      drop((dnorm(outer(u_h - 1.96, x_grid, "-")) + dnorm(outer(-u_h - 1.96, x_grid, "-"))) %*% a)
    }
  }
  
  # By design: histogram, curvature check with its bootstrap bounds, and the closest honest t-curve if asked for
  set.seed(seed)
  breaks <- seq(0, max_t, by = binwidth)
  res <- lapply(setNames(meth, meth), function(m) {
    s      <- dat$method == m
    td     <- tt[s]
    counts <- tabulate(findInterval(abs(td), breaks), nbins = length(breaks))[seq_len(length(breaks) - 1)]
    S  <- rowsum(contrib(abs(td)), dat$title[s])                     #articles x grid
    N  <- drop(rowsum(rep(1, length(td)), dat$title[s]))
    D  <- colSums(S) / sum(N)
    G  <- nrow(S)
    C  <- t(vapply(seq_len(boots), function(b) tabulate(sample.int(G, G, replace = TRUE), nbins = G), numeric(G)))
    Ds <- (C %*% S) / drop(C %*% N)                                  #boots x grid
    sdv <- pmax(apply(Ds, 2, sd), 1e-12)
    dev <- sweep(-sweep(Ds, 2, D), 2, sdv, "/")                       #(D - D*) / sd
    c_w <- quantile(apply(dev[, win, drop = FALSE], 1, max), level, names = FALSE)
    U_sim <- D + c_w * sdv
    U_pw  <- apply(Ds, 2, quantile, probs = level, names = FALSE)
    list(density = counts / (length(td) * binwidth), D = D,
         sim = stretches(win & U_sim < 0), pw = stretches(U_pw < 0),
         honest = if (show_honest) honest_density(td) else NULL,
         n = sum(s), articles = length(unique(dat$title[s])))
  })
  
  # The projection test's p-value and B-hat under this option, from the full run
  tests <- tryCatch({
    r   <- readRDS(file.path(results_dir, paste0("MM_weights_results_20_J0_", key, tag)))
    num <- function(v) vapply(v, function(z) as.numeric(z)[1], numeric(1))
    data.frame(p = num(r$pval), B = num(r$breakdown), row.names = as.character(r$methods))
  }, error = function(e) NULL)
  if (is.null(tests)) message("No results file MM_weights_results_20_J0_", key, tag, ": the panels go without p-values")
  
  # The figure
  dark  <- grDevices::adjustcolor("#D55E00", alpha.f = 0.45)
  light <- grDevices::adjustcolor("#E69F00", alpha.f = 0.25)
  file  <- file.path(results_dir, paste0("Figure_tcurve_shape_", key, tag))
  ymax  <- 1.08 * max(unlist(lapply(res, function(z) c(z$density, z$honest))))
  for (ext in c("pdf", "png")) {
    if (ext == "pdf") {
      grDevices::pdf(paste0(file, ".pdf"), width = 6.5, height = 5.2, pointsize = 9)
    } else {
      grDevices::png(paste0(file, ".png"), width = 6.5, height = 5.2, units = "in", res = 300, pointsize = 9)
    }
    op <- par(mfrow = c(2, 2), mar = c(3.4, 3.6, 2.8, 0.8), mgp = c(2.1, 0.6, 0), las = 1)
    for (m in meth) {
      z <- res[[m]]
      plot(NA, xlim = c(0, max_t), ylim = c(0, ymax), xaxs = "i", yaxs = "i", xlab = "|t|", ylab = "Density")
      if (show_pointwise && nrow(z$pw)) rect(z$pw[, 1] - 0.01, 0, z$pw[, 2] + 0.01, ymax, col = light, border = NA)
      if (nrow(z$sim)) rect(z$sim[, 1] - 0.01, 0, z$sim[, 2] + 0.01, ymax, col = dark, border = NA)
      rect(head(breaks, -1), 0, breaks[-1], z$density, col = "grey75", border = "white", lwd = 0.4)
      if (show_honest) lines(u_h, z$honest, lwd = 1.8)
      abline(v = 1.96, lty = 3, col = "grey30")
      title(main = sprintf("%s (n = %s, %d articles)", m, format(z$n, big.mark = ","), as.integer(z$articles)),
            line = 1.3, cex.main = 1)
      if (!is.null(tests) && m %in% rownames(tests)) {
        mtext(sprintf("projection test p = %s; B-hat %.1f x 10^-3", sub("^0\\.", ".", sprintf("%.2f", tests[m, "p"])),
                      1e3 * tests[m, "B"]), side = 3, line = 0.3, cex = 0.75)
      }
      if (m == meth[1]) {
        keep <- c(TRUE, TRUE, show_pointwise, show_honest)
        
        if(show_pointwise){
          legend("topright", bty = "n", cex = 0.8,
                 legend = c("Reported |t|", sprintf("Too peaked for honest reporting (%g%%, %g-%g)", 100 * level, window[1], window[2]),
                            sprintf("Too peaked for honest reporting (%g%% pointwise)", 100 * level), "Closest honest t-curve")[keep],
                 fill = c("grey75", dark, light, NA)[keep], border = c("white", NA, NA, NA)[keep],
                 lty = c(NA, NA, NA, 1)[keep], lwd = c(NA, NA, NA, 1.8)[keep], col = c(NA, NA, NA, "black")[keep])
          
        }
        else{
          legend("topright", bty = "n", cex = 0.8,
                 legend = c("Reported |t|", "Closest honest t-curve")[keep],
                 fill = c("grey75", NA)[keep], border = c("white",  NA)[keep],
                 lty = c(NA,  1)[keep], lwd = c(NA, 1.8)[keep], col = c(NA,  "black")[keep])
          
        }
        }
    }
    par(op)
    grDevices::dev.off()
  }
  summary <- data.frame(design = meth, n = vapply(res, `[[`, 0, "n"),
                        too_peaked_simultaneous = vapply(res, function(z) as_text(z$sim), ""),
                        too_peaked_pointwise = vapply(res, function(z) as_text(z$pw), ""),
                        min_curvature = round(vapply(res, function(z) min(z$D), 0), 3),
                        at_t = vapply(res, function(z) u_c[which.min(z$D)], 0), row.names = NULL)
  cat("\nStretches of |t| where the reported t-curve is too peaked for honest reporting (", labs[option], ", h = ",
      bandwidth, "; simultaneous over ", window[1], "-", window[2], "):\n", sep = "")
  print(summary, row.names = FALSE)
  write.csv(summary, paste0(file, "_check.csv"), row.names = FALSE)
  saveRDS(list(breaks = breaks, u = u_c, results = res, summary = summary, option = key, bandwidth = bandwidth,
               window = window), paste0(file, "_data"))
  cat("Figure written to ", file, ".pdf/.png\n", sep = "")
  invisible(list(results = res, summary = summary, file = file))
}
fig_out <- fig_tcurve_vs_honest(option = 4)