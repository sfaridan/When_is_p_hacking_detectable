#### =====================================================================================================
#### cochrane_application.R
#### The projection test, Kudrin's (2025) tests, the EKW tests and Stouffer's p-curve test on randomized trials in
#### the Cochrane Database of Systematic Reviews, from Schwab's data set (Open Science Framework, osf.io/xjv9g,
#### doi:10.17605/OSF.IO/XJV9G, CC-BY 4.0), the data behind van Zwet, Schwab and Senn (2021).
#### By default every efficacy outcome of every randomized trial with enough participants, with repeated copies of
#### the same outcome data removed, continuous outcomes printed too coarsely (kappa above kappa_cap) left out, and
#### outcomes clustered by trial; optionally continuous outcomes only. Continuous outcomes are de-rounded as in the
#### application to Brodeur et al.: none, uniform, smallest consistent SE (s_min), and s_min with noise. Writes,
#### each file tagged by the sample's settings:
####   Cochrane_sample_report<tag>.txt   every step that drops outcomes, the final sample and its subsamples, and
####                                     how the sample differs from van Zwet, Schwab and Senn's (2021)
####   Cochrane_results<tag>.csv/.rds    each test's p-value (averaged over de-rounding draws) and its standard deviation
####   Cochrane_stability<tag>.csv       each test's p-value across de-rounding draws: mean, standard deviation, range,
####                                     and the share of draws that reject at 5%
####   Figure_cochrane_tcurve<tag>       reported |z| and the closest honest t-curve
####   Cochrane_sample<tag>.rds          the sample itself, which cochrane_power_experiment.R reads
#### Run first with inspect_only <- TRUE: it writes the sample report, prints the sample, then stops.
#### =====================================================================================================
rm(list = ls())

#### 0. Settings -------------------------------------------------------------------------------------------
root      <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
data_file <- "C:/Users/sfaridani6/Documents/Data/xjv9g-osfstorage-archive/CDSR.RData"   #from osf.io/xjv9g
res_dir   <- file.path(root, "results", "cochrane")
run_label <- ""                     #added to every output file name
inspect_only <- TRUE               #TRUE: write the sample report, print the sample, then stop

# The file's column names, by role. NULL: detect from common names (printed in step 1)
COL <- list(review = NULL,          #review identifier (CD number or DOI)
            study = NULL,           #trial label within a review, e.g. "Smith 2001"
            year = NULL,            #trial year; if absent, parsed from the trial label
            comparison = NULL, outcome = NULL, subgroup = NULL,   #numbers ordering a review's analyses
            measure = NULL,         #effect measure: RR, OR, RD, MD, SMD, HR, ...
            est = NULL, se = NULL,  #effect estimate and its standard error, if the file has them
            n1 = NULL, n2 = NULL,   #participants per arm
            n_total = NULL,         #total participants, used where the arms' sizes are missing
            events1 = NULL, events2 = NULL,                       #events per arm (dichotomous outcomes)
            mean1 = NULL, mean2 = NULL, sd1 = NULL, sd2 = NULL)   #means and SDs per arm (continuous outcomes)

# Sample
rct_only           <- TRUE        #keep only rows the file flags as randomized (its RCT column; table in step 1)
outcome_groups     <- "efficacy"  #keep only these values of outcome.group; NULL keeps all (table in step 1)
continuous_only    <- TRUE       #TRUE: only continuous outcomes (z from means and SDs), in all analyses: leaves out
#dichotomous outcomes and outcomes reported only as an estimate and SE
outcomes_per_trial <- "all"       #"all": every outcome of a trial, repeated copies of the same outcome data removed;
#"first": only each trial's first comparison and outcome (the earlier specification)
n_min              <- 100         #minimum participants per outcome (both arms): small trials have Student-t tails or discrete z
events_min         <- 10          #dichotomous outcomes: minimum events and minimum non-events, summed over both arms
z_from             <- "file"      #"file": the file's own z (van Zwet et al.'s; for dichotomous outcomes, the log-odds-ratio
#z from the counts); "counts": z from the arm-level data here (log OR; mean difference)
use_generic        <- "auto"      #with z_from = "counts": outcomes with only an estimate and SE are used if est/se
#reproduces z from the counts
cluster_by         <- "trial"     #"trial": one cluster per trial label, across reviews; "trial_in_review": one per trial
#within a review; "review": all trials of a review in one cluster (most conservative)
year_split         <- 2005        #compare trials published by this year with later ones
drop_zero_z        <- FALSE       #TRUE: drop outcomes whose z is exactly 0, as Kudrin (2025) drops coefficients printed as
#zero. A diagnostic only: it selects on the estimate and leaves a hole around 0, while
#uniform de-rounding spreads those zeros over their rounding intervals
kappa_cap          <- 0.5         #drop continuous outcomes whose rounding unit in z-units, kappa = 2a / SE_min, exceeds
#this (NULL: no cap), as the application to Brodeur et al. drops t-statistics with kappa
#above its cap. This selects on how results are printed, never on estimates. Dichotomous
#outcomes come from exact counts and outcomes reported as an estimate and SE have no known
#unit, so neither is capped

# De-rounding
treatments <- c("none", "uniform", "smin", "smin_noise")
#for continuous outcomes: "none" the printed z; "uniform" means and SDs drawn uniformly
#within their rounding cells; "smin" z rescaled from the printed SE to SE_min, the
#smallest consistent SE; "smin_noise" the same plus uniform noise at each mean's last
#digit (the same draws as "uniform"). Under the two noisy ones p-values and B-hat are
#averaged over n_draws draws, and their standard deviation across draws is recorded.
#Under "smin" and "smin_noise" a projection-test rejection is certified when B-hat
#exceeds the rounding allowance C mean(kappa^6). Dichotomous outcomes come from exact
#counts and are never de-rounded
n_draws    <- 20

# Tests: as in application_MM_analysis.R
cv <- 1.96; shift <- cv; L <- 6.5; numgrid <- 3000; boots <- 1000; J <- 40; p <- 2
J0_list      <- c(20)
run_noshift  <- FALSE      #the projection test without the 1.96 shift (J0 = 20)
run_kudrin   <- TRUE
run_ekw      <- TRUE
run_stouffer <- TRUE
seed         <- 20261031

# Honest replicas of the binary outcomes: could non-normality alone produce their B-hat?
binary_replicas <- 0          #replicas (e.g. 20; 0: none). Each keeps every binary outcome's arm sizes, draws new event counts
#from binomials at its rates, computes z as for the data, applies the events rule and runs the
#projection test; the binary outcomes' observed B-hat is compared with the replicas'. Needs
#continuous_only <- FALSE
#s_match <- sqrt(max(0, 1 - 1 / var(zb)))  
replica_shrink  <- c(1, 0.86)  #effects in the replicas: 1 = the observed rates; 0.5 = each log odds ratio shrunk halfway
#toward no effect, which concentrates the effects and makes the check harder to pass
replica_events  <- c(10,20,30,40, 50)  #events rules compared: the sample's (events_min) and one where non-normality is negligible
run_main_tests  <- FALSE       #FALSE: stop after the replicas (with continuous_only <- FALSE, a run of the check alone)
stopifnot(outcomes_per_trial %in% c("all", "first"), cluster_by %in% c("trial", "trial_in_review", "review"),
          all(treatments %in% c("none", "uniform", "smin", "smin_noise")), 20 %in% J0_list,
          binary_replicas >= 0, all(replica_shrink >= 0 & replica_shrink <= 1), all(replica_events >= 1))
noisy_tr <- c("uniform", "smin_noise")                  #treatments drawn n_draws times
C_allow  <- c("20" = 0.13, "30" = 0.31)                 #the allowance constant C by J0, for p = 2 (Propositions 3 and 4)

#### 1. Read the file and find its columns --------------------------------------------------------------
read_any <- function(f) {
  if (!file.exists(f)) stop("No file ", f, ": download the data set from osf.io/xjv9g and set data_file")
  ext <- tolower(tools::file_ext(f))
  if (ext %in% c("csv", "tsv", "txt")) {
    if (requireNamespace("data.table", quietly = TRUE)) return(as.data.frame(data.table::fread(f)))
    return(utils::read.csv(f, stringsAsFactors = FALSE, sep = if (ext == "tsv") "\t" else ","))
  }
  if (ext == "rds") return(as.data.frame(readRDS(f)))
  if (ext %in% c("rdata", "rda")) {
    e <- new.env(); nm <- load(f, envir = e)
    if (length(nm) > 1) message("The file holds ", paste(nm, collapse = ", "), "; using ", nm[1])
    return(as.data.frame(get(nm[1], envir = e)))
  }
  if (ext %in% c("xlsx", "xls")) return(as.data.frame(readxl::read_excel(f)))
  if (ext == "dta") return(as.data.frame(haven::read_dta(f)))
  stop("Unknown file type .", ext)
}
raw <- read_any(data_file)
cat(sprintf("Read %s rows and %d columns from %s\nColumns: %s\n", format(nrow(raw), big.mark = ","), ncol(raw),
            basename(data_file), paste(names(raw), collapse = ", ")))

for (v in intersect(c("RCT", "outcome.group", "outcome.flag", "phase", "comparator", "specialty"), names(raw))) {
  cat("\n", v, " (most common values):\n", sep = "")
  print(head(sort(table(raw[[v]], useNA = "ifany"), decreasing = TRUE), 12))
}

norm  <- function(x) gsub("[^a-z0-9]", "", tolower(x))      #"Study.Name" -> "studyname"
names_like <- list(                                          #common names, after norm()
  review     = c("reviewid", "review", "reviewdoi", "doi", "cdnumber", "cdno", "cd", "reviewno", "id"),
  study      = c("studyname", "study", "studlab", "studyid", "studylabel", "trial", "trialname"),
  year       = c("studyyear", "year", "pubyear", "publicationyear"),
  comparison = c("comparisonnr", "comparisonno", "comparisonnumber", "compnr", "comparisonid"),
  outcome    = c("outcomenr", "outcomeno", "outcomenumber", "outcomeid"),
  subgroup   = c("subgroupnr", "subgroupno", "subgroupnumber", "subgroupid"),
  measure    = c("outcomemeasure", "effectmeasure", "measure", "sm", "esmeasure", "effecttype"),
  est        = c("effect", "estimate", "te", "yi", "effectsize", "est", "effectestimate", "effectes"),
  se         = c("se", "seeffect", "sete", "stderr", "standarderror", "sei", "effectse"),
  n1         = c("total1", "ne", "n1", "nexp", "totalexp", "totalintervention", "ntreat", "nt"),
  n2         = c("total2", "nc", "n2", "nctrl", "totalctrl", "totalcontrol", "ncontrol"),
  n_total    = c("effectn", "ntotal", "totaln", "totalparticipants"),
  events1    = c("events1", "evente", "eventse", "eventsexp", "event1", "eventsintervention"),
  events2    = c("events2", "eventc", "eventsc", "eventsctrl", "event2", "eventscontrol"),
  mean1      = c("mean1", "meane", "meanexp", "meanintervention"),
  mean2      = c("mean2", "meanc", "meanctrl", "meancontrol"),
  sd1        = c("sd1", "sde", "sdexp", "sdintervention"),
  sd2        = c("sd2", "sdc", "sdctrl", "sdcontrol"))
find_col <- function(role) {
  if (!is.null(COL[[role]])) {
    if (!COL[[role]] %in% names(raw)) stop("COL$", role, " = '", COL[[role]], "' is not a column of the file")
    return(COL[[role]])
  }
  hit <- names(raw)[match(names_like[[role]], norm(names(raw)))]
  hit <- hit[!is.na(hit)]
  if (length(hit)) hit[1] else NULL
}
cols <- setNames(lapply(names(names_like), find_col), names(names_like))
cat("\nColumns used (role = column):\n")
print(noquote(vapply(cols, function(x) if (is.null(x)) "-" else x, "")))
has <- function(...) all(!vapply(cols[c(...)], is.null, TRUE))
if (!has("study") || !(has("est", "se") || has("n1", "n2", "events1", "events2") || has("n1", "n2", "mean1", "mean2", "sd1", "sd2")))
  stop("Need the trial label and either (est, se) or arm-level data. Set the missing names in COL; the file's columns are ",
       paste(names(raw), collapse = ", "))
if (!has("n1", "n2") && !has("n_total")) warning("No sample sizes found: the large-trial restriction cannot be applied", call. = FALSE)

#### 2. One z per outcome ---------------------------------------------------------------------------------------
g   <- function(role) if (is.null(cols[[role]])) rep(NA, nrow(raw)) else raw[[cols[[role]]]]
num <- function(v) suppressWarnings(as.numeric(as.character(v)))
col_or_na <- function(v) if (v %in% names(raw)) as.character(raw[[v]]) else rep(NA_character_, nrow(raw))
d <- data.frame(row = seq_len(nrow(raw)),
                review = if (is.null(cols$review)) "one review" else as.character(g("review")),
                study = as.character(g("study")),
                comparison = num(g("comparison")), outcome = num(g("outcome")), subgroup = num(g("subgroup")),
                measure = toupper(as.character(g("measure"))),
                est = num(g("est")), se = num(g("se")),
                n1 = num(g("n1")), n2 = num(g("n2")), e1 = num(g("events1")), e2 = num(g("events2")),
                m1 = num(g("mean1")), m2 = num(g("mean2")), s1 = num(g("sd1")), s2 = num(g("sd2")),
                group = col_or_na("outcome.group"), flag = col_or_na("outcome.flag"), specialty = col_or_na("specialty"),
                comparator = col_or_na("comparator"), phase = col_or_na("phase"),
                stringsAsFactors = FALSE)
d$n <- ifelse(is.finite(d$n1 + d$n2), d$n1 + d$n2, num(g("n_total")))
d$rct <- if ("RCT" %in% names(raw)) tolower(trimws(as.character(raw$RCT))) %in% c("yes", "y", "true", "1", "rct") else NA
if (rct_only && !"RCT" %in% names(raw)) stop("rct_only = TRUE, but the file has no RCT column")
if (rct_only && !any(d$rct)) stop("No row counts as randomized: check how the RCT column is coded (table above)")
if (!is.null(outcome_groups) && !"outcome.group" %in% names(raw)) stop("outcome_groups is set, but the file has no outcome.group column")

# Year: the file's, or the 4-digit year in the trial label
yr <- regexpr("(19|20)[0-9]{2}", d$study)
year_label <- rep(NA_real_, nrow(d)); year_label[yr > 0] <- as.numeric(regmatches(d$study, yr))
d$year <- ifelse(is.finite(num(g("year"))), num(g("year")), year_label)

# z from the arm-level data, on the review's own measure: log RR, log OR (also Peto OR) or RD for dichotomous
# outcomes, and the mean difference for continuous ones (an SMD tests the same null, with nearly the same z)
binary <- is.finite(d$e1) & is.finite(d$e2) & is.finite(d$n1) & is.finite(d$n2) & d$n1 > 0 & d$n2 > 0 &
  d$e1 >= 0 & d$e2 >= 0 & d$e1 <= d$n1 & d$e2 <= d$n2
contin <- !binary & is.finite(d$m1) & is.finite(d$m2) & is.finite(d$s1) & is.finite(d$s2) &
  is.finite(d$n1) & is.finite(d$n2) & d$n1 > 1 & d$n2 > 1
d$z_raw <- NA_real_
if (any(binary)) {
  a <- d$e1[binary]; b <- d$n1[binary] - a; c <- d$e2[binary]; dd <- d$n2[binary] - c
  zero <- a == 0 | b == 0 | c == 0 | dd == 0                    #0.5 added to every cell of a table with a zero
  a <- a + 0.5 * zero; b <- b + 0.5 * zero; c <- c + 0.5 * zero; dd <- dd + 0.5 * zero
  ms <- d$measure[binary]
  is_rd <- grepl("(^|[^A-Z])RD([^A-Z]|$)|RISK DIFF", ms)
  is_rr <- !is_rd & grepl("(^|[^A-Z])RR([^A-Z]|$)|RISK RATIO|RELATIVE RISK", ms)
  p1 <- a / (a + b); p2 <- c / (c + dd)
  zb <- ifelse(is_rd, (p1 - p2) / sqrt(p1 * (1 - p1) / (a + b) + p2 * (1 - p2) / (c + dd)),
               ifelse(is_rr, log(p1 / p2) / sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + dd)),
                      log(a * dd / (b * c)) / sqrt(1 / a + 1 / b + 1 / c + 1 / dd)))
  d$z_raw[binary] <- zb
}
if (any(contin)) {
  i <- which(contin)
  d$z_raw[i] <- (d$m1[i] - d$m2[i]) / sqrt(d$s1[i]^2 / d$n1[i] + d$s2[i]^2 / d$n2[i])
}
d$z_raw[!is.finite(d$z_raw)] <- NA
d$z_est <- ifelse(is.finite(d$est) & is.finite(d$se) & d$se > 0, d$est / d$se, NA)

# Does est/se reproduce the z from the counts? It must, for trials with only an estimate and SE to be usable
both <- is.finite(d$z_raw) & is.finite(d$z_est)
gap  <- if (any(both)) median(abs(abs(d$z_raw[both]) - abs(d$z_est[both]))) else NA
if (any(both)) cat(sprintf("\nRows with z from both the counts and est/se: %s; median gap in |z|: %.3f; gap above 0.1: %.1f%%\n",
                           format(sum(both), big.mark = ","), gap, 100 * mean(abs(abs(d$z_raw[both]) - abs(d$z_est[both])) > 0.1)))
generic_ok <- switch(use_generic, auto = isTRUE(gap < 0.05), yes = TRUE, no = FALSE)
if (z_from == "counts" && use_generic == "auto" && any(is.finite(d$z_est) & !is.finite(d$z_raw)))
  cat(if (generic_ok) "est/se agrees with the counts: trials with only an estimate and SE are used\n" else
    "est/se does not clearly agree with the counts (scale?): trials with only an estimate and SE are left out\n")
d$z      <- ifelse(is.finite(d$z_raw), d$z_raw, if (generic_ok) d$z_est else NA)
if ("z" %in% names(raw)) {                                    #the file's own z (van Zwet et al.'s), as a check
  z_file <- num(raw$z)
  bothf  <- is.finite(d$z_raw) & is.finite(z_file)
  cat(sprintf("Rows with z from the counts and the file's z: %s; median gap in |z|: %.3f; gap above 0.1: %.1f%%\n",
              format(sum(bothf), big.mark = ","), median(abs(abs(d$z_raw[bothf]) - abs(z_file[bothf]))),
              100 * mean(abs(abs(d$z_raw[bothf]) - abs(z_file[bothf])) > 0.1)))
  if (z_from == "file") d$z <- ifelse(is.finite(z_file), z_file, NA)
} else if (z_from == "file") stop("z_from = \"file\", but the file has no z column")
d$source <- ifelse(binary & is.finite(d$z_raw), "dichotomous", ifelse(contin & is.finite(d$z_raw), "continuous",
                                                                      ifelse(is.finite(d$z), "estimate and SE", NA)))

# Rounding units of continuous outcomes' means and SDs, from their decimals as stored. A dropped trailing zero
# only makes a unit too large, so the finer of the two arms' units is used; an exact 0 says nothing about decimals
n_dec <- function(x) {
  out <- rep(NA_real_, length(x)); ok <- is.finite(x)
  s <- sub("0+$", "", formatC(abs(x[ok]), format = "f", digits = 8))
  out[ok] <- nchar(s) - regexpr(".", s, fixed = TRUE)
  out
}
unit_of  <- function(x) { u <- 10^-n_dec(x); u[is.finite(x) & x == 0] <- NA; u }
a_m      <- pmin(unit_of(d$m1), unit_of(d$m2), na.rm = TRUE)
b_s      <- pmin(unit_of(d$s1), unit_of(d$s2), na.rm = TRUE)
d$a_mean <- ifelse(is.finite(a_m), a_m, b_s)
d$b_sd   <- ifelse(is.finite(b_s), b_s, a_m)

# kappa: the rounding unit of a continuous outcome's mean difference in z-units. A difference of two means each
# rounded to a can be off by up to a, so its rounding interval is 2a wide; SE_min, the standard error from the
# smallest SDs consistent with the printed ones (s - b/2), is the smallest consistent standard error. (The
# paper's kappa = a / s_min is the same thing for a coefficient printed to a.) SE_min = 0 gives kappa = Inf
d$se     <- sqrt(d$s1^2 / d$n1 + d$s2^2 / d$n2)                                     #from the printed SDs
d$se_min <- sqrt(pmax(d$s1 - d$b_sd / 2, 0)^2 / d$n1 + pmax(d$s2 - d$b_sd / 2, 0)^2 / d$n2)
d$kappa  <- ifelse(d$source %in% "continuous", 2 * d$a_mean / d$se_min, NA)


#### 3. The sample, step by step, with every drop recorded for the sample report ----------------------------
d$trial    <- gsub("\\s+", " ", tolower(trimws(d$study)))
src_levels <- c("dichotomous", "continuous", "estimate and SE", "no z")
steps <- NULL
note  <- function(step, rule, x) {                             #record the sample after each step
  src <- table(factor(ifelse(is.na(x$source), "no z", x$source), levels = src_levels))
  steps <<- rbind(steps, data.frame(step = step, rule = rule, outcomes = nrow(x), trials = length(unique(x$trial)),
                                    reviews = length(unique(x$review)), dichotomous = as.vector(src[1]),
                                    continuous = as.vector(src[2]), estimate_and_SE = as.vector(src[3]),
                                    no_z = as.vector(src[4]), stringsAsFactors = FALSE))
}
inf_last <- function(x) ifelse(is.finite(x), x, Inf)
by_order <- function(x) x[order(x$review, inf_last(x$comparison), inf_last(x$outcome), inf_last(x$subgroup), x$row), ]
d0 <- d                                                        #every row, for the comparison with van Zwet et al.
note("file", "every row of the file: one trial's outcome in one analysis of one review", d)
if (rct_only) {
  d <- d[d$rct %in% TRUE, ]
  note("randomized", "the file's RCT column is \"yes\"", d)
}
if (!is.null(outcome_groups)) {
  d <- d[d$group %in% outcome_groups, ]
  note("outcome group", paste0("outcome.group is ", paste(outcome_groups, collapse = " or ")), d)
}
d <- d[is.finite(d$z), ]
note("z available", if (z_from == "file") "the file's z column is a number" else "z can be computed from the arm-level data (or est/se)", d)
if (continuous_only) {
  d <- d[d$source == "continuous", ]
  note("continuous only", "continuous outcomes, with z from the arms' means and SDs", d)
}
d <- by_order(d)
if (outcomes_per_trial == "all") {
  key <- paste(d$trial, d$n1, d$n2, d$e1, d$e2, d$m1, d$m2, d$s1, d$s2, signif(d$z, 8), sep = "|")
  d   <- d[!duplicated(key), ]
  note("distinct outcome data", paste("the first copy of each trial's outcome data (same trial label, arm sizes, events or",
                                      "means and SDs, and z): a review's subgroup and sensitivity analyses repeat outcome",
                                      "data, and so does a trial met again in another review"), d)
} else {
  d <- d[!duplicated(paste(d$review, d$trial, sep = "\r")), ]
  note("first outcome", "each trial's first comparison and outcome within a review", d)
  d <- d[!(duplicated(paste(d$trial, d$n, sep = "\r")) & is.finite(d$n)), ]
  note("unique trials", "the first review's copy of a trial met again (same label and participants)", d)
}
if (has("n1", "n2") || has("n_total")) {
  d <- d[is.finite(d$n) & d$n >= n_min, ]
  note("participants", sprintf("at least %d participants in the two arms together", n_min), d)
}
d <- d[d$source != "dichotomous" | pmin(d$e1 + d$e2, d$n - d$e1 - d$e2) >= events_min, ]
note("events", sprintf("dichotomous outcomes: at least %d events and %d non-events in the two arms together", events_min, events_min), d)
cont_pre_cap <- d[d$source %in% "continuous", c("kappa", "z")]  #for the report: continuous outcomes before the cap
if (!is.null(kappa_cap)) {
  d <- d[!(d$source %in% "continuous") | (is.finite(d$kappa) & d$kappa <= kappa_cap), ]
  note("kappa", sprintf("continuous outcomes: kappa = 2a / SE_min, the rounding unit of the mean difference in z-units, at most %s",
                        format(kappa_cap)), d)
}
if (drop_zero_z) {
  d <- d[abs(d$z) > 1e-12, ]
  note("nonzero z", "z is not exactly 0", d)
}
steps$dropped <- c(NA, -diff(steps$outcomes))
flow   <- setNames(steps$outcomes, steps$step)
d$unit <- switch(cluster_by, trial = d$trial, trial_in_review = paste(d$review, d$trial, sep = " | "), review = d$review)

# Subsamples of the results table
groups <- list(all = rep(TRUE, nrow(d)),
               early = is.finite(d$year) & d$year <= year_split,
               late = is.finite(d$year) & d$year > year_split,
               dichotomous = d$source == "dichotomous",
               continuous = d$source == "continuous")
group_label <- c(all = "All outcomes", early = sprintf("Published by %d", year_split),
                 late = sprintf("Published after %d", year_split), dichotomous = "Dichotomous outcomes",
                 continuous = "Continuous outcomes")
group_rule  <- c(all = "every outcome in the final sample",
                 early = sprintf("trials published by %d (study.year, else the trial label's year)", year_split),
                 late = sprintf("trials published after %d", year_split),
                 dichotomous = "dichotomous outcomes: z from event counts",
                 continuous = "continuous outcomes: z from the arms' means and SDs")
too_small   <- vapply(groups, sum, 0) < 200
same_as_all <- vapply(names(groups), function(k) k != "all" && identical(groups[[k]], groups$all), TRUE)
skipped     <- names(groups)[too_small | same_as_all]
skip_reason <- ifelse(too_small[skipped], "fewer than 200 outcomes", "identical to all outcomes")
groups      <- groups[!(too_small | same_as_all)]

#### Sample report, written before any test runs ------------------------------------------------------------
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
tag <- paste0("_n", n_min, "_", outcomes_per_trial, "_", cluster_by, if (!is.null(kappa_cap)) paste0("_kappa", kappa_cap) else "",
              if (drop_zero_z) "_nozero" else "",
              if (continuous_only) "_cont" else "", run_label)
show  <- function(x, ...) {                                  #a printed object as report lines: full width, no blank lines
  op <- options(width = 250); on.exit(options(op))
  out <- capture.output(print(x, ...))
  paste0("    ", out[nzchar(trimws(out))])
}
nfmt  <- function(x) format(x, big.mark = ",")
per_unit <- table(d$unit)

# Our reading of van Zwet, Schwab and Senn's (2021) selection, applied to this file: randomized trials, efficacy
# outcomes, each trial's first comparison and outcome, one per trial label
v  <- by_order(d0[d0$rct %in% TRUE & d0$group %in% "efficacy" & is.finite(d0$z), ])
v  <- v[!duplicated(paste(v$review, v$trial, sep = "\r")), ]
v1 <- v[!duplicated(v$trial), ]
vb <- by_order(d0[d0$rct %in% TRUE & d0$group %in% "efficacy" & is.finite(d0$z), ])   #reading B: each review's first
ma <- paste(vb$review, vb$comparison, vb$outcome, sep = "\r")                           #comparison and outcome only
vb <- vb[ma %in% ma[!duplicated(vb$review)], ]
vb <- vb[!duplicated(vb$trial), ]
vc <- by_order(d0[d0$rct %in% TRUE & d0$group %in% "efficacy" & is.finite(d0$z), ])   #reading C: each comparison's
fo <- ave(inf_last(vc$outcome), paste(vc$review, vc$comparison, sep = "\r"), FUN = min) #first outcome (one key: no empty groups)
vc <- vc[inf_last(vc$outcome) == fo, ]
vc <- vc[!duplicated(vc$trial), ]
cont <- d$source == "continuous"
kap_all  <- cont_pre_cap$kappa                                #continuous outcomes before the cap
kap_over <- if (is.null(kappa_cap)) rep(FALSE, length(kap_all)) else !(is.finite(kap_all) & kap_all <= kappa_cap)
kap      <- d$kappa[cont]                                     #in the final sample

rp <- c(
  "SAMPLE REPORT: tests for selective reporting in the Cochrane Database of Systematic Reviews",
  sprintf("Written %s by cochrane_application.R. Data: %s (%s rows), Schwab's data set at osf.io/xjv9g.",
          format(Sys.time(), "%Y-%m-%d %H:%M"), basename(data_file), nfmt(nrow(raw))),
  "A row of the file is one trial's outcome in one analysis (comparison, outcome, subgroup) of one review.",
  "",
  "1. SETTINGS THAT DEFINE THE SAMPLE",
  sprintf("  rct_only = %s; outcome_groups = %s; continuous_only = %s; outcomes_per_trial = %s", rct_only,
          if (is.null(outcome_groups)) "all" else paste(outcome_groups, collapse = ", "), continuous_only, outcomes_per_trial),
  sprintf("  n_min = %d participants; events_min = %d events and non-events (dichotomous); kappa_cap = %s; drop_zero_z = %s",
          n_min, events_min, if (is.null(kappa_cap)) "none" else format(kappa_cap), drop_zero_z),
  sprintf("  z_from = %s (%s); cluster_by = %s; year_split = %d", z_from,
          if (z_from == "file") "the file's own z column" else "z recomputed from the arm-level data", cluster_by, year_split),
  sprintf("  De-rounding treatments: %s%s", paste(treatments, collapse = ", "),
          if (any(noisy_tr %in% treatments)) sprintf(" (%d draws for %s)", n_draws, paste(intersect(noisy_tr, treatments), collapse = " and ")) else ""),
  if (binary_replicas > 0) sprintf("  Honest binary replicas: %d, with effects shrunk by %s; events rules %s", binary_replicas,
                                   paste(replica_shrink, collapse = " and "), paste(replica_events, collapse = " and ")) else NULL,
  "",
  "2. EVERY STEP, IN ORDER (counts after the step; 'dropped' is outcomes removed by the step)",
  show(transform(steps[, c("step", "dropped", "outcomes", "trials", "reviews", "dichotomous", "continuous", "estimate_and_SE", "no_z")],
                 dropped = ifelse(is.na(dropped), "-", format(dropped)))),
  "  (trials: distinct trial labels; reviews: distinct reviews)",
  "  Rules:",
  sprintf("    %-22s %s", steps$step, steps$rule),
  "",
  "3. THE FINAL SAMPLE",
  sprintf("  %s outcomes from %s trials in %s reviews; %s clusters (%s); outcomes per cluster: mean %.1f, median %d, max %d.",
          nfmt(nrow(d)), nfmt(length(unique(d$trial))), nfmt(length(unique(d$review))), nfmt(length(per_unit)),
          cluster_by, mean(per_unit), as.integer(median(per_unit)), as.integer(max(per_unit))),
  sprintf("  Participants per outcome: median %s (quartiles %s-%s). Share with |z| > 1.96: %.3f. Exact zeros: %s.",
          format(median(d$n, na.rm = TRUE)), format(quantile(d$n, .25, na.rm = TRUE)), format(quantile(d$n, .75, na.rm = TRUE)),
          mean(abs(d$z) > 1.96), nfmt(sum(d$z == 0))),
  "  How z is computed:", show(table(d$source)),
  "  Cochrane data type (outcome.flag):", show(table(d$flag, useNA = "ifany")),
  "  Period:", show(table(ifelse(is.finite(d$year), ifelse(d$year <= year_split, sprintf("by %d", year_split),
                                                           sprintf("after %d", year_split)), "year unknown"))),
  "  Comparator:", show(table(d$comparator, useNA = "ifany")),
  "  Phase:", show(table(d$phase, useNA = "ifany")),
  "  Specialty (10 most common):",
  show(setNames(as.data.frame(head(sort(table(d$specialty), decreasing = TRUE), 10), stringsAsFactors = FALSE),
                c("specialty", "outcomes")), row.names = FALSE),
  "",
  "4. SUBSAMPLES IN THE RESULTS TABLE",
  sprintf("  %-12s %-75s %8s outcomes %7s clusters %6s reviews", names(groups), group_rule[names(groups)],
          nfmt(vapply(groups, sum, 0)), nfmt(vapply(groups, function(s) length(unique(d$unit[s])), 0)),
          nfmt(vapply(groups, function(s) length(unique(d$review[s])), 0))),
  if (length(skipped)) sprintf("  Not run: %s (%s).", skipped, skip_reason) else NULL,
  if (length(treatments) > 1) paste("  De-rounding changes only continuous outcomes; under a noisy treatment, a subsample without",
                                    "them is run once instead of n_draws times.") else NULL,
  "",
  "5. HOW THIS SAMPLE DIFFERS FROM VAN ZWET, SCHWAB AND SENN (2021)",
  "  Their sample: one (b, s) pair per trial, from the trial's first or primary efficacy outcome, randomized trials of any",
  "  size; their paper reports 23,551 trials (23,747 in its arXiv version).",
  "  Three readings of that selection, applied to this file (randomized trials, efficacy outcomes, one per trial label):",
  sprintf("    A. each trial's first comparison and outcome in its review:          %s trials (%s with at least %d participants)",
          nfmt(nrow(v1)), nfmt(sum(is.finite(v1$n) & v1$n >= n_min)), n_min),
  sprintf("    B. only each review's first comparison and outcome (its primary meta-analysis): %s trials (%s with at least %d)",
          nfmt(nrow(vb)), nfmt(sum(is.finite(vb$n) & vb$n >= n_min)), n_min),
  sprintf("    C. each comparison's first outcome (a comparison usually lists its primary outcome first): %s trials (%s with at least %d)",
          nfmt(nrow(vc)), nfmt(sum(is.finite(vc$n) & vc$n >= n_min)), n_min),
  "  Differences in our final sample:",
  if (outcomes_per_trial == "all") "    - every efficacy outcome of a trial, not only the first; outcomes are clustered by trial" else
    "    - one outcome per trial, as theirs",
  sprintf("    - only outcomes with at least %d participants%s", n_min,
          if (!continuous_only) sprintf(", and dichotomous outcomes with at least %d events and %d non-events", events_min, events_min) else ""),
  if (continuous_only) "    - continuous outcomes only" else NULL,
  if (!is.null(kappa_cap)) sprintf("    - continuous outcomes whose rounding unit in z-units, kappa = 2a / SE_min, exceeds %s are dropped",
                                   format(kappa_cap)) else NULL,
  if (drop_zero_z) "    - outcomes with z exactly 0 dropped" else NULL,
  sprintf("    - z: %s", if (z_from == "file") "the file's own z column, as theirs" else "recomputed from the arm-level data"),
  "    - repeated copies of the same outcome data are dropped (same trial label, arm data and z)",
  "",
  "6. ROUNDING OF CONTINUOUS OUTCOMES",
  if (length(kap_all)) c(
    "  kappa = 2a / SE_min. a: the finer of the two means' rounding units; SE_min: the standard error from the smallest SDs",
    "  consistent with the printed ones (s - b/2). A difference of two means each rounded to a can be off by up to a, so its",
    "  rounding interval is 2a wide (the paper's kappa = a / s_min, for a coefficient printed to a).",
    sprintf("  Before the cap: %s continuous outcomes, %s of them with z exactly 0 (equal rounded means). kappa, quantiles 50, 75, 90, 99%%:",
            nfmt(length(kap_all)), nfmt(sum(cont_pre_cap$z == 0))),
    show(round(quantile(kap_all[is.finite(kap_all)], c(.5, .75, .9, .99)), 3)),
    if (!is.null(kappa_cap)) sprintf("  Dropped by the cap of %s (kappa above it, or not finite): %s outcomes (%.1f%%), %s of them with z exactly 0.",
                                     format(kappa_cap), nfmt(sum(kap_over)), 100 * mean(kap_over),
                                     nfmt(sum(kap_over & cont_pre_cap$z == 0))) else "  No cap.",
    if (length(kap)) sprintf("  Final sample: %s continuous outcomes, %s of them with z exactly 0; kappa median %.3f, maximum %.3f; rounding allowance 0.13 x mean(kappa^6) = %.2e (J0 = 20; Propositions 3 and 4).",
                             nfmt(sum(cont)), nfmt(sum(d$z[cont] == 0)), median(kap), max(kap), 0.13 * mean(kap^6)) else
                               "  Final sample: no continuous outcomes.",
    if ("uniform" %in% treatments) paste("  Uniform de-rounding draws each mean within +-1/2 of its unit and each SD from its",
                                         "rounding interval; the change in z is added to the outcome's z.") else NULL,
    if (any(c("smin", "smin_noise") %in% treatments))
      c("  Under smin each continuous outcome's z is rescaled from its printed SE to SE_min; under smin + noise, uniform noise at",
        "  each mean's last digit is added to the numerator (the same draws as uniform). Under both, a projection-test rejection",
        "  is certified when B-hat exceeds the subsample's rounding allowance C mean(kappa^6), C = 0.13 for J0 = 20 and 0.31 for",
        "  J0 = 30; the mean counts dichotomous outcomes as exact (kappa = 0) and does not cover outcomes reported as an estimate",
        "  and SE, whose rounding unit is unknown (they stay as printed under every treatment).") else NULL)
  else "  No continuous outcomes in the sample."
)
report_file <- file.path(res_dir, paste0("Cochrane_sample_report", tag, ".txt"))
writeLines(rp, report_file)
cat(paste(rp, collapse = "\n"), "\n\nSample report written to ", report_file, "\n", sep = "")
sample_file <- file.path(res_dir, paste0("Cochrane_sample", tag, ".rds"))   #the sample itself, for cochrane_power_experiment.R
saveRDS(list(sample = d, steps = steps, saved = format(Sys.time()), data_file = data_file,
             settings = list(rct_only = rct_only, outcome_groups = outcome_groups, continuous_only = continuous_only,
                             outcomes_per_trial = outcomes_per_trial, n_min = n_min, events_min = events_min,
                             kappa_cap = kappa_cap, z_from = z_from, cluster_by = cluster_by, drop_zero_z = drop_zero_z)),
        sample_file)
cat("Sample saved to ", sample_file, "\n", sep = "")
if (inspect_only) {
  cat("\n|z| by source (quantiles 0.5, 0.9, 0.99): a source far above the others has z on the wrong scale\n")
  print(round(t(sapply(split(abs(d$z), d$source), quantile, probs = c(.5, .9, .99))), 2))
  print(head(d[, c("review", "study", "year", "source", "n", "z", "unit")], 10))
  stop("Stopped after inspection (inspect_only = TRUE)", call. = FALSE)
}

#### 4. Tests -------------------------------------------------------------------------------------------------
source(file.path(root, "code", "functions_all_weighted.R"))   #run_test_weight_fast(), compute_epsilons()
epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U
if (run_ekw)    source(file.path(root, "code", "Tests.R"))
if (run_kudrin) source(file.path(root, "code", "kudrin_tests.R"))
first <- function(v) as.numeric(v)[1]
getv  <- function(x, nm) if (nm %in% names(x)) x[[nm]] else NA_real_

proj <- function(z, unit, J0, shift) {                        #p-value and B-hat of the projection test
  set.seed(seed)
  r <- run_test_weight_fast(data.frame(t = z, title = as.character(unit)), numcoeffs = J + 1, epsilon_U_pre = epsilon_U_pre,
                            sigma_Y = 1, shift = shift, L = L, numgrid = numgrid, boots = boots, J0 = J0, p = p)
  c(p = first(r$pval), B = first(r$breakdown))
}
run_ewk <- function(ps, ids = seq_along(ps), pmax = 0.15) {    #as in application_MM_analysis.R
  pmin <- min(ps[ps > 0])
  out <- c(LCM = LCM(ps, pmin, pmax), Fisher = Fisher(ps, pmin, pmax), Disc = Discontinuity_test(ps, .05),
           Binomial = Binomial(ps, .04, .05, "c"), CS1 = CoxShi(ps, ids, pmin, pmax, 15, 1, 0),
           CS2B = CoxShi(ps, ids, pmin, pmax, 15, 2, 1))
  vapply(out, first, 0)
}
stouffer_p <- function(t) {                                    #p-curve test for left skew (as in the analysis script)
  a   <- abs(t)
  sig <- is.finite(a) & 2 * pnorm(-a) < 0.05
  if (!any(sig)) return(NA_real_)
  zq  <- qnorm(log(40) + pnorm(-a[sig], log.p = TRUE), log.p = TRUE)
  pnorm(sum(zq) / sqrt(sum(sig)), lower.tail = FALSE)
}
run_tests <- function(z, unit) {                               #every test on one sample; clustered by unit
  out <- c()
  for (J0 in J0_list) {
    r <- proj(z, unit, J0, shift)
    out[paste0("proj_p_J0_", J0)] <- r[["p"]]; out[paste0("proj_B_J0_", J0)] <- r[["B"]]
  }
  if (run_noshift) { r <- proj(z, unit, 20, 0); out["proj_p_noshift"] <- r[["p"]]; out["proj_B_noshift"] <- r[["B"]] }
  if (run_kudrin) {
    set.seed(seed)
    kt <- kudrin_tests(z, cluster = unit, boots = boots)
    out["T_inf"] <- first(kt$pval[["T_inf"]]); out["T_H"] <- first(kt$pval[["T_H"]])
    out["T_H_hybrid"] <- first(kt$pval[["T_H_hybrid"]]); out["hybrid_restricted"] <- as.numeric(isTRUE(kt$restricted))
  }
  if (run_ekw) {
    set.seed(seed)
    e <- run_ewk(2 * pnorm(-abs(z)), ids = as.integer(factor(unit)))
    out[names(e)] <- e
  }
  if (run_stouffer) out["Stouffer"] <- stouffer_p(z)
  out
}

# Honest replicas of the binary outcomes. Each replica keeps every binary outcome's arm sizes and draws new event counts
# from binomials, arm 2 being the control arm: at the observed rates (shrink = 1), or at the observed control rate with
# the log odds ratio shrunk toward 0 (shrink < 1). z is computed as for the data (Wald log odds ratio, 0.5 added to every
# cell of a table with a zero), the events rule is applied to the replicated counts, and the projection test (J0 = 20)
# is run, clustered by trial. There is no selection in a replica, only non-normality: if the observed B-hat is far above
# the replicas', non-normality cannot be its source. Common random numbers: within a replica every shrink factor uses
# the same uniforms
if (binary_replicas > 0) {
  bin <- which(d$source == "dichotomous")
  if (!length(bin)) {
    message("binary_replicas: the sample has no dichotomous outcomes; set continuous_only <- FALSE. Skipped.")
  } else {
    wald_z <- function(a, n1, c, n2) {
      b <- n1 - a; dd <- n2 - c
      zero <- a == 0 | b == 0 | c == 0 | dd == 0
      a <- a + 0.5 * zero; b <- b + 0.5 * zero; c <- c + 0.5 * zero; dd <- dd + 0.5 * zero
      log(a * dd / (b * c)) / sqrt(1 / a + 1 / b + 1 / c + 1 / dd)
    }
    passes <- function(e1, n1, e2, n2, E) pmin(e1 + e2, n1 + n2 - e1 - e2) >= E
    Bd    <- d[bin, ]
    rate2 <- (Bd$e2 + 0.5) / (Bd$n2 + 1)                                   #control rate, kept off 0 and 1
    beta  <- log((Bd$e1 + 0.5) * (Bd$n2 - Bd$e2 + 0.5) / ((Bd$n1 - Bd$e1 + 0.5) * (Bd$e2 + 0.5)))
    rate1 <- function(s) plogis(qlogis(rate2) + s * beta)
    z_obs <- ifelse(is.finite(Bd$z_raw), Bd$z_raw, Bd$z)                    #the same formula as the replicas'
    obs <- do.call(rbind, lapply(replica_events, function(E) {
      keep <- passes(Bd$e1, Bd$n1, Bd$e2, Bd$n2, E)
      r0   <- proj(z_obs[keep], Bd$unit[keep], 20, shift)
      data.frame(events_min = E, n_obs = sum(keep), p_obs = r0[["p"]], B_obs = r0[["B"]])
    }))
    cat("\nBinary outcomes, observed: B-hat x 10^3 by events rule:", paste(sprintf("%d: %.3f (p = %.3f)", obs$events_min,
                                                                                   1e3 * obs$B_obs, obs$p_obs), collapse = "; "), "\n")
    rep_rows <- NULL
    for (r in seq_len(binary_replicas)) {
      t0 <- Sys.time()
      set.seed(seed + 5000 + r)
      U1  <- runif(nrow(Bd)); U2 <- runif(nrow(Bd))
      e2r <- qbinom(U2, Bd$n2, rate2)
      for (sh in replica_shrink) {
        e1r <- qbinom(U1, Bd$n1, rate1(sh))
        zr  <- wald_z(e1r, Bd$n1, e2r, Bd$n2)
        for (E in replica_events) {
          keep <- passes(e1r, Bd$n1, e2r, Bd$n2, E) & is.finite(zr)
          rr   <- proj(zr[keep], Bd$unit[keep], 20, shift)
          rep_rows <- rbind(rep_rows, data.frame(replica = r, shrink = sh, events_min = E, n = sum(keep),
                                                 p = rr[["p"]], B = rr[["B"]]))
        }
      }
      cat(sprintf("Binary replica %d of %d: B-hat x 10^3 %s (%.1f min)\n", r, binary_replicas,
                  paste(sprintf("%.3f", 1e3 * tail(rep_rows$B, length(replica_shrink) * length(replica_events))), collapse = ", "),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
    replica_summary <- do.call(rbind, lapply(split(rep_rows, list(rep_rows$events_min, rep_rows$shrink), drop = TRUE), function(x) {
      o <- obs[obs$events_min == x$events_min[1], ]
      data.frame(events_min = x$events_min[1], shrink = x$shrink[1], n_obs = o$n_obs, B_obs = o$B_obs, p_obs = o$p_obs,
                 replicas = nrow(x), n_replica = mean(x$n), B_mean = mean(x$B), B_q95 = quantile(x$B, 0.95, names = FALSE),
                 B_max = max(x$B), share_B_at_least_obs = mean(x$B >= o$B_obs), share_rejecting = mean(x$p <= 0.05),
                 stringsAsFactors = FALSE)
    }))
    replica_summary <- replica_summary[order(replica_summary$events_min, -replica_summary$shrink), ]
    rownames(replica_summary) <- NULL
    dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(replica_summary, file.path(res_dir, paste0("Cochrane_binary_replicas", tag, ".csv")), row.names = FALSE)
    saveRDS(list(summary = replica_summary, replicas = rep_rows, observed = obs,
                 settings = list(binary_replicas = binary_replicas, replica_shrink = replica_shrink,
                                 replica_events = replica_events, seed = seed)),
            file.path(res_dir, paste0("Cochrane_binary_replicas", tag, ".rds")))
    cat("\nHonest replicas of the binary outcomes (B-hat x 10^3; observed against the replicas'):\n")
    show_rep <- replica_summary
    for (cn in c("B_obs", "B_mean", "B_q95", "B_max")) show_rep[[cn]] <- round(1e3 * show_rep[[cn]], 3)
    print(show_rep, row.names = FALSE, digits = 3)
  }
}
if (!run_main_tests) stop("run_main_tests = FALSE: stopped after the binary replicas", call. = FALSE)

# The four treatments of rounding, for draw r (continuous outcomes only; the others stay as printed):
#   "none"        the printed z
#   "uniform"     each mean drawn within +-1/2 of its unit and each SD from its rounding interval (kept above 0);
#                 the change this makes to the outcome's z from the counts is added to its z
#   "smin"        z rescaled from the printed SE to SE_min, the smallest consistent SE (Section 5.2)
#   "smin_noise"  the same, plus uniform noise at each mean's last digit in the numerator, a (U1 - U2) / SE_min
# Draw r of "uniform" and of "smin_noise" use the same uniforms, so the two differ only in the denominator. With
# z_from = "file" the file's z keeps its own definition. Without a kappa cap SE_min can be 0; such outcomes stay
# as printed under "smin" and "smin_noise"
cont_i <- which(d$source == "continuous" & is.finite(d$z_raw))
z_treat <- function(tr, r = 1) {
  z <- d$z
  if (tr == "none" || !length(cont_i)) return(z)
  i <- cont_i
  if (tr %in% noisy_tr) {
    set.seed(seed + r)
    U <- matrix(runif(4 * length(i)), ncol = 4)
  }
  if (tr == "uniform") {
    a <- d$a_mean[i]; b <- d$b_sd[i]
    sd_draw <- function(s, u) { lo <- pmax(s - b / 2, 0); lo + (s + b / 2 - lo) * u }
    m1 <- d$m1[i] + a * (U[, 1] - 0.5); m2 <- d$m2[i] + a * (U[, 2] - 0.5)
    s1 <- sd_draw(d$s1[i], U[, 3]);     s2 <- sd_draw(d$s2[i], U[, 4])
    zc <- (m1 - m2) / sqrt(s1^2 / d$n1[i] + s2^2 / d$n2[i])
    z[i] <- ifelse(is.finite(zc), z[i] + (zc - d$z_raw[i]), z[i])
    return(z)
  }
  ok      <- d$se_min[i] > 0
  stretch <- ifelse(ok, d$se[i] / d$se_min[i], 1)
  z[i]    <- z[i] * stretch
  if (tr == "smin_noise") {
    noise <- ifelse(ok, d$a_mean[i] * (U[, 1] - U[, 2]) / d$se_min[i], 0)
    z[i]  <- z[i] + ifelse(is.finite(noise), noise, 0)
  }
  z
}

# The rounding allowance of a subsample s at J0: C mean(kappa^6) over its outcomes, dichotomous ones exact
kap6     <- ifelse(is.finite(d$kappa), d$kappa^6, 0)
allow_of <- function(s, J0) {
  if (p == 2 && as.character(J0) %in% names(C_allow)) C_allow[[as.character(J0)]] * mean(kap6[s]) else NA_real_
}


draw_rows <- list()
for (tr in treatments) {
  for (r in if (tr %in% noisy_tr) seq_len(n_draws) else 1L) {
    z_all <- z_treat(tr, r)
    for (k in names(groups)) {
      s     <- groups[[k]]
      noisy <- tr %in% noisy_tr && any(s[cont_i])
      if (!noisy && r > 1) next                               #nothing to de-round: one run is enough
      t0  <- Sys.time()
      res <- run_tests(z_all[s], d$unit[s])
      mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      draw_rows[[length(draw_rows) + 1]] <- data.frame(treatment = tr, draw = r, group = k, as.list(res), minutes = mins,
                                                       stringsAsFactors = FALSE)
      cat(sprintf("[%s%s] %-22s projection p = %.3f (B-hat %.2f x 10^-3%s)%s | Kudrin T_inf p = %.3f, T_H p = %.3f, hybrid p = %.3f | %.1f min\n",
                  tr, if (noisy) sprintf(", draw %d of %d", r, n_draws) else "", group_label[[k]],
                  getv(res, "proj_p_J0_20"), 1e3 * getv(res, "proj_B_J0_20"),
                  if (tr %in% c("smin", "smin_noise")) sprintf(", allowance %.2f", 1e3 * allow_of(s, 20)) else "",
                  if (30 %in% J0_list) sprintf(", J0 = 30 p = %.3f", getv(res, "proj_p_J0_30")) else "",
                  getv(res, "T_inf"), getv(res, "T_H"), getv(res, "T_H_hybrid"), mins))
    }
    if (tr %in% noisy_tr && r == 1 && n_draws > 1) {
      done <- vapply(draw_rows, function(x) x$treatment == tr, TRUE)
      cat(sprintf("  (about %.0f more minutes for the remaining %d draws)\n",
                  sum(vapply(draw_rows[done], function(x) x$minutes, 0)) * (n_draws - 1), n_draws - 1))
    }
  }
}

# Averages and standard deviations over draws (one draw under "none" and "smin", and for subsamples with nothing to
# de-round); the rounding allowance by J0 and, under "smin" and "smin_noise", whether the averaged B-hat exceeds it
draws_tab <- do.call(rbind, draw_rows)
num_cols  <- setdiff(names(draws_tab), c("treatment", "draw", "group"))
p_cols    <- intersect(c(grep("^proj_p_", num_cols, value = TRUE), "T_inf", "T_H", "T_H_hybrid", "LCM", "Fisher", "Disc",
                         "Binomial", "CS1", "CS2B", "Stouffer"), num_cols)
b_cols    <- grep("^proj_B_", num_cols, value = TRUE)
sd_or_na  <- function(v) if (length(v) > 1) stats::sd(v) else NA_real_
by_cell   <- split(draws_tab, paste(draws_tab$treatment, draws_tab$group))
results <- do.call(rbind, lapply(by_cell, function(x) {
  s <- groups[[x$group[1]]]
  data.frame(treatment = x$treatment[1], group = x$group[1], label = group_label[[x$group[1]]], draws = nrow(x),
             n = sum(s), clusters = length(unique(d$unit[s])), share_sig = mean(abs(d$z[s]) > 1.96),
             as.list(colMeans(x[, num_cols, drop = FALSE])),
             as.list(setNames(vapply(p_cols, function(cn) sd_or_na(x[[cn]]), 0), paste0("sd_", p_cols))),
             as.list(setNames(vapply(J0_list, function(J0) allow_of(s, J0), 0), paste0("allowance_J0_", J0_list))),
             as.list(setNames(vapply(J0_list, function(J0) {
               if (!(x$treatment[1] %in% c("smin", "smin_noise"))) return(NA_real_)
               as.numeric(mean(x[[paste0("proj_B_J0_", J0)]]) > allow_of(s, J0))
             }, 0), paste0("certified_J0_", J0_list))),
             stringsAsFactors = FALSE)
}))
results <- results[order(match(results$treatment, treatments), match(results$group, names(groups))), ]
rownames(results) <- NULL

# Stability across draws: for each test, its p-value's mean, standard deviation and range, and the share of draws that
# reject at 5%. For B-hat: the share of draws with B-hat > 0 (the projection test rejecting) or, under "smin" and
# "smin_noise", with B-hat above the rounding allowance (a certified rejection)
stability <- do.call(rbind, lapply(by_cell, function(x) do.call(rbind, lapply(c(p_cols, b_cols), function(cn) {
  v     <- x[[cn]]
  allow <- if (x$treatment[1] %in% c("smin", "smin_noise") && grepl("^proj_B_J0_", cn)) {
    allow_of(groups[[x$group[1]]], as.numeric(sub("^proj_B_J0_", "", cn)))
  } else 0
  data.frame(treatment = x$treatment[1], group = x$group[1], statistic = cn, draws = length(v), mean = mean(v),
             sd = sd_or_na(v), min = min(v), max = max(v),
             share_rejecting_5pct = if (cn %in% p_cols) mean(v < 0.05) else mean(v > allow), stringsAsFactors = FALSE)
}))))
stability <- stability[order(match(stability$treatment, treatments), match(stability$group, names(groups)),
                             match(stability$statistic, c(p_cols, b_cols))), ]
rownames(stability) <- NULL

#### 5. Save --------------------------------------------------------------------------------------------------------
settings <- list(saved = format(Sys.time()), data_file = data_file, columns = cols, rct_only = rct_only,
                 outcome_groups = outcome_groups, continuous_only = continuous_only, outcomes_per_trial = outcomes_per_trial,
                 n_min = n_min, events_min = events_min, z_from = z_from, use_generic = use_generic,
                 generic_used = generic_ok, est_se_gap = gap, cluster_by = cluster_by, year_split = year_split,
                 drop_zero_z = drop_zero_z, kappa_cap = kappa_cap, treatments = treatments, n_draws = n_draws, shift = shift, L = L,
                 numgrid = numgrid, boots = boots, J = J, p = p, J0_list = J0_list, epsilon_U_pre = epsilon_U_pre,
                 seed = seed, steps = steps, flow = flow)
saveRDS(list(results = results, draws = draws_tab, stability = stability, settings = settings,
             sample = d[, c("review", "study", "year", "source", "n", "z", "unit")]),
        file.path(res_dir, paste0("Cochrane_results", tag, ".rds")))
write.csv(results, file.path(res_dir, paste0("Cochrane_results", tag, ".csv")), row.names = FALSE)
write.csv(stability, file.path(res_dir, paste0("Cochrane_stability", tag, ".csv")), row.names = FALSE)
cat("\nResults (averages over the draws of the noisy treatments):\n"); print(format(results, digits = 3), row.names = FALSE)
for (trn in intersect(noisy_tr, treatments)) {
  u <- stability[stability$treatment == trn & stability$draws > 1 & stability$statistic %in% p_cols, ]
  if (!nrow(u)) next
  tab <- function(col) {
    m <- tapply(u[[col]], list(u$group, u$statistic), identity)
    m[intersect(names(groups), rownames(m)), intersect(p_cols, colnames(m)), drop = FALSE]
  }
  cat(sprintf("\nStandard deviation of each test's p-value across the %s draws:\n", trn)); print(round(tab("sd"), 3))
  cat(sprintf("\nShare of %s draws in which each test rejects at 5%%:\n", trn)); print(round(tab("share_rejecting_5pct"), 2))
}

#### 6. Figure: reported |z| and the closest honest t-curve (no de-rounding), as in the application -------------
hf <- function(x, J) {                                         #psi_j(x) phi(x), j = 0..J
  G <- matrix(0, length(x), J + 1)
  G[, 1] <- dnorm(x); G[, 2] <- x * G[, 1]
  for (j in 1:(J - 1)) G[, j + 2] <- (x * G[, j + 1] - sqrt(j) * G[, j]) / sqrt(j + 1)
  G
}
w      <- pmin(1, 2^(20 - (0:40)))                              #J0 = 20
x_grid <- seq(-L, L, length.out = numgrid)
B_grid <- cbind(t(hf(x_grid / sqrt(2), 40) / sqrt(2)) * 2^(-(0:40) / 2), 0)   #the test's grid, plus 9999
u_h    <- seq(0, 5, by = 0.01)
honest_density <- function(z) {
  theta <- colMeans(0.5 * (hf(z - 1.96, 40) + hf(-z - 1.96, 40)))
  a     <- nnls::nnls(rbind(w * B_grid, 1e3), c(w * theta, 1e3))$x[seq_along(x_grid)]
  drop((dnorm(outer(u_h - 1.96, x_grid, "-")) + dnorm(outer(-u_h - 1.96, x_grid, "-"))) %*% a)
}
res0       <- results[results$treatment == "none", ]
fig_groups <- intersect(c("all", "early", "late"), res0$group)
if (requireNamespace("nnls", quietly = TRUE) && length(fig_groups)) {
  breaks <- seq(0, 5, by = 0.1)
  panels <- lapply(setNames(fig_groups, fig_groups), function(k) {
    z <- abs(d$z[groups[[k]]])
    list(density = tabulate(findInterval(z, breaks), nbins = length(breaks))[seq_len(length(breaks) - 1)] / (length(z) * 0.1),
         honest = honest_density(z), n = length(z))
  })
  ymax <- 1.08 * max(unlist(lapply(panels, function(z) c(z$density, z$honest))))
  file <- file.path(res_dir, paste0("Figure_cochrane_tcurve", tag))
  for (ext in c("pdf", "png")) {
    if (ext == "pdf") grDevices::pdf(paste0(file, ".pdf"), width = 2.4 * length(fig_groups) + 0.6, height = 3.2, pointsize = 9)
    else grDevices::png(paste0(file, ".png"), width = 2.4 * length(fig_groups) + 0.6, height = 3.2, units = "in", res = 300, pointsize = 9)
    op <- par(mfrow = c(1, length(fig_groups)), mar = c(3.4, 3.6, 2.8, 0.8), mgp = c(2.1, 0.6, 0), las = 1)
    for (k in fig_groups) {
      z <- panels[[k]]; r <- res0[res0$group == k, ]
      plot(NA, xlim = c(0, 5), ylim = c(0, ymax), xaxs = "i", yaxs = "i", xlab = "|z|", ylab = "Density")
      rect(head(breaks, -1), 0, breaks[-1], z$density, col = "grey75", border = "white", lwd = 0.4)
      lines(u_h, z$honest, lwd = 1.8)
      abline(v = 1.96, lty = 3, col = "grey30")
      title(main = sprintf("%s (n = %s)", group_label[[k]], format(z$n, big.mark = ",")), line = 1.3, cex.main = 1)
      mtext(sprintf("projection test p = %s; B-hat %.1f x 10^-3", sub("^0\\.", ".", sprintf("%.2f", r$proj_p_J0_20)),
                    1e3 * r$proj_B_J0_20), side = 3, line = 0.3, cex = 0.75)
      if (k == fig_groups[1]) legend("topright", bty = "n", cex = 0.8, legend = c("Reported |z|", "Closest honest t-curve"),
                                     fill = c("grey75", NA), border = c("white", NA), lty = c(NA, 1), lwd = c(NA, 1.8),
                                     col = c(NA, "black"))
    }
    par(op)
    grDevices::dev.off()
  }
  cat("Figure written to ", file, ".pdf/.png\n", sep = "")
}
cat("Results written to ", file.path(res_dir, paste0("Cochrane_results", tag, ".csv/.rds")), "\n", sep = "")