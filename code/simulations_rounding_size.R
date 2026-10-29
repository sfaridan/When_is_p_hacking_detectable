################################################################################
#  simulations_rounding_size.R
#
#  Size of the projection test, Kudrin's (2025) three tests (T_inf, T_H, T_H^hybrid) and the
#  EKW (2022) tests when NOTHING is selectively reported, but the coefficients and standard
#  errors behind the t-statistics are rounded, as in published regression tables.
#  Companion of simulations_var_dial_plus_gammahat.R: same engine (kudrin_figures_engine.R),
#  same tests called with the same arguments, same panels (distributions of true effects Pi0)
#  and the same pool design. What changes is how a reported t-statistic is generated (rounding
#  instead of p-hacking) and the horizontal axis, which is now the severity of rounding.
#
#  Rounding DGP (no selective reporting anywhere)
#    Each study has an honest t-statistic T* = h + e (h ~ Pi0, e ~ Student-t(nu), as in
#    run_sims()) and reports a coefficient b and a standard error se with b / se = T*. A table
#    shows b and se to the same number of decimal places (the convention Kranz and Puetz,
#    2022, use for these data). Measured in units of the last reported digit, se equals x and
#    b equals T* x, so the table shows the integers
#        s = round(x)       significand of the reported se    (0.012 -> 12)
#        m = round(T* x)    the reported coefficient in the same units,
#    and the meta-analyst computes t = m / s from the table. Hence |t - T*| <= (1 + T*) / (2 s):
#    the smaller the significand, the coarser the rounding relative to the noise.
#    Rows of the figures (CFG$rows), the options of application_MM_all_four_rounding.R, in units
#    of the last digit (which b and se share), with s_min = s - 1/2 the smallest se printed as s:
#      "raw"         t = m / s, as printed                                       (option 1)
#      "deround"     t = (m + V) / W, V ~ U(-1/2, 1/2), W ~ U(s_min, s_max)    (2: the usual one)
#      "smin"        t = m / s_min                                     (3: paper, Section 5.2.1)
#      "smin_noise"  t = (m + V) / s_min, with the same V as "deround"          (4)
#    As in the application, every row analyzes the same results: those with kappa = a / s_min
#    <= CFG$kappa_cap (the paper's cap is 0.3: significands of 4 and up). Proposition 4: the true
#    se is at least s_min, so dividing by s_min only stretches T by lambda = se / s_min >= 1,
#    which is absorbed, and leaves the coefficient's rounding, to multiples of kappa.
#
#  Calibration to Brodeur, Cook and Heyes (2020, AER, "Methods matter")
#    x = s_BCH + U(-1/2, 1/2), where s_BCH is drawn from the significands of the standard
#    errors reported in their data, by default given the study's effect: from the results whose
#    reported |t| lies in the same decile bin as |h| (CFG$brodeur$condition_on = "effect"). Then
#    precision depends on the effect but not on the noise Z, as Proposition 4 assumes, and the
#    cap selects on H only, which keeps every row's null. "t" conditions on the realized
#    |t| = |h + Z| instead, so the cap also selects on Z: a stress test outside Proposition 4.
#    Decimal places are recovered as in Kranz and Puetz (2022): the larger of the decimal
#    places of the coefficient and of the standard error, which restores most trailing zeros.
#
#  Severity of rounding (horizontal axis): coarseness c
#    x = x_BCH / c. At coarseness c every rounding unit is c times as large, relative to the
#    standard error, as in Brodeur et al. (2020): c = 0 is no rounding, c = 1 (vertical line) is
#    the rounding in their data, c = 10 is every table showing one decimal place fewer, and
#    c = 1/2 is rounding twice as fine. A standard error that would show as zero is shown with
#    as many more decimals as it takes to be nonzero (CFG$zero_se_rule = "more_decimals"), or
#    the result is left out ("drop": a zero se gives no t-statistic). At c = 0 all rows report
#    the honest T*, consume the same random numbers and therefore agree exactly.
#
#  Uses your code only (all in CODE_DIR):
#    kudrin_tests.R                            kudrin_tests()           T_inf, T_H, T_H^hybrid
#    functions_all_weighted.R                  run_test_weight_fast()   projection test
#    reproduce_kudrin_simulations_functions.R  sourced as in the companion script
#    kudrin_figures_engine.R                   test calls, cells, resuming, parallel runs
#    Tests.R (EKW)                             only if CFG$run_ekw = TRUE
#  and the Brodeur et al. (2020) data, "MM data.dta" from their replication package
#  (openICPSR 120246): columns mu and sd, as in application_MM_all_four_rounding.R (CFG$brodeur).
#  This script replaces only the engine's draw_meta() (how a meta-sample is drawn) and draws
#  its own figures; every test is the engine's one_rep().
#
#  Outputs in OUT_DIR (NAME = CFG$name):
#    NAME_S1000.pdf/.png, NAME_S5000.pdf/.png   size of all tests against coarseness, one per S
#    NAME_S1000_EKW.pdf/.png, ...               the same for each EKW test separately
#    NAME_tcurves.pdf/.png                      what rounding does to the t-curve (first panel)
#    NAME_S1000.tex, NAME_S5000.tex             LaTeX tables (booktabs)
#    results_NAME.csv/.rds                      rejection rates by S, row, panel and c, with the
#                                               rounding statistics of each pool
#    rounding_summary_NAME.csv                  rounding statistics of each pool and c
#    brodeur_rounding_profile.rds               significands of the Brodeur et al. data (cache)
#  Finished cells are saved in OUT_DIR/cells and pools in OUT_DIR/pools, so an interrupted run
#  resumes where it stopped.
#
#  Run once with QUICK <- TRUE (minutes) to check that everything works, then with FALSE.
#  Without the data, CFG$brodeur$synthetic <- TRUE runs the pipeline on made-up significands
#  (output folder suffix "_synthetic"); those results say nothing about Brodeur et al.
################################################################################

rm(list = ls())

## ---- 0. Settings -------------------------------------------------------------------------

ROOT     <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
CODE_DIR <- file.path(ROOT, "code")

QUICK <- TRUE   # TRUE: small run that checks the pipeline; FALSE: the settings below

CFG <- list(
  name           = "rounding_size",       # names the output folder and files
  S_values       = c(5000),         # one figure per S (the four Brodeur et al. subsamples
  # have 3,148 to 7,569 t-statistics)
  severity       = c(0, 1/2, 1, 2),   # coarseness c of rounding (see header)
  x_scale        = "log2",                # "log2" (c = 0 drawn one step left of the smallest
  # c > 0) or "linear"
  reps           = 1000,                  # Monte Carlo repetitions per (S, row, panel, c)
  alpha          = 0.05,                  # test level
  
  # Rows, top to bottom: any of "raw", "deround", "smin", "smin_noise"
  rows           = c("raw", "deround", "smin_noise"),
  row_titles     = NULL,                  # optional, e.g. c(raw = "Rounded tables")
  
  # Panels, left to right: distributions of true effects (as in the companion script)
  panels         = list(
    list(id = "pihat",    pi0 = "kudrin_pi_hat"),
    list(id = "twopoint", pi0 = list(shape = "twopoint", h_center = 0)) #,
    #list(id = "N07",      pi0 = list(shape = "normal", h_center = 2, sigma_h = 0.7))
  ),
  
  n_pool         = 1e6,                   # studies per pool
  nu             = 9999999,               # degrees of freedom of the noise e (run_sims' nu)
  
  # Rounding
  zero_se_rule   = "more_decimals",       # se that would show as zero: "more_decimals" or "drop"
  kappa_cap      = 0.3,                   # every row drops results with kappa = a / s_min above
  # this (paper, Section 5.2.1: 0.3; the application script
  # has kappa_cap <- 0.5); Inf keeps them all
  
  # Brodeur, Cook and Heyes (2020) data: calibration and the vertical line
  brodeur        = list(
    file           = "C:/Users/sfaridani6/Documents/Data/120246-V1/MM data.dta",                # NULL: the first of `candidates` that exists; or a path
    # to "MM data.dta" (.csv and .rds also work)
    candidates     = c("C:/Users/sfaridani6/Documents/Data/120246-V1/MM data.dta",  # as in the
                       file.path(ROOT, c("data/MM data.dta", "data/MM Data.dta",     # application
                                         "Data/MM data.dta", "MM data.dta", "code/MM data.dta"))),
    mu_var         = "mu",                # reported coefficient (columns as in the application)
    se_var         = "sd",                # reported standard error
    method_var     = "method",            # DID / IV / RCT / RDD (only for the printed summary)
    methods        = NULL,                # NULL: all results; or e.g. c("DID", "IV")
    report_var     = "ireport",           # how a result was reported (se, t, p, CI); its table
    # is printed on the first read, to check sd
    report_keep    = NULL,                # NULL: every row with a coefficient and an se; or
    # the report_var values to keep (those for an se)
    condition_on   = "effect",            # draw significands given the effect |h| ("effect"),
    # the realized |t| ("t") or nothing ("none")
    t_bins         = 10,                  # number of |t| quantile bins for that
    storage_digits = 6,                   # significant digits kept before counting decimals
    # (6, as n_dec() in the application)
    synthetic      = FALSE                # TRUE: made-up significands, ONLY to test the pipeline
  ),
  
  # Kudrin tests: your kudrin_tests()
  kudrin_boots   = 199,                   # T_inf bootstrap draws
  eig_floor      = 1e-8,                  # hybrid metric floor
  
  # Projection test: your run_test_weight_fast(), with the paper's settings
  run_projection = TRUE,
  num_coeffs     = 41,                    # J = 40
  J0             = 20,                    # weights w_j = min(1, (1/p)^(j - J0))
  p_weight       = 2,
  shift          = 1.96,
  L              = 6.5,
  numgrid        = 3000,
  proj_boots     = 100,
  proj_mode      = "full",                # "full" as in the paper, or "warp" (pooled critical value)
  epsilon_U_pre  = NULL,                  # NULL: compute_epsilons(L, numgrid) once and cache it
  # (slow); or 0.00063, the value in the paper
  
  # EKW tests (Tests.R)
  run_ekw        = TRUE,
  ekw_J          = 30,
  
  # Figures
  ylim           = NULL,                  # NULL: per row, from the rejection rates; or c(0, 1)
  mc_band        = TRUE,                  # shade the 95% simulation band around the level
  figure_note    = TRUE,                  # two lines of text above the panels
  
  # Computing
  n_cores        = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  seed           = 20261028,
  resume         = TRUE,
  run_selfchecks = TRUE
)

OUT_DIR <- file.path(ROOT, "simulations", CFG$name)

if (QUICK) {
  CFG$reps          <- 100
  CFG$n_pool        <- 2e5
  CFG$kudrin_boots  <- 99
  CFG$proj_mode     <- "warp"
  CFG$epsilon_U_pre <- 0.00063
  OUT_DIR           <- paste0(OUT_DIR, "_quick")
}
if (isTRUE(CFG$brodeur$synthetic)) OUT_DIR <- paste0(OUT_DIR, "_synthetic")

if (!CFG$zero_se_rule %in% c("more_decimals", "drop")) {
  stop("CFG$zero_se_rule must be \"more_decimals\" or \"drop\".")
}
if (!CFG$brodeur$condition_on %in% c("effect", "t", "none")) {
  stop("CFG$brodeur$condition_on must be \"effect\", \"t\" or \"none\".")
}
if (!is.numeric(CFG$kappa_cap) || CFG$kappa_cap <= 0) stop("CFG$kappa_cap must be a number > 0 (or Inf).")
SEV_GRID <- sort(unique(as.numeric(CFG$severity)))
if (length(SEV_GRID) == 0 || any(!is.finite(SEV_GRID)) || any(SEV_GRID < 0)) {
  stop("CFG$severity must be finite numbers >= 0.")
}


## ---- 1. Rows: what the meta-analyst does with the rounded tables ------------------------
## column: the t-statistic a row reports (draw_meta() below). All rows of a panel draw from
## the same pool with the same seed, so they test the same studies at every c.

ROUNDING_ROWS <- list(
  raw        = list(column = "raw", short = "None",
                    title  = "No de-rounding: t = m / s",
                    latex  = "No de-rounding: $t = m/s$"),
  deround    = list(column = "deround", short = "Uniform",
                    title  = "Uniform de-rounding: t = (m + V) / W",
                    latex  = "Uniform de-rounding: $t = (m + V)/W$, $W \\sim U[s_{\\min}, s_{\\max}]$"),
  smin       = list(column = "smin", short = "s_min",
                    title  = "Smallest consistent SE: t = m / s_min",
                    latex  = "Smallest consistent standard error: $t = m/s_{\\min}$"),
  smin_noise = list(column = "smin_noise", short = "s_min + noise",
                    title  = "Smallest consistent SE + numerator noise: t = (m + V) / s_min",
                    latex  = "Smallest consistent standard error with numerator noise: $t = (m + V)/s_{\\min}$")
)
bad <- setdiff(CFG$rows, names(ROUNDING_ROWS))
if (length(bad) > 0) {
  stop("Unknown row(s): ", paste(bad, collapse = ", "), ". Choose from ",
       paste(names(ROUNDING_ROWS), collapse = ", "), ".")
}
ROWS <- lapply(CFG$rows, function(m) {
  rw <- c(list(id = m), ROUNDING_ROWS[[m]])
  if (!is.null(CFG$row_titles) && m %in% names(CFG$row_titles)) rw$title <- CFG$row_titles[[m]]
  rw
})
names(ROWS) <- CFG$rows


## ---- 2. Your code ------------------------------------------------------------------------

CODE_FILES  <- file.path(CODE_DIR, c(if (CFG$run_ekw) "Tests.R",
                                     "reproduce_kudrin_simulations_functions.R",
                                     "functions_all_weighted.R",
                                     "kudrin_tests.R"))
ENGINE_FILE <- file.path(CODE_DIR, "kudrin_figures_engine.R")
missing_files <- c(CODE_FILES, ENGINE_FILE)[!file.exists(c(CODE_FILES, ENGINE_FILE))]
if (length(missing_files) > 0) stop("Not found: ", paste(missing_files, collapse = ", "))
for (f in CODE_FILES) source(f)
source(ENGINE_FILE)
stopifnot(exists("kudrin_tests"), exists("run_test_weight_fast"), exists("create_basis_stable"),
          exists("compute_epsilons"), exists("one_rep"), exists("run_cell"),
          exists("summarise_cell"), exists("start_cluster"), exists("draw_series"))
dir.create(file.path(OUT_DIR, "cells"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "pools"), recursive = TRUE, showWarnings = FALSE)

if (CFG$run_selfchecks && exists("kudrin_tests_selfcheck")) {
  kudrin_tests_selfcheck()                # Kudrin tests: Si, I(a), E[G_tilde], QP, full runs
}


## ---- 3. Distributions of true effects -----------------------------------------------------

PI0_ALIASES <- list(                     # run_sims() pi0_shape names
  null            = list(shape = "point",   h_center = 0),
  onehalf         = list(shape = "point",   h_center = 0.5),
  two             = list(shape = "point",   h_center = 2),
  standard_normal = list(shape = "normal",  h_center = 0, sigma_h = 1),
  normal_11       = list(shape = "normal",  h_center = 1, sigma_h = 1),
  hard            = list(shape = "normal",  h_center = 2, sigma_h = 0.7),
  unif_half       = list(shape = "uniform", h_center = 0, sigma_h = 1),
  "poisson(2)"    = list(shape = "poisson", h_center = 2),
  "chi2(2)"       = list(shape = "chi2",    h_center = 2),
  fitted          = list(shape = "fitted"),
  kudrin_gamma    = list(shape = "gamma",   gamma_shape = 1.4134, gamma_rate = 1.2980),
  kudrin_pi_hat   = list(shape = "gamma",   gamma_shape = 1.4134, gamma_scale = 1.2980,
                         label = quote(italic(h) %~% hat(Pi)))
)

as_pi0 <- function(x) {
  if (is.function(x)) return(list(shape = "custom", fun = x))
  if (is.character(x) && length(x) == 1) {
    if (is.null(PI0_ALIASES[[x]])) {
      stop("Unknown pi0 name: ", x, ". Known names: ", paste(names(PI0_ALIASES), collapse = ", "))
    }
    return(PI0_ALIASES[[x]])
  }
  if (is.list(x) && !is.null(x$shape)) return(x)
  stop("A pi0 must be a run_sims() pi0_shape name, a list with a shape, or a function(n).")
}

# n draws of h, with run_sims()' formulas
draw_pi0 <- function(spec, n) {
  s <- spec$shape
  if (s == "custom") {
    h <- spec$fun(n)
  } else if (s == "point") {
    h <- rep(spec$h_center, n)
  }else if (s == "twopoint") {
    h <- rep(spec$h_center, n)+(runif(n)>=0.5)*2
  } else if (s == "normal") {
    h <- rnorm(n) * spec$sigma_h + spec$h_center
  } else if (s == "uniform") {
    h <- (runif(n) - 0.5) * spec$sigma_h + spec$h_center
  } else if (s == "chi2") {
    h <- rchisq(n, spec$h_center)
  } else if (s == "poisson") {
    h <- rpois(n, spec$h_center)
  } else if (s == "gamma") {
    if (is.null(spec$gamma_scale)) {
      h <- rgamma(n, shape = spec$gamma_shape, rate = spec$gamma_rate)
    } else {
      h <- rgamma(n, shape = spec$gamma_shape, scale = spec$gamma_scale)
    }
  } else if (s == "double_normal") {
    h  <- rnorm(n) * spec$sigma_h + spec$h_center
    rr <- runif(n)
    h[rr >= 0.5] <- h[rr >= 0.5] - 2
  } else if (s == "fitted") {             # fitted by zcurve, Faridani (2026)
    h <- sample(c(0, 1, 2, 3, 4, 5, 6), size = n, replace = TRUE,
                prob = c(0.02, 0.47, 0.01, 0.27, 0.14, 0.00, 0.09))
  } else {
    stop("Unknown pi0 shape: ", s)
  }
  h <- as.numeric(h)
  if (length(h) != n || any(!is.finite(h))) stop("pi0 '", s, "' did not return ", n, " finite draws.")
  h
}

# Panel title (plotmath)
pi0_label <- function(spec) {
  if (!is.null(spec$label)) return(as.expression(spec$label))
  s  <- spec$shape
  g4 <- function(x) signif(x, 4)
  lab <- if (s == "point") {
    bquote(italic(h) == .(g4(spec$h_center)))
  } else if (s == "normal") {
    bquote(italic(h) %~% italic(N)*"("*.(g4(spec$h_center))*","~.(g4(spec$sigma_h))^2*")")
  } else if (s == "twopoint") {
    bquote(italic(h) %~% "0 and 2 mix" )
  } else if (s == "uniform") {
    bquote(italic(h) %~% italic(U)*"["*.(g4(spec$h_center - spec$sigma_h / 2))*","~
             .(g4(spec$h_center + spec$sigma_h / 2))*"]")
  } else if (s == "chi2") {
    bquote(italic(h) %~% chi[.(g4(spec$h_center))]^2)
  } else if (s == "poisson") {
    bquote(italic(h) %~% "Poisson("*.(g4(spec$h_center))*")")
  } else if (s == "gamma" && is.null(spec$gamma_scale)) {
    bquote(italic(h) %~% "Gamma("*.(g4(spec$gamma_shape))*", rate"~.(g4(spec$gamma_rate))*")")
  } else if (s == "gamma") {
    bquote(italic(h) %~% "Gamma("*.(g4(spec$gamma_shape))*", scale"~.(g4(spec$gamma_scale))*")")
  } else if (s == "double_normal") {
    bquote(italic(h) %~% "mixture of"~italic(N)*"("*.(g4(spec$h_center))*","~.(g4(spec$sigma_h))^2*")"~
             "and"~italic(N)*"("*.(g4(spec$h_center - 2))*","~.(g4(spec$sigma_h))^2*")")
  } else if (s == "fitted") {
    quote(italic(h) %~% "z-curve fit")
  } else {
    quote(italic(h) %~% "custom")
  }
  as.expression(lab)
}

PANELS <- lapply(seq_along(CFG$panels), function(k) {
  p    <- CFG$panels[[k]]
  spec <- as_pi0(p$pi0)
  list(id    = if (is.null(p$id)) sprintf("pi%d", k) else p$id,
       spec  = spec,
       title = if (is.null(p$label)) pi0_label(spec) else as.expression(p$label))
})
PANEL_IDS <- vapply(PANELS, function(pn) pn$id, "")
if (anyDuplicated(PANEL_IDS) > 0) stop("Panel ids must be distinct: ", paste(PANEL_IDS, collapse = ", "))
names(PANELS) <- PANEL_IDS



## ---- 4. Rounding in Brodeur, Cook and Heyes (2020) ----------------------------------------
## The significands of the reported standard errors, with decimal places as in Kranz and Puetz
## (2022). This is what calibrates the simulated tables at c = 1.

# Decimal places of each element of x. Numbers are first rounded to the significant digits the
# file can store (Stata floats hold about 7), so float noise such as 0.0120000001788 counts as
# 0.012. Text keeps its trailing zeros and is counted as written.
num_decimals <- function(x, storage_digits = 7, max_d = 20) {
  out <- rep(NA_integer_, length(x))
  if (is.character(x)) {
    x     <- trimws(x)
    plain <- !is.na(x) & grepl("^[-+]?[0-9]*\\.?[0-9]*$", x) & grepl("[0-9]", x)
    point <- plain & grepl(".", x, fixed = TRUE)
    out[plain] <- 0L
    out[point] <- nchar(sub("^[^.]*\\.", "", x[point]))
    other <- !plain & !is.na(suppressWarnings(as.numeric(x)))      # e.g. "1e-04"
    if (any(other)) out[other] <- num_decimals(as.numeric(x[other]), storage_digits, max_d)
    return(out)
  }
  x  <- signif(abs(as.numeric(x)), storage_digits)
  ok <- is.finite(x)
  for (d in 0:max_d) {
    todo <- which(ok & is.na(out))
    if (length(todo) == 0) break
    y   <- x[todo] * 10^d
    hit <- abs(y - round(y)) < 5e-8
    out[todo[hit]] <- d
  }
  out
}

# Significand of the reported se, with the decimal places of Kranz and Puetz (2022): the larger
# of those of the coefficient and of the se. A table shows both to the same decimal place, so
# this restores trailing zeros of the se that the data do not store: 0.023 (0.01) -> 10.
kp_significand <- function(mu_raw, se_raw, storage_digits = 7) {
  d  <- pmax(num_decimals(mu_raw, storage_digits), num_decimals(se_raw, storage_digits),
             na.rm = TRUE)
  se <- signif(abs(suppressWarnings(as.numeric(se_raw))), storage_digits)
  round(se * 10^d)
}

strip_labels <- function(x) {
  if (inherits(x, "haven_labelled") && requireNamespace("haven", quietly = TRUE)) {
    x <- haven::zap_labels(x)
  }
  if (is.factor(x)) x <- as.character(x)
  as.vector(x)
}

as_text <- function(x) {
  if (inherits(x, "haven_labelled") && requireNamespace("haven", quietly = TRUE)) {
    return(as.character(haven::as_factor(x)))
  }
  as.character(x)
}

read_brodeur <- function(file) {
  B   <- CFG$brodeur
  ext <- tolower(tools::file_ext(file))
  if (ext == "dta") {
    if (requireNamespace("haven", quietly = TRUE)) return(haven::read_dta(file))
    if (requireNamespace("readstata13", quietly = TRUE)) {
      return(readstata13::read.dta13(file, convert.factors = FALSE))
    }
    stop("Reading ", file, " needs the haven package: install.packages(\"haven\").")
  }
  if (ext == "csv") {                     # read b and se as text, so trailing zeros survive
    return(utils::read.csv(file, stringsAsFactors = FALSE,
                           colClasses = setNames(rep("character", 2), c(B$mu_var, B$se_var))))
  }
  if (ext == "rds") return(readRDS(file))
  stop("Unknown file type (use .dta, .csv or .rds): ", file)
}

# Edges of nb quantile bins of the reported |t|, from 0 to Inf; bins without data are merged
t_bin_edges <- function(absT, nb) {
  if (is.null(nb) || nb <= 1) return(c(0, Inf))
  e <- unique(quantile(absT, probs = seq(0, 1, length.out = nb + 1), type = 1, names = FALSE))
  e <- unique(c(0, e[-c(1, length(e))], Inf))
  repeat {
    cnt <- tabulate(findInterval(absT, e), nbins = length(e) - 1)
    if (all(cnt > 0) || length(e) <= 2) break
    j <- which(cnt == 0)[1]
    e <- e[-(if (j < length(e) - 1) j + 1 else j)]
  }
  e
}

finish_profile <- function(pr) {
  pr$t_edges <- t_bin_edges(pr$absT, if (CFG$brodeur$condition_on != "none") CFG$brodeur$t_bins else 1)
  pr
}

brodeur_profile <- function() {
  B <- CFG$brodeur
  if (isTRUE(B$synthetic)) {
    warning("CFG$brodeur$synthetic = TRUE: made-up significands, NOT the Brodeur et al. data.")
    set.seed(CFG$seed + 7)
    n <- 20000
    return(finish_profile(list(key = paste("synthetic", CFG$seed), file = "SYNTHETIC",
                               s = pmax(1, round(10^runif(n, 0, 3.3))), absT = abs(rnorm(n, 1.5, 2)),
                               method = rep("synthetic", n), n_file = n, n_coef_se = n)))
  }
  file <- B$file
  if (is.null(file)) {
    cand <- B$candidates
    file <- cand[file.exists(cand)][1]
  }
  if (is.na(file) || !file.exists(file)) {
    stop("The Brodeur, Cook and Heyes (2020) data were not found. Set CFG$brodeur$file to your copy ",
         "of \"MM data.dta\" (their replication package, openICPSR project 120246), or set ",
         "CFG$brodeur$synthetic <- TRUE to test the pipeline without it.")
  }
  info  <- file.info(file)
  key   <- paste(normalizePath(file), info$size, format(info$mtime), B$mu_var, B$se_var,
                 B$method_var, paste(B$methods, collapse = ","), B$storage_digits,
                 B$report_var, paste(B$report_keep, collapse = ","), "mu != 0")
  cache <- file.path(OUT_DIR, "brodeur_rounding_profile.rds")
  if (file.exists(cache)) {
    pr <- readRDS(cache)
    if (identical(pr$key, key)) return(finish_profile(pr))
  }
  df   <- read_brodeur(file)
  miss <- setdiff(c(B$mu_var, B$se_var), names(df))
  if (length(miss) > 0) {
    stop("Column(s) ", paste(miss, collapse = ", "), " not in ", file, ". Its columns are: ",
         paste(names(df), collapse = ", "), ". Set CFG$brodeur$mu_var and CFG$brodeur$se_var.")
  }
  mu_raw <- strip_labels(df[[B$mu_var]])
  se_raw <- strip_labels(df[[B$se_var]])
  meth   <- if (!is.null(B$method_var) && B$method_var %in% names(df)) {
    as_text(df[[B$method_var]])
  } else {
    rep("all", nrow(df))
  }
  mu   <- suppressWarnings(as.numeric(mu_raw))
  se   <- suppressWarnings(as.numeric(se_raw))
  keep <- is.finite(mu) & is.finite(se) & se > 0 & mu != 0   # a coefficient and se; mu = 0:
  # Kudrin's drop, as in the application
  if (!is.null(B$methods)) keep <- keep & meth %in% B$methods
  if (!is.null(B$report_var) && B$report_var %in% names(df)) {
    rep_how <- as_text(df[[B$report_var]])
    message("Rows by ", B$report_var, " (columns: has a usable ", B$mu_var, " and ", B$se_var,
            "). If some values mean a t-statistic or p-value was reported, check what ",
            B$se_var, " holds for them, or set CFG$brodeur$report_keep:")
    print(table(rep_how, usable = keep, useNA = "ifany"))
    if (!is.null(B$report_keep)) keep <- keep & rep_how %in% as.character(B$report_keep)
  } else if (!is.null(B$report_keep)) {
    stop("CFG$brodeur$report_keep is set, but column ", B$report_var, " is not in ", file, ".")
  }
  s    <- kp_significand(mu_raw[keep], se_raw[keep], B$storage_digits)
  ok   <- is.finite(s) & s >= 1
  pr   <- list(key = key, file = file, s = s[ok], absT = abs(mu[keep] / se[keep])[ok],
               method = meth[keep][ok], n_file = nrow(df), n_coef_se = sum(keep))
  saveRDS(pr, cache)
  finish_profile(pr)
}

profile_stats <- function(s) {
  data.frame(n = length(s), one_digit = mean(s <= 9), below_37 = mean(s < 37),
             kappa_drop = mean(!kappa_ok(s)), median = median(s), mean_unit = mean(1 / s))
}

# one_digit: share of SEs with a single significant digit; below_37: share with significand
# below 37 (the observations Kranz and Puetz 2022 drop); kappa_drop: share above the cap
# CFG$kappa_cap (with shared decimals); mean_unit: mean rounding unit / se
print_profile <- function(pr) {
  message(sprintf("Rounding in Brodeur et al. (2020), %s: %s of %s rows report a coefficient and a standard error.",
                  basename(pr$file), format(length(pr$s), big.mark = ","),
                  format(pr$n_file, big.mark = ",")))
  groups <- split(pr$s, pr$method)
  tab    <- do.call(rbind, c(list(profile_stats(pr$s)), lapply(groups, profile_stats)))
  rownames(tab) <- c("all", names(groups))
  print(tab, digits = 3)
  if (length(pr$t_edges) > 2) {
    nb  <- length(pr$t_edges) - 1
    b   <- factor(findInterval(pr$absT, pr$t_edges), levels = seq_len(nb))
    byt <- data.frame(abs_t_from = pr$t_edges[-(nb + 1)], abs_t_to = pr$t_edges[-1],
                      n = as.vector(table(b)),
                      one_digit = as.vector(tapply(pr$s <= 9, b, mean)),
                      below_37  = as.vector(tapply(pr$s < 37, b, mean)))
    message("The same by bin of the reported |t| (a simulated study draws its significand from the bin ",
            "holding its ", if (CFG$brodeur$condition_on == "t") "realized |t|" else "effect |h|", "):")
    print(byt, digits = 3, row.names = FALSE)
  }
  invisible(tab)
}


## ---- 5. Rounded tables: the DGP and the pools ----------------------------------------------

# Round the tables of studies with honest |t| = t and standard error xB (in units of the last
# reported digit at c = 1) at coarseness sev > 0. Returns the significand s of the reported se,
# the reported coefficient m in the same units, x, and keep (FALSE: the se shows as zero, which
# can happen only under the rule "drop").
round_results <- function(t, xB, sev, rule) {
  x <- xB / sev
  if (identical(rule, "more_decimals")) {        # a se that would show as 0 gets more decimals
    low <- round(x) < 1
    while (any(low)) {
      x[low] <- 10 * x[low]
      low    <- round(x) < 1
    }
  }
  s <- round(x)
  list(s = s, m = round(t * x), x = x, keep = s >= 1)
}

# Kept in every row (as in the application, the sample does not depend on the option): a se
# printed with significand s >= 1 and kappa = a / s_min = 1 / (s - 1/2) <= CFG$kappa_cap, where
# a = b = one unit of the last digit, which b and se share
kappa_ok <- function(s) s >= 1 & 1 / (s - 0.5) <= CFG$kappa_cap

# Reported |t| of row `column` from rounded tables r (round_results()), with the de-rounding
# noise u1 (numerator, V) and u2 (denominator) uniform on (-1/2, 1/2)
row_t <- function(column, r, u1, u2) {
  switch(column,
         raw        = r$m / r$s,
         deround    = abs((r$m + u1) / (r$s + u2)),
         smin       = r$m / (r$s - 0.5),
         smin_noise = abs((r$m + u1) / (r$s - 0.5)),
         stop("Unknown row column: ", column))
}

# Rows of the pool that every row analyzes at coarseness sev, cached for the current pool and sev
ELIG <- new.env()
eligible_rows <- function(pool, sev) {
  key <- paste(attr(pool, "key"), sev, CFG$kappa_cap, CFG$zero_se_rule)
  if (!identical(ELIG$key, key)) {
    s        <- round_results(numeric(nrow(pool)), pool[, "xB"], sev, CFG$zero_se_rule)$s
    ELIG$idx <- which(kappa_ok(s))
    ELIG$key <- key
  }
  ELIG$idx
}

# Replaces the engine's draw_meta(): S studies from those kept at coarseness tau (= c), each
# reporting the t-statistic of row `column`. Every row draws the same studies and consumes the
# same random numbers (the de-rounding noise is drawn in every row), so the rows of a cell test
# the same meta-samples, and at c = 0 (no rounding) they agree exactly. Returns |t|.
draw_meta <- function(pool, column, tau, S) {
  sev <- tau
  if (sev > 0 && (is.finite(CFG$kappa_cap) || identical(CFG$zero_se_rule, "drop"))) {
    elig <- eligible_rows(pool, sev)
    if (length(elig) == 0) stop("No result is kept at c = ", sev, ".")
    idx  <- elig[sample.int(length(elig), S, replace = TRUE)]
  } else {
    idx  <- sample.int(nrow(pool), S, replace = TRUE)
  }
  u <- runif(2 * S, -0.5, 0.5)                    # the meta-analyst's de-rounding noise
  t <- pool[idx, "honest"]
  if (sev == 0 || identical(column, "honest")) return(t)
  r <- round_results(t, pool[idx, "xB"], sev, CFG$zero_se_rule)
  row_t(column, r, u[seq_len(S)], u[S + seq_len(S)])
}

pool_id <- function(panel_id) paste0("round__", panel_id)

# Significands for simulated studies, drawn from the Brodeur et al. significands; absT is what
# they are drawn given (the effect |h|, or the realized |t|): a study draws from the results
# whose reported |t| lies in the same bin
draw_significands <- function(absT, prof) {
  if (identical(CFG$brodeur$condition_on, "none") || length(prof$t_edges) <= 2) {
    return(prof$s[sample.int(length(prof$s), length(absT), replace = TRUE)])
  }
  b_bch <- findInterval(prof$absT, prof$t_edges)
  b_sim <- findInterval(absT, prof$t_edges)
  out   <- numeric(length(absT))
  for (b in sort(unique(b_sim))) {
    pick      <- which(b_sim == b)
    src       <- which(b_bch == b)
    out[pick] <- prof$s[src[sample.int(length(src), length(pick), replace = TRUE)]]
  }
  out
}

# One pool per panel, cached in OUT_DIR/pools: the honest |t| of every study and its standard
# error in units of the last reported digit at c = 1, xB = s_BCH + U(-1/2, 1/2), so that
# round(xB) = s_BCH. Its seed depends only on the panel's position.
build_pool_rounding <- function(k, prof) {
  pn   <- PANELS[[k]]
  key  <- paste("rounding", paste(deparse(pn$spec), collapse = " "), CFG$n_pool, CFG$seed, k,
                CFG$nu, prof$key, CFG$brodeur$condition_on, CFG$brodeur$t_bins)
  file <- file.path(OUT_DIR, "pools", paste0(safe_id(pool_id(pn$id)), ".rds"))
  if (file.exists(file)) {
    pool <- readRDS(file)
    if (identical(attr(pool, "key"), key)) return(pool)
  }
  set.seed(CFG$seed + 1000 * k + 3)
  message(sprintf("Simulating the rounding pool for panel %s (%s studies) ...", pn$id,
                  format(CFG$n_pool, big.mark = ",", scientific = FALSE)))
  secs <- system.time({
    h     <- draw_pi0(pn$spec, CFG$n_pool)
    tstar <- abs(rt(CFG$n_pool, CFG$nu) + h)        # honest t, as ts1 in run_sims()
    s     <- draw_significands(if (CFG$brodeur$condition_on == "t") tstar else abs(h), prof)
    pool  <- cbind(honest = tstar, xB = s + runif(CFG$n_pool, -0.5, 0.5))
  })[["elapsed"]]
  attr(pool, "key") <- key
  saveRDS(pool, file)
  message(sprintf("  done in %.0f s", secs))
  pool
}

# What rounding does to a pool at each c (first n_max studies): significands of the reported
# SEs, the share left out by the cap, and among the kept studies the mean stretch se / s_min,
# the share of t = m/s exactly 2, mean |m/s - T*|, and the share pushed across 1.96 either way
rounding_summary <- function(pool, sevs, n_max = 2e5) {
  n  <- min(nrow(pool), n_max)
  t  <- pool[seq_len(n), "honest"]
  xB <- pool[seq_len(n), "xB"]
  do.call(rbind, lapply(sevs, function(sev) {
    if (sev == 0) {
      return(data.frame(severity = 0, one_digit = 0, below_37 = 0, median_s = Inf, mean_unit = 0,
                        dropped = 0, kappa_dropped = 0, mean_lambda = 1, at_2 = 0,
                        mean_abs_err = 0, crosses_1.96 = 0))
    }
    r  <- round_results(t, xB, sev, CFG$zero_se_rule)
    k  <- r$keep
    s  <- r$s[k]
    kk <- kappa_ok(r$s)                               # the studies every row analyzes
    tr <- r$m[kk] / r$s[kk]
    data.frame(severity = sev, one_digit = mean(s <= 9), below_37 = mean(s < 37),
               median_s = median(s), mean_unit = mean(1 / s), dropped = mean(!k),
               kappa_dropped = mean(!kk), mean_lambda = mean(r$x[kk] / (r$s[kk] - 0.5)),
               at_2 = mean(abs(tr - 2) < 1e-9), mean_abs_err = mean(abs(tr - t[kk])),
               crosses_1.96 = mean((tr > 1.96) != (t[kk] > 1.96)))
  }))
}

rounding_selfcheck <- function() {
  chk <- function(ok, what) if (!isTRUE(all(ok))) stop("Rounding self-check failed: ", what)
  chk(identical(num_decimals(c(0.012, 1.5, 2, 1e-4, 123.456, 0.0120000001788139,
                               1234.5670166015625, NA)), c(3L, 1L, 0L, 4L, 3L, 3L, 3L, NA)),
      "decimal places of numbers (including Stata float noise)")
  chk(identical(num_decimals(c("0.010", "2", "1.50", "-0.0300", NA)), c(3L, 0L, 2L, 4L, NA)),
      "decimal places of text")
  chk(identical(kp_significand(c(0.023, 1.25, 5), c(0.01, 0.3, 2)), c(10, 30, 2)),
      "Kranz-Puetz significands")
  pool <- POOLS[[1]]
  n    <- min(nrow(pool), 1e5)
  t    <- pool[seq_len(n), "honest"]
  xB   <- pool[seq_len(n), "xB"]
  r1   <- round_results(t, xB, 1, CFG$zero_se_rule)
  chk(r1$s == round(xB) & r1$s >= 1, "at c = 1 the significands are those of Brodeur et al.")
  for (sev in SEV_GRID[SEV_GRID > 0]) {
    r <- round_results(t, xB, sev, CFG$zero_se_rule)
    k <- r$keep
    chk(r$s[k] >= 1, sprintf("significands >= 1 at c = %g", sev))
    chk(abs(r$m[k] / r$s[k] - t[k]) <= (1 + t[k]) / (2 * r$s[k]) + 1e-9,
        sprintf("|t - T*| <= (1 + T*) / (2 s) at c = %g", sev))
    chk(r$x[k] >= r$s[k] - 0.5 - 1e-12,
        sprintf("true se >= s_min (lambda >= 1, Proposition 4) at c = %g", sev))
    sp <- round_results(numeric(nrow(pool)), pool[, "xB"], sev, CFG$zero_se_rule)$s
    e  <- eligible_rows(pool, sev)
    chk(length(e) == sum(kappa_ok(sp)) && all(kappa_ok(sp[e])),
        sprintf("the kept studies are those under the cap at c = %g", sev))
  }
  for (sev in unique(c(0, SEV_GRID))) {
    out <- lapply(ROWS, function(rw) {
      set.seed(1)
      list(t = draw_meta(pool, rw$column, sev, 200), u = runif(1))
    })
    for (o in out) {
      chk(identical(o$u, out[[1]]$u), sprintf("rows consume the same random numbers at c = %g", sev))
      chk(is.finite(o$t) & o$t >= 0, sprintf("finite |t| at c = %g", sev))
      if (sev == 0) chk(identical(o$t, out[[1]]$t), "rows agree at c = 0")
    }
  }
  message("Rounding self-checks passed.")
}


## ---- 6. Figures ----------------------------------------------------------------------------
## The engine's look (colours, markers, line types, legend strip), with coarseness c on the
## horizontal axis, the rounding of Brodeur et al. (2020) as a vertical line at c = 1, and the
## 95% simulation band around the level shaded.

STYLE$bch         <- list(col = "#009E73", lty = "solid", lwd = 1.2)
LEG_LABEL$bch     <- "Rounding in Brodeur et al. (2020)"
LEG_LABEL$EKW_max <- "Highest rejection rate of EKW tests"

EKW_TESTS <- c("EKW_Binomial", "EKW_Discontinuity", "EKW_Fisher", "EKW_LCM", "EKW_CS1",
               "EKW_CSUB", "EKW_CS2B")
EKW_STYLE <- list(
  EKW_Binomial      = list(col = "#E69F00", lty = "solid",    lwd = 1.1, pch = 1,  cex = 0.60, plwd = 0.9),
  EKW_Discontinuity = list(col = "#0072B2", lty = "dashed",   lwd = 1.1, pch = 2,  cex = 0.60, plwd = 0.9),
  EKW_Fisher        = list(col = "#CC79A7", lty = "dotdash",  lwd = 1.1, pch = 5,  cex = 0.60, plwd = 0.9),
  EKW_LCM           = list(col = "#D55E00", lty = "longdash", lwd = 1.1, pch = 6,  cex = 0.60, plwd = 0.9),
  EKW_CS1           = list(col = "#56B4E9", lty = "twodash",  lwd = 1.1, pch = 3,  cex = 0.60, plwd = 0.9),
  EKW_CSUB          = list(col = "#8C8C8C", lty = "dotted",   lwd = 1.4, pch = 4,  cex = 0.60, plwd = 0.9),
  EKW_CS2B          = list(col = "#000000", lty = "solid",    lwd = 1.4, pch = 16, cex = 0.50, plwd = 1.0)
)
EKW_LABEL <- list(EKW_Binomial = "Binomial", EKW_Discontinuity = "Discontinuity",
                  EKW_Fisher = "Fisher", EKW_LCM = "LCM", EKW_CS1 = "CS1", EKW_CSUB = "CSUB",
                  EKW_CS2B = "CS2B")

# Horizontal position of coarseness v: log2(v), with c = 0 one step left of the smallest c > 0
sev_pos <- function(v) {
  if (!identical(CFG$x_scale, "log2")) return(v)
  pos  <- SEV_GRID[SEV_GRID > 0]
  zero <- if (length(pos) > 0) log2(min(pos)) - 1 else 0
  ifelse(v > 0, log2(pmax(v, 1e-300)), zero)
}

sev_label <- function(v) {
  vapply(v, function(a) {
    if (a == 0) return("0")
    if (a < 1 && abs(1 / a - round(1 / a)) < 1e-8) return(paste0("1/", round(1 / a)))
    format(signif(a, 4), trim = TRUE, scientific = FALSE)
  }, "")
}

sev_tag <- function(v) gsub("[^0-9A-Za-z]", "p", format(signif(v, 6), scientific = FALSE, trim = TRUE))

y_ticks <- function(ylim) {
  yt <- pretty(ylim, n = 4)
  yt[yt >= ylim[1] - 1e-9 & yt <= ylim[2] + 1e-9]
}

y_labels <- function(yt) {
  lab <- sub("^0\\.", ".", formatC(yt, format = "fg", digits = 3))
  lab[abs(yt) < 1e-12] <- "0"
  lab
}

# Vertical range of one row of panels: CFG$ylim, or from 0 to just above the highest rate (at
# least .2, so that the level stays visible)
row_ylim <- function(rr, series) {
  if (!is.null(CFG$ylim)) return(CFG$ylim)
  v   <- unlist(rr[intersect(series, names(rr))], use.names = FALSE)
  v   <- v[is.finite(v)]
  top <- max(c(0.2, 1.1 * v))
  c(0, min(1, ceiling(20 * top) / 20))
}

draw_panel_rounding <- function(r, title, series, styles, xlim, ylim) {
  r  <- r[order(r$severity), , drop = FALSE]
  xs <- sev_pos(SEV_GRID)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i")
  if (isTRUE(CFG$mc_band) && nrow(r) > 0) {
    half <- qnorm(0.975) * sqrt(CFG$alpha * (1 - CFG$alpha) / min(r$reps))
    rect(xlim[1], max(ylim[1], CFG$alpha - half), xlim[2], min(ylim[2], CFG$alpha + half),
         col = "grey92", border = NA)
  }
  yt <- y_ticks(ylim)
  abline(h = yt[yt > ylim[1] & yt < ylim[2]], v = xs, col = GRID_COL, lwd = 0.5)
  axis(1, at = xs, labels = FALSE, lwd = 0, lwd.ticks = 0.5, col.ticks = AXIS_COL)
  mtext(sev_label(SEV_GRID), side = 1, at = xs, line = 0.3, cex = 0.66, col = AXIS_COL)
  axis(2, at = yt, labels = y_labels(yt), lwd = 0, lwd.ticks = 0.5, col.ticks = AXIS_COL,
       col.axis = AXIS_COL, cex.axis = 0.8)
  box(lwd = 0.6, col = AXIS_COL)
  mtext(expression("Coarseness of rounding ("*italic(c)*")"), side = 1, line = 1.35, cex = 0.78)
  mtext("Rejection rate", side = 2, line = 1.9, cex = 0.78, las = 0)
  mtext(title, side = 3, line = 0.25, cex = 0.8)
  abline(h = CFG$alpha, col = STYLE$sig$col, lty = STYLE$sig$lty, lwd = STYLE$sig$lwd)
  abline(v = sev_pos(1), col = STYLE$bch$col, lty = STYLE$bch$lty, lwd = STYLE$bch$lwd)
  for (tn in series) {
    if (tn %in% names(r)) draw_series(sev_pos(r$severity), r[[tn]], styles[[tn]])
  }
  invisible(NULL)
}

# As the engine's draw_legend_strip(), for any list of styles and labels
draw_legend_rounding <- function(keys, styles, labels, ncol) {
  sty  <- styles[keys]
  labs <- do.call(c, lapply(keys, function(k) as.expression(labels[[k]])))
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

# One figure: rows x panels at sample size S; series: result columns to draw
plot_rounding_figure <- function(res, S, rows, panels, stem, series, styles, labels, keys,
                                 note = NULL) {
  nr        <- length(rows)
  np        <- length(panels)
  mid       <- ceiling(np / 2)
  keys      <- keys[keys %in% names(styles)]
  ncol      <- if (np >= 2) 3 else 2
  leg_lines <- 0.9 + 1.1 * ceiling(length(keys) / ncol)
  top_lines <- if (is.null(note)) 0 else 0.3 + 1.1 * length(note)
  line_in   <- 8 * 1.2 / 72                         # height of one margin line at pointsize 8
  height    <- 2.05 * nr + (leg_lines + top_lines) * line_in
  width     <- max(3.3, 2.2 * np)
  xs        <- sev_pos(SEV_GRID)
  pad       <- if (identical(CFG$x_scale, "log2")) 0.45 else 0.04 * max(1, diff(range(xs)))
  xlim      <- range(xs) + c(-pad, pad)
  for (ext in c("pdf", "png")) {
    open_device(paste0(stem, ".", ext), ext, height = height, width = width)
    par(mfrow = c(nr, np))                          # separate call: mfrow resets cex
    par(cex = 1, mar = c(2.6, 3.0, 2.5, 0.8), oma = c(leg_lines, 0, top_lines, 0),
        mgp = c(1.5, 0.35, 0), tcl = -0.22, las = 1, family = "sans")
    for (i in seq_len(nr)) {
      rr   <- res[res$S == S & res$row == rows[[i]]$id, , drop = FALSE]
      ylim <- row_ylim(rr, series)
      for (k in seq_len(np)) {
        r <- rr[rr$panel == panels[[k]]$id, , drop = FALSE]
        draw_panel_rounding(r, panels[[k]]$title, series, styles, xlim, ylim)
        if (k == mid) mtext(rows[[i]]$title, side = 3, line = 1.2, font = 2, cex = 1.05)
      }
    }
    for (j in seq_along(note)) {
      mtext(note[j], side = 3, outer = TRUE, line = 0.3 + 1.1 * (length(note) - j), cex = 0.72)
    }
    draw_legend_rounding(keys, styles, labels, ncol)
    dev.off()
  }
  invisible(paste0(stem, c(".pdf", ".png")))
}

figure_note <- function(S, res) {
  if (!isTRUE(CFG$figure_note)) return(NULL)
  reps <- paste(format(sort(unique(res$reps[res$S == S])), big.mark = ","), collapse = "/")
  c(sprintf("S = %s, %s repetitions per point. Nothing is selectively reported, so every rejection is false. Shaded: 95%% simulation band around %s%%.",
            format(S, big.mark = ",", scientific = FALSE), reps, format(100 * CFG$alpha)),
    sprintf("c: rounding unit / SE as a multiple of Brodeur et al. (2020) (green line, c = 1; c = 0: no rounding)%s.%s",
            if (is.finite(CFG$kappa_cap)) sprintf("; every row drops kappa > %s", format(CFG$kappa_cap)) else "",
            if (isTRUE(CFG$brodeur$synthetic)) " [SYNTHETIC significands]" else ""))
}

# Main figure (projection test, Kudrin's three tests, highest EKW rate) and the same for each
# EKW test
make_figures <- function(res, S) {
  stem   <- file.path(OUT_DIR, sprintf("%s_S%d", CFG$name, as.integer(S)))
  series <- intersect(c("EKW_max", "T_inf", "T_H", "T_H_hybrid", "proj"), names(res))
  plot_rounding_figure(res, S, ROWS, PANELS, stem, series, STYLE, LEG_LABEL,
                       keys = c("sig", "bch", series), note = figure_note(S, res))
  ekw <- intersect(EKW_TESTS, names(res))
  if (length(ekw) > 0) {
    plot_rounding_figure(res, S, ROWS, PANELS, paste0(stem, "_EKW"), ekw,
                         styles = c(STYLE[c("sig", "bch")], EKW_STYLE),
                         labels = c(LEG_LABEL[c("sig", "bch")], EKW_LABEL),
                         keys = c("sig", "bch", ekw), note = figure_note(S, res))
  }
  invisible(stem)
}

# What rounding does to the t-curve: histograms of each row's reported |t| (bins of .02) against
# the honest |t| of the same kept studies (black), at the coarseness values in sevs
plot_tcurves <- function(pool, stem, sevs, panel_title, n_max = 2e5) {
  sevs <- sevs[sevs > 0]
  if (length(sevs) == 0) return(invisible(NULL))
  n    <- min(nrow(pool), n_max)
  t    <- pool[seq_len(n), "honest"]
  xB   <- pool[seq_len(n), "xB"]
  set.seed(CFG$seed + 99)
  u1   <- runif(n, -0.5, 0.5)
  u2   <- runif(n, -0.5, 0.5)
  br   <- seq(0, 5, by = 0.02)
  mids <- br[-1] - 0.01
  dens <- function(v) {
    nv <- length(v)
    v  <- v[is.finite(v) & v < 5]
    hist(v, breaks = br, plot = FALSE)$counts / (nv * 0.02)
  }
  cols <- c("#0072B2", "#CC79A7", "#009E73", "#D55E00")
  nc   <- length(ROWS)
  for (ext in c("pdf", "png")) {
    open_device(paste0(stem, ".", ext), ext, height = 1.9 * length(sevs) + 0.3,
                width = max(6.6, 2.2 * nc))
    par(mfrow = c(length(sevs), nc))
    par(cex = 1, mar = c(2.6, 3.0, 1.8, 0.8), oma = c(0, 0, 1.3, 0), mgp = c(1.5, 0.35, 0),
        tcl = -0.22, las = 1, family = "sans")
    for (sev in sevs) {
      r  <- round_results(t, xB, sev, CFG$zero_se_rule)
      kk <- kappa_ok(r$s)
      if (sum(kk) < 100) {                           # (nearly) everything above the cap
        for (j in seq_len(nc)) plot.new()
        next
      }
      hon <- dens(t[kk])
      for (j in seq_len(nc)) {
        v   <- row_t(ROWS[[j]]$column, r, u1, u2)[kk]
        d   <- dens(v)
        top <- 1.05 * min(max(d), 3 * max(hon))
        plot(mids, d, type = "h", lwd = 0.7, col = cols[(j - 1) %% length(cols) + 1],
             xlim = c(0, 5), ylim = c(0, top), xaxs = "i", yaxs = "i", xlab = "", ylab = "",
             axes = FALSE)
        lines(mids, hon, col = "black", lwd = 0.9)
        abline(v = 1.96, lty = 3, lwd = 0.8)
        axis(1, lwd = 0, lwd.ticks = 0.5, col.ticks = AXIS_COL, col.axis = AXIS_COL, cex.axis = 0.8)
        axis(2, lwd = 0, lwd.ticks = 0.5, col.ticks = AXIS_COL, col.axis = AXIS_COL, cex.axis = 0.8)
        box(lwd = 0.6, col = AXIS_COL)
        mtext("|t|", side = 1, line = 1.35, cex = 0.78)
        mtext("Density", side = 2, line = 1.9, cex = 0.78, las = 0)
        mtext(sprintf("%s, c = %s", ROWS[[j]]$short, sev_label(sev)), side = 3, line = 0.25, cex = 0.8)
      }
    }
    mtext(panel_title, side = 3, outer = TRUE, line = 0.2, cex = 0.85)
    dev.off()
  }
  invisible(paste0(stem, c(".pdf", ".png")))
}


## ---- 7. Run --------------------------------------------------------------------------------

# One cell per (S, row, panel, c); the engine's tau holds c. The seed depends on S, the panel
# and c but not on the row (common random numbers across rows).
make_cells_rounding <- function(rows, panels, sevs) {
  cells <- list()
  for (fi in seq_along(CFG$S_values)) {
    for (rw in rows) {
      for (k in seq_along(panels)) {
        for (j in seq_along(sevs)) {
          S   <- CFG$S_values[fi]
          pn  <- panels[[k]]
          sev <- sevs[j]
          cells[[length(cells) + 1]] <- list(
            S = S, row = rw$id, panel = pn$id, tau = sev, severity = sev,
            pool = pool_id(pn$id), column = rw$column,
            seed = CFG$seed + 1000000 * fi + 10000 * k + 100 * j,
            file = file.path(OUT_DIR, "cells",
                             sprintf("cell_S%d_%s_%s_c%s.rds", as.integer(S), safe_id(rw$id),
                                     safe_id(pn$id), sev_tag(sev))))
        }
      }
    }
  }
  cells
}

run_cells_rounding <- function(cells, cl) {
  acc <- NULL
  for (cell in cells) {
    s          <- summarise_cell(run_cell(cell, cl))
    s$severity <- cell$severity
    acc        <- rbind(acc, s)
    message(sprintf("S = %d  %-8s %-9s c = %-5s  T_inf %.3f   T_H %.3f   hybrid %.3f%s%s   (%.0f s)",
                    as.integer(cell$S), cell$row, cell$panel, sev_label(cell$severity),
                    s$T_inf, s$T_H, s$T_H_hybrid,
                    if (CFG$run_projection) sprintf("   projection %.3f", s$proj) else "",
                    if (CFG$run_ekw) sprintf("   EKW max %.3f", s$EKW_max) else "", s$secs))
  }
  acc
}

PROFILE <- brodeur_profile()
print_profile(PROFILE)

POOLS     <- list()
POOL_KEYS <- list()
for (k in seq_along(PANELS)) {
  id              <- pool_id(PANELS[[k]]$id)
  POOLS[[id]]     <- build_pool_rounding(k, PROFILE)
  POOL_KEYS[[id]] <- paste(attr(POOLS[[id]], "key"), "| zero se:", CFG$zero_se_rule,
                           "| kappa cap:", CFG$kappa_cap, "| rows v2")   # cells are redone if these change
}

ROUND_SUMMARY <- do.call(rbind, lapply(PANELS, function(pn) {
  data.frame(panel = pn$id, rounding_summary(POOLS[[pool_id(pn$id)]], SEV_GRID),
             stringsAsFactors = FALSE)
}))
write.csv(ROUND_SUMMARY, file.path(OUT_DIR, paste0("rounding_summary_", CFG$name, ".csv")),
          row.names = FALSE)
message("Rounding in the simulated tables (one_digit: SEs with one significant digit; below_37: ",
        "significand below 37; kappa_dropped: left out by the cap; and among the kept studies ",
        "mean_lambda: mean se / s_min; at_2: m/s exactly 2; mean_abs_err: mean |m/s - T*|; ",
        "crosses_1.96: m/s and T* on opposite sides of 1.96):")
print(ROUND_SUMMARY, digits = 3, row.names = FALSE)

if (CFG$run_selfchecks) rounding_selfcheck()

plot_tcurves(POOLS[[pool_id(PANELS[[1]]$id)]], file.path(OUT_DIR, paste0(CFG$name, "_tcurves")),
             sevs = unique(c(1, max(SEV_GRID))), panel_title = PANELS[[1]]$title)

PROJ <- if (CFG$run_projection) proj_constants(file.path(OUT_DIR, "proj_constants.rds")) else NULL
if (!is.null(PROJ)) {
  message(sprintf("Projection test: epsilon_U = %.5f per sqrt(n) (paper: 0.00063).", PROJ$epsilon_U_pre))
}

CELLS <- make_cells_rounding(ROWS, PANELS, SEV_GRID)
pilot_estimate(CELLS)

cl  <- start_cluster()
res <- tryCatch({
  if (!is.null(cl)) {
    parallel::clusterExport(cl, c("round_results", "kappa_ok", "row_t", "eligible_rows", "ELIG"),
                            envir = globalenv())
  }
  acc <- NULL
  for (S in CFG$S_values) {
    acc <- rbind(acc, run_cells_rounding(Filter(function(cc) cc$S == S, CELLS), cl))
    saveRDS(acc, file.path(OUT_DIR, paste0("results_", CFG$name, ".rds")))
    make_figures(acc, S)
    message(sprintf("Figures for S = %d written to %s", as.integer(S), OUT_DIR))
  }
  acc
}, finally = {
  if (!is.null(cl)) parallel::stopCluster(cl)
})

res <- merge(res, ROUND_SUMMARY, by = c("panel", "severity"), all.x = TRUE, sort = FALSE)
res <- res[order(res$S, match(res$row, CFG$rows), match(res$panel, PANEL_IDS), res$severity), ,
           drop = FALSE]
rownames(res) <- NULL
saveRDS(res, file.path(OUT_DIR, paste0("results_", CFG$name, ".rds")))
write.csv(res, file.path(OUT_DIR, paste0("results_", CFG$name, ".csv")), row.names = FALSE)
tau0_check(res)                           # c = 0: rows test the same honest meta-samples
message("Figures and results written to ", OUT_DIR)


## ---- 8. LaTeX tables -------------------------------------------------------------------------
## One table per S, laid out like the figure: a block of rows per treatment of the rounded
## numbers, c down each block, and a column group per panel with one column per test. Written
## to OUT_DIR/NAME_S1000.tex, ...; \input{} them in the paper, which needs \usepackage{booktabs}.
## To rebuild figures and tables from saved results without rerunning the simulations, run
## sections 0-6 and the function definitions of section 7, then
## res <- readRDS(file.path(OUT_DIR, paste0("results_", CFG$name, ".rds"))),
## make_figures(res, S) for each S, and this section.

TEST_LATEX <- c(proj = "Proj.", T_inf = "$T_\\infty$", T_H = "$T_H$",
                T_H_hybrid = "$T_H^{\\mathrm{hyb}}$", EKW_max = "EKW")

# Plain-text titles: escape the characters LaTeX treats as special (math in $...$ is kept)
latex_text <- function(x) {
  x <- gsub("([&%#_])", "\\\\\\1", x)
  x <- gsub("<=", "$\\\\le$", x)
  gsub(">=", "$\\\\ge$", x)
}

# Panel header from a pi0 spec (as pi0_label(), in LaTeX)
pi0_latex <- function(spec) {
  g <- function(x) format(signif(x, 4))
  s <- spec$shape
  if (identical(spec, PI0_ALIASES$kudrin_pi_hat)) return("$h \\sim \\hat{\\Pi}$")
  if (s == "point")   return(sprintf("$h = %s$", g(spec$h_center)))
  if (s == "normal")  return(sprintf("$h \\sim N(%s, %s^2)$", g(spec$h_center), g(spec$sigma_h)))
  if (s == "uniform") return(sprintf("$h \\sim U[%s, %s]$", g(spec$h_center - spec$sigma_h / 2),
                                     g(spec$h_center + spec$sigma_h / 2)))
  if (s == "chi2")    return(sprintf("$h \\sim \\chi^2_{%s}$", g(spec$h_center)))
  if (s == "poisson") return(sprintf("$h \\sim \\mathrm{Poisson}(%s)$", g(spec$h_center)))
  if (s == "gamma") {
    if (is.null(spec$gamma_scale)) {
      return(sprintf("$h \\sim \\Gamma(%s, \\mathrm{rate}\\ %s)$", g(spec$gamma_shape), g(spec$gamma_rate)))
    }
    return(sprintf("$h \\sim \\Gamma(%s, \\mathrm{scale}\\ %s)$", g(spec$gamma_shape), g(spec$gamma_scale)))
  }
  if (s == "double_normal") {
    return(sprintf("$h \\sim \\frac{1}{2} N(%s, %s^2) + \\frac{1}{2} N(%s, %s^2)$", g(spec$h_center),
                   g(spec$sigma_h), g(spec$h_center - 2), g(spec$sigma_h)))
  }
  if (s == "fitted") return("$h \\sim$ z-curve fit")
  "$h \\sim$ custom"
}

# The twopoint panel in LaTeX (pi0_latex() above has no case for it)
pi0_latex_rounding <- function(spec) {
  if (identical(spec$shape, "twopoint")) {
    g <- function(x) format(signif(x, 4))
    return(sprintf("$h \\in \\{%s, %s\\}$", g(spec$h_center), g(spec$h_center + 2)))
  }
  pi0_latex(spec)
}

TABLE_OPTS <- list(
  severities    = NULL,     # NULL: every c in the results; or a subset such as c(0, 1, 4, 16)
  digits        = 3,        # decimals of the rejection rates
  bold_oversize = TRUE,     # bold rates above the 95% simulation band around the level
  fragment      = FALSE     # TRUE: write only the tabular, without the table environment and notes
)

write_latex_tables_rounding <- function(res, opts = TABLE_OPTS) {
  tests <- intersect(c("proj", "T_inf", "T_H", "T_H_hybrid", "EKW_max"), names(res))
  kt    <- intersect(c("T_inf", "T_H", "T_H_hybrid"), tests)
  nt    <- length(tests)
  np    <- length(PANELS)
  ncols <- 1 + nt * np
  plab  <- vapply(seq_along(PANELS), function(k) {
    l <- CFG$panels[[k]]$latex
    if (is.null(l)) pi0_latex_rounding(PANELS[[k]]$spec) else l
  }, "")
  rlab  <- vapply(ROWS, function(rw) {
    if (!is.null(CFG$row_latex) && rw$id %in% names(CFG$row_latex)) return(CFG$row_latex[[rw$id]])
    if (!is.null(CFG$row_titles) && rw$id %in% names(CFG$row_titles)) return(latex_text(rw$title))
    if (!is.null(rw$latex)) rw$latex else latex_text(rw$title)
  }, "")
  head <- c(paste0(" & ", paste(sprintf("\\multicolumn{%d}{c}{%s}", nt, plab), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + (seq_len(np) - 1) * nt, 1 + seq_len(np) * nt),
                  collapse = " "),
            paste0("$c$ & ", paste(rep(TEST_LATEX[tests], np), collapse = " & "), " \\\\"))
  files <- character(0)
  
  for (S in sort(unique(res$S))) {
    rs    <- res[res$S == S, , drop = FALSE]
    sevs  <- sort(unique(rs$severity))
    if (!is.null(opts$severities)) sevs <- sevs[round(sevs, 8) %in% round(opts$severities, 8)]
    dag_k <- FALSE
    dag_p <- FALSE
    bold  <- FALSE
    body  <- character(0)
    for (i in seq_along(ROWS)) {
      if (i > 1) body <- c(body, "\\midrule")
      body <- c(body, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", ncols, rlab[i]))
      for (sev in sevs) {
        entries <- character(0)
        for (k in seq_along(PANELS)) {
          r <- rs[rs$row == ROWS[[i]]$id & rs$panel == PANELS[[k]]$id & abs(rs$severity - sev) < 1e-8, ,
                  drop = FALSE]
          if (nrow(r) == 0) {
            entries <- c(entries, rep("", nt))
            next
          }
          v  <- vapply(tests, function(tn) as.numeric(r[[tn]][1]), 0)
          tx <- ifelse(is.finite(v), formatC(v, format = "f", digits = opts$digits), "--")
          if (isTRUE(opts$bold_oversize)) {
            top   <- CFG$alpha + qnorm(0.975) * sqrt(CFG$alpha * (1 - CFG$alpha) / r$reps[1])
            b     <- is.finite(v) & v > top
            tx[b] <- paste0("\\textbf{", tx[b], "}")
            bold  <- bold || any(b)
          }
          if (length(kt) > 0 && isTRUE(r$kudrin_failed[1] > 0)) {
            j <- match(kt, tests)
            tx[j] <- paste0(tx[j], "$^\\dagger$")
            dag_k <- TRUE
          }
          if ("proj" %in% tests && isTRUE(r$proj_failed[1] > 0)) {
            j <- match("proj", tests)
            tx[j] <- paste0(tx[j], "$^\\ddagger$")
            dag_p <- TRUE
          }
          entries <- c(entries, tx)
        }
        lab <- sev_label(sev)
        if (abs(sev - 1) < 1e-8) lab <- paste0(lab, "$^\\ast$")
        body <- c(body, paste0(lab, " & ", paste(entries, collapse = " & "), " \\\\"))
      }
    }
    
    tab   <- c(sprintf("\\begin{tabular}{l%s}", strrep(strrep("c", nt), np)),
               "\\toprule", head, "\\midrule", body, "\\bottomrule", "\\end{tabular}")
    S_tex <- gsub(",", "{,}", format(S, big.mark = ",", scientific = FALSE), fixed = TRUE)
    notes <- paste0(
      sprintf("Rejection rates of %s\\%%-level tests over %s Monte Carlo repetitions; each meta-sample has ",
              format(100 * CFG$alpha), paste(format(unique(rs$reps), big.mark = ","), collapse = "/")),
      sprintf("$S = %s$ $t$-statistics and nothing is selectively reported, so every rejection is a false positive. ",
              S_tex),
      "Each $t$-statistic is computed from a coefficient $b$ and a standard error rounded to the same decimal place. ",
      "$c$: coarseness of rounding, the rounding unit relative to the standard error as a multiple of that in ",
      "Brodeur et al. (2020) ($^\\ast$: as in their data); $c = 0$: no rounding. ",
      "Blocks of rows: what the meta-analyst does with the rounded numbers; $s_{\\min} = s - b/2$ is the ",
      "smallest standard error printed as $s$ and $V$ is uniform on $\\pm 1/2$ unit of the coefficient's last digit. ",
      if (is.finite(CFG$kappa_cap)) sprintf("Every block uses the same results: those with $\\kappa = a/s_{\\min} \\le %s$. ",
                                            format(CFG$kappa_cap)) else "",
      "Column groups: distribution of ",
      "true effects $h$. Proj.: projection test; $T_\\infty$, $T_H$, $T_H^{\\mathrm{hyb}}$: Kudrin (2025)",
      if ("EKW_max" %in% tests) "; EKW: highest rejection rate among the EKW (2022) tests" else "", ".",
      if (bold) " Bold: above the 95\\% simulation band around the nominal level." else "",
      if (isTRUE(CFG$brodeur$synthetic)) " SYNTHETIC significands, not calibrated to Brodeur et al. (2020)." else "",
      if (dag_k) " $^\\dagger$Kudrin's tests could not be computed in some repetitions; the rate is among the others." else "",
      if (dag_p) " $^\\ddagger$The projection test could not be computed in some repetitions; the rate is among the others." else "")
    caption <- sprintf("Size when coefficients and standard errors are rounded, $S = %s$", S_tex)
    out <- if (isTRUE(opts$fragment)) tab else c(
      "\\begin{table}[htbp]", "\\centering",
      sprintf("\\caption{%s}", caption),
      sprintf("\\label{tab:%s_S%d}", CFG$name, as.integer(S)),
      "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
      tab,
      "\\par\\smallskip",
      sprintf("\\begin{minipage}{\\linewidth}\\footnotesize\\textit{Notes:} %s\\end{minipage}", notes),
      "\\end{table}")
    f <- file.path(OUT_DIR, sprintf("%s_S%d.tex", CFG$name, as.integer(S)))
    writeLines(out, f)
    files <- c(files, f)
  }
  invisible(files)
}

TABLE_FILES <- write_latex_tables_rounding(res)
message("LaTeX tables written: ", paste(basename(TABLE_FILES), collapse = ", "))