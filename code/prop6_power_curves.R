################################################################################
#  prop6_power_curves.R
#
#  Power curves under the DGP of Proposition 6 (Appendix C): smooth publication bias that is
#  detectable but departs from the honest t-curves at low frequencies, where Kudrin's (2025)
#  tests do not look.
#
#  DGP. True effects H ~ N(0, var_h); an honest t-statistic is T = H + Z with Z ~ N(0, 1); a
#  study is published with probability P(R | T = t) = 1 - exp(-t^2 / filter_scale). Proposition 6
#  has var_h = 2 and filter_scale = 12. A meta-sample has n studies with one t-statistic each: a
#  study went through the filter with probability tau (its t is a draw of T | R) and is honest
#  otherwise (its t is a draw of T). tau = 0 gives size; tau = 1 is the DGP of Proposition 6, whose
#  t-curve is zero at t = 0. Nothing is rounded or de-rounded: every test sees the exact t.
#
#  Which tau are detectable. With rho = P(not R) = sqrt(2/3) and c = tau rho / (1 - rho), the
#  reported t-curve is (1 + c) phi_3 - c phi_2, phi_v the N(0, v) density. Deconvolving the N(0, 1)
#  noise leaves the signed measure (1 + c) N(0, 2) - c N(0, 1), a probability distribution iff
#  c <= 1 + sqrt(2). So H0 holds exactly, and the selection is undetectable (Definition 1), iff
#  tau <= tau* = (1 + sqrt(2)) (1 - rho) / rho = 0.543 (TAU_STAR below, for any var_h and
#  filter_scale). The figure shades tau <= tau*, where no valid test has power above size.
#
#  Tests, on every meta-sample (each study is its own cluster):
#    Projection test  run_test_weight_fast() with the paper's settings: J = 40, J0 = 20, p = 2,
#                     centered at 1.96, a grid of 3,000 points on [-6.5, 6.5] plus 9999
#    Kudrin (2025)    kudrin_tests(): T_inf, T_H and T_H^hybrid (sinc kernel, kappa = 3)
#    EKW (2022)       the seven tests of ekw_rejections() in kudrin_figures_engine.R (binomial,
#                     density discontinuity, Fisher, LCM, CS1, CSUB, CS2B), keeping their p-values.
#                     CSUB is CoxShi(..., K = 1, B = 1), monotonicity with bounds: Tests.R's CoxShi()
#                     cannot run K = 0 (4J - 2 constraint rows against 2J right-hand sides), so the
#                     engine's CoxShi(..., 0, 1) errors and counts as a non-rejection every time
#    p-curve          the test for left skew by Stouffer's method, stouffer_p() copied from
#                     cochrane_power_experiment.R
#  The grey line rejects when the smallest of the eight EKW and p-curve p-values is at most alpha,
#  as the bottom row of Figure 3; it is not a size-controlled test.
#
#  Common random numbers: within a repetition every tau shares the honest draws, the published
#  draws and the uniforms that decide which studies went through the filter.
#
#  Outputs in OUT_DIR, with TAG = <CFG$name>_n<CFG$n>:
#    TAG.pdf/.png      the power curves, redrawn after every batch of repetitions
#    TAG_rates.csv     rejection rates by tau and test (every EKW test separately), the share of
#                      failed tests, and the population distance d_w of the t-curve from H0
#    TAG.tex           LaTeX table of the rejection rates at CFG$table_taus
#    reps_n<n>/        one file per repetition, so an interrupted run resumes
#
#  Run once with QUICK <- TRUE to check that everything works, then with FALSE.
################################################################################

rm(list = ls())

## ---- 0. Settings -------------------------------------------------------------------------

ROOT     <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
CODE_DIR <- file.path(ROOT, "code")
QUICK    <- FALSE   # TRUE: small run that checks the pipeline; FALSE: the settings below

CFG <- list(
  name          = "prop6_power",
  n             = 5000,                  # t-statistics per meta-sample (Kudrin's Fig. 3: 5,000)
  taus          = c(0,0.5,0.8, 0.85, 0.9, 0.95, 1),   # share through the filter;
  # extra points above 0.8, where the t-curve moves away from the honest set
  table_taus    = c(0, 0.5, 0.8, 0.9, 1),  # columns of the LaTeX table (must be among taus)
  
  # DGP of Proposition 6
  var_h         = 2,                     # H ~ N(0, var_h)
  filter_scale  = 12,                    # P(R | T = t) = 1 - exp(-t^2 / filter_scale)
  
  # Projection test: run_test_weight_fast(), with the paper's settings
  num_coeffs    = 41,                    # J = 40
  J0            = 20,                    # weights w_j = min(1, (1/p)^(j - J0))
  p_weight      = 2,
  shift         = 1.96,
  L             = 6.5,
  numgrid       = 3000,
  proj_boots    = 200,
  epsilon_U_pre = NULL,                  # NULL: compute_epsilons(L, numgrid) once and cache it
  # (slow); or 0.0006333757, the value in the paper
  
  # Kudrin (2025): kudrin_tests() with its defaults (sinc kernel, kappa = 3)
  kudrin_boots  = 199,                   # T_inf bootstrap draws
  
  # EKW (2022) tests, as in kudrin_figures_engine.R
  ekw_J         = 30,                    # Cox-Shi bins on [p_min, 0.15]
  
  # Computing
  reps          = 500,
  alpha         = 0.05,
  n_cores       = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  seed          = 20261103,
  resume        = TRUE
)
OUT_DIR <- file.path(ROOT, "simulations", CFG$name)
if (QUICK) {
  CFG$reps          <- 40
  CFG$proj_boots    <- 99
  CFG$kudrin_boots  <- 99
  CFG$epsilon_U_pre <- 0.0006333757
  OUT_DIR           <- paste0(OUT_DIR, "_quick")
}
TAG     <- sprintf("%s_n%d", CFG$name, as.integer(CFG$n))
REP_DIR <- file.path(OUT_DIR, sprintf("reps_n%d", as.integer(CFG$n)))
dir.create(REP_DIR, recursive = TRUE, showWarnings = FALSE)
stopifnot(CFG$n >= 100, CFG$var_h >= 0, CFG$filter_scale > 0, any(CFG$taus == 0),
          all(CFG$taus >= 0 & CFG$taus <= 1), !anyDuplicated(CFG$taus),
          all(round(CFG$table_taus, 8) %in% round(CFG$taus, 8)))


## ---- 1. Your code ------------------------------------------------------------------------

CODE_FILES <- file.path(CODE_DIR, c("functions_all_weighted.R", "Tests.R", "kudrin_tests.R"))
missing_files <- CODE_FILES[!file.exists(CODE_FILES)]
if (length(missing_files) > 0) stop("Not found: ", paste(missing_files, collapse = ", "))
for (f in CODE_FILES) source(f)
stopifnot(exists("run_test_weight_fast"), exists("create_basis_stable"), exists("compute_epsilons"),
          exists("project_simplex_nnls"), exists("kudrin_tests"), exists("CoxShi"), exists("LCM"))


## ---- 2. The DGP of Proposition 6 -----------------------------------------------------------

# f_T(t) exp(-t^2 / filter_scale) = P_NOT phi_{V_F}(t), with phi_v the N(0, v) density. So
# P(not R) = P_NOT (rho = sqrt(2/3) in Proposition 6) and f_{T|R} = (phi_{V_T} - P_NOT phi_{V_F}) / (1 - P_NOT)
V_T   <- 1 + CFG$var_h                                                 # Var(T): 3
V_F   <- V_T * (CFG$filter_scale / 2) / (V_T + CFG$filter_scale / 2)   # 2
P_NOT <- sqrt(V_F / V_T)
F_PUB <- function(x) (pnorm(x / sqrt(V_T)) - P_NOT * pnorm(x / sqrt(V_F))) / (1 - P_NOT)   # CDF of T | R

# Reported t-curve when a share tau went through the filter
f_tau <- function(t, tau) {
  f_hon <- dnorm(t, sd = sqrt(V_T))
  (1 - tau) * f_hon + tau * (f_hon - P_NOT * dnorm(t, sd = sqrt(V_F))) / (1 - P_NOT)
}

# Largest tau at which H0 holds: with c = tau P_NOT / (1 - P_NOT), the deconvolved t-curve is the
# signed measure (1 + c) N(0, V_T - 1) - c N(0, V_F - 1). The ratio of its two densities is smallest
# at h = 0, where it is r = sqrt((V_F - 1) / (V_T - 1)), so it is a probability distribution iff
# c <= r / (1 - r). If V_F <= 1 (or var_h = 0) there is nothing to deconvolve into and any tau > 0
# is detectable.
TAU_STAR <- if (CFG$var_h <= 0 || V_F <= 1) 0 else {
  r <- sqrt((V_F - 1) / (V_T - 1))
  min(1, r / (1 - r) * (1 - P_NOT) / P_NOT)
}

draw_honest <- function(n) rnorm(n, sd = sqrt(CFG$var_h)) + rnorm(n)    # T = H + Z

# n draws of T | R: studies T = H + Z, each published with probability 1 - exp(-T^2 / filter_scale),
# until n are published (the first n, in order)
draw_published <- function(n) {
  out <- numeric(0)
  while (length(out) < n) {
    m    <- ceiling(1.25 * (n - length(out)) / (1 - P_NOT)) + 100
    t    <- draw_honest(m)
    keep <- runif(m) < 1 - exp(-t^2 / CFG$filter_scale)
    out  <- c(out, t[keep])
  }
  out[seq_len(n)]
}

# Check the sampler against the exact distribution of T | R
local({
  set.seed(CFG$seed - 1)
  t0  <- draw_honest(1e6)
  pr  <- mean(runif(1e6) < 1 - exp(-t0^2 / CFG$filter_scale))
  x   <- draw_published(2e5)
  ks  <- suppressWarnings(ks.test(x, F_PUB))
  message(sprintf(paste("DGP check: P(R) = %.4f (exact %.4f); published draws against the exact CDF of T | R:",
                        "KS p = %.3f; P(|t| < 0.25) = %.5f (exact %.5f)."),
                  pr, 1 - P_NOT, ks$p.value, mean(abs(x) < 0.25), F_PUB(0.25) - F_PUB(-0.25)))
  if (ks$p.value < 1e-3) warning("The published draws do not match the exact distribution of T | R.")
})
message(sprintf("tau* = %.4f: H0 holds exactly, so the selection is undetectable, for tau <= tau*.", TAU_STAR))


## ---- 3. Tests ------------------------------------------------------------------------------

# Projection test, as in cochrane_power_experiment.R (PROJ: cached basis and epsilon_U, section 6)
proj_test <- function(t) {
  pr <- tryCatch(run_test_weight_fast(data.frame(t = t, title = as.character(seq_along(t)),
                                                 stringsAsFactors = FALSE),
                                      numcoeffs = CFG$num_coeffs, sigma_Y = 1, shift = CFG$shift,
                                      L = CFG$L, numgrid = CFG$numgrid, boots = CFG$proj_boots,
                                      U = PROJ$U, epsilon_U_pre = PROJ$epsilon_U_pre,
                                      p = CFG$p_weight, J0 = CFG$J0),
                 error = function(e) NULL)
  one <- function(v) if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)[1]
  if (is.null(pr)) return(c(proj = NA_real_, proj_stat = NA_real_, proj_breakdown = NA_real_))
  c(proj = one(pr$pval), proj_stat = one(pr$resid), proj_breakdown = one(pr$breakdown))
}

KUDRIN_NAMES <- c("T_inf", "T_H", "T_H_hybrid")
kudrin_p <- function(t) {
  kt <- tryCatch(kudrin_tests(t, boots = CFG$kudrin_boots, alpha = CFG$alpha), error = function(e) NULL)
  if (is.null(kt)) return(setNames(rep(NA_real_, 3), KUDRIN_NAMES))
  setNames(as.numeric(kt$pval[KUDRIN_NAMES]), KUDRIN_NAMES)
}

# The seven EKW tests on (0, 0.15], called as in ekw_rejections() of kudrin_figures_engine.R, except
# CSUB (see the header). A failed test (an error, or Cox-Shi's codes 888 and 999) gives NA.
EKW_NAMES <- c("EKW_Binomial", "EKW_Discontinuity", "EKW_Fisher", "EKW_LCM", "EKW_CS1", "EKW_CSUB",
               "EKW_CS2B")
ekw_p <- function(t) {
  ps   <- 2 * pnorm(-abs(t))
  pmn  <- min(ps[ps > 0])
  pmx  <- 0.15
  J    <- CFG$ekw_J
  safe <- function(expr) {
    v <- tryCatch(suppressWarnings(expr), error = function(e) NA_real_)
    if (length(v) != 1 || !is.numeric(v) || !is.finite(v) || v < 0 || v > 1) NA_real_ else as.numeric(v)
  }
  c(EKW_Binomial      = safe(Binomial(ps, 0.04, 0.05, "c")),
    EKW_Discontinuity = safe(Discontinuity_test(ps, 0.05)),
    EKW_Fisher        = safe(Fisher(ps, pmn, pmx)),
    EKW_LCM           = safe(LCM(ps, pmn, pmx)),
    EKW_CS1           = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 1, 0)),
    EKW_CSUB          = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 1, 1)),
    EKW_CS2B          = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 2, 1)))
}

# The p-curve test for left skew (Simonsohn, Nelson and Simmons 2014; Simonsohn's report, point 5):
# every p < .05 times 20, through the normal quantile function, combined by Stouffer's method,
# Z = sum(qnorm(20 p)) / sqrt(k); one-sided p-value 1 - Phi(Z) for left skew. Like p-curve, it
# treats the p-values as independent. qnorm(20 p) is computed on the log scale, so huge |t| do
# not underflow to -Inf. (Copied from cochrane_power_experiment.R.)
stouffer_p <- function(t) {
  a   <- abs(t)
  sig <- is.finite(a) & 2 * pnorm(-a) < 0.05
  k   <- sum(sig)
  if (k == 0) return(c(stouffer = NA_real_, stouffer_k = 0))
  z   <- qnorm(log(40) + pnorm(-a[sig], log.p = TRUE), log.p = TRUE)
  c(stouffer = pnorm(sum(z) / sqrt(k), lower.tail = FALSE), stouffer_k = k)
}

# Every test on one meta-sample; min_p is the smallest of the eight EKW and p-curve p-values
run_tests <- function(t) {
  e  <- ekw_p(t)
  s  <- stouffer_p(t)
  pv <- c(e, s[["stouffer"]])
  c(proj_test(t), kudrin_p(t), e, s, min_p = if (all(is.na(pv))) NA_real_ else min(pv, na.rm = TRUE))
}


## ---- 4. One repetition ---------------------------------------------------------------------

# Repetition r: n honest draws, n published draws and n uniforms, shared by every tau
one_rep <- function(r) {
  set.seed(CFG$seed + 1000 * r)
  t_hon <- draw_honest(CFG$n)
  t_pub <- draw_published(CFG$n)
  u     <- runif(CFG$n)
  out   <- lapply(CFG$taus, function(tau) {
    through <- u < tau
    data.frame(rep = r, tau = tau, filtered = mean(through),
               as.list(run_tests(ifelse(through, t_pub, t_hon))), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

rep_file <- function(r) file.path(REP_DIR, sprintf("rep_%05d.rds", as.integer(r)))

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


## ---- 5. Rejection rates, figure and table ---------------------------------------------------

TESTS   <- c("proj", KUDRIN_NAMES, EKW_NAMES, "stouffer", "min_p")
NA_IS_0 <- c(EKW_NAMES, "stouffer", "min_p")   # failed EKW test, or no p < .05 for p-curve: no rejection

# Rejection rates by tau. A failed EKW test counts as a non-rejection (as in kudrin_figures_engine.R);
# the projection and Kudrin tests' rates are among the repetitions in which they ran.
summarise_rates <- function(raw) {
  do.call(rbind, lapply(sort(unique(raw$tau)), function(tau) {
    d <- raw[abs(raw$tau - tau) < 1e-9, , drop = FALSE]
    s <- data.frame(tau = tau, reps = nrow(d), filtered = mean(d$filtered))
    for (tn in TESTS) {
      p <- d[[tn]]
      s[[tn]] <- if (tn %in% NA_IS_0) mean(!is.na(p) & p <= CFG$alpha) else mean(p <= CFG$alpha, na.rm = TRUE)
    }
    s$d_w_pop       <- D_POP[which.min(abs(CFG$taus - tau))]
    s$proj_stat     <- mean(d$proj_stat, na.rm = TRUE)
    s$proj_failed   <- mean(is.na(d$proj))
    s$kudrin_failed <- mean(is.na(d$T_inf))
    for (e in EKW_NAMES) s[[paste0(e, "_failed")]] <- mean(is.na(d[[e]]))
    s
  }))
}

# Style of kudrin_figures_engine.R; the grey line is the smallest EKW / p-curve p-value here
STYLE <- list(
  sig        = list(col = "#000000", lty = "dotted",   lwd = 0.8),
  min_p      = list(col = "#8C8C8C", lty = "longdash", lwd = 1.1, pch = 0,  cex = 0.6,  plwd = 0.9),
  T_inf      = list(col = "#0072B2", lty = "solid",    lwd = 1.1, pch = 2,  cex = 0.62, plwd = 0.9),
  T_H        = list(col = "#E69F00", lty = "dashed",   lwd = 1.1, pch = 3,  cex = 0.62, plwd = 0.9),
  T_H_hybrid = list(col = "#CC79A7", lty = "dotdash",  lwd = 1.1, pch = 5,  cex = 0.62, plwd = 0.9),
  proj       = list(col = "#000000", lty = "solid",    lwd = 1.6, pch = 16, cex = 0.55, plwd = 1.0)
)
LEG_LABEL <- list(sig        = sprintf("%s%% level", format(100 * CFG$alpha)),
                  min_p      = "Any EKW test or p-curve (smallest p)",
                  T_inf      = expression(italic(T)[infinity]),
                  T_H        = expression(italic(T)[italic(H)]),
                  T_H_hybrid = expression(italic(T)[italic(H)]^hybrid),
                  proj       = "Projection test")
LEG_ORDER <- c("sig", "min_p", "T_inf", "T_H", "T_H_hybrid", "proj")
AXIS_COL  <- "grey25"
GRID_COL  <- "grey88"
SHADE_COL <- "grey95"   # tau <= tau*: H0 holds
BAND_COL  <- "grey85"   # 95% Monte Carlo band around alpha
XLAB      <- expression("Share of studies through the filter ("*tau*"); at "*tau == 0*": size")

draw_series <- function(x, y, st) {
  ok <- is.finite(y)
  lines(x[ok], y[ok], col = st$col, lty = st$lty, lwd = st$lwd)
  if (!is.null(st$pch)) points(x[ok], y[ok], col = st$col, pch = st$pch, cex = st$cex, lwd = st$plwd, xpd = TRUE)
  invisible(NULL)
}

open_device <- function(file, ext, height, width) {
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

# The power curves of every test in one panel
plot_curves <- function(sm, stem) {
  sm   <- sm[order(sm$tau), , drop = FALSE]
  reps <- min(sm$reps)
  note <- c(sprintf("H ~ N(0, %s), T = H + Z; a study is published with probability 1 - exp(-t^2/%s); nothing is rounded.",
                    format(CFG$var_h), format(CFG$filter_scale)),
            sprintf("n = %s t-statistics per meta-sample; %d repetitions; dark band: 95%% Monte Carlo band around %s%%.",
                    format(CFG$n, big.mark = ","), as.integer(reps), format(100 * CFG$alpha)))
  ncol      <- 3
  line_in   <- 8 * 1.2 / 72                        # height of one margin line at pointsize 8
  leg_lines <- 0.9 + 1.1 * ceiling(length(LEG_ORDER) / ncol)
  top_lines <- 1.4 + 1.1 * length(note)
  height    <- 3.3 + (leg_lines + top_lines) * line_in
  ts        <- min(max(TAU_STAR, 0), 1)
  half      <- qnorm(0.975) * sqrt(CFG$alpha * (1 - CFG$alpha) / max(1, reps))
  for (ext in c("pdf", "png")) {
    open_device(paste0(stem, ".", ext), ext, height = height, width = 4.8)
    par(cex = 1, mar = c(2.8, 3.0, 1.5, 0.8), oma = c(leg_lines, 0, top_lines, 0),
        mgp = c(1.5, 0.35, 0), tcl = -0.22, las = 1, family = "sans")
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
    if (ts > 0) rect(0, 0, ts, 1, col = SHADE_COL, border = NA)
    rect(0, max(0, CFG$alpha - half), 1, CFG$alpha + half, col = BAND_COL, border = NA)
    abline(h = c(0.25, 0.5, 0.75), v = c(0.25, 0.5, 0.75), col = GRID_COL, lwd = 0.5)
    abline(h = CFG$alpha, col = STYLE$sig$col, lty = STYLE$sig$lty, lwd = STYLE$sig$lwd)
    if (ts > 0 && ts < 1) {
      abline(v = ts, col = "grey35", lty = "22", lwd = 0.8)
      text(ts + 0.012, 0.97, bquote(tau^"*" == .(sprintf("%.3f", ts))), adj = c(0, 1), cex = 0.8, col = "grey20")
    }
    if (ts >= 0.25) {
      mtext(expression(italic(H)[0]~"holds: undetectable"), side = 3, line = 0.2, at = ts / 2,
            cex = 0.72, col = "grey25")
    }
    if (ts <= 0.75) {
      mtext(expression(italic(H)[0]~"fails: detectable"), side = 3, line = 0.2, at = (1 + ts) / 2,
            cex = 0.72, col = "grey25")
    }
    at   <- c(0, 0.25, 0.5, 0.75, 1)
    labs <- c("0", ".25", ".5", ".75", "1")
    for (side in 1:2) {
      axis(side, at = at, labels = labs, lwd = 0, lwd.ticks = 0.5, col.ticks = AXIS_COL,
           col.axis = AXIS_COL, cex.axis = 0.8)
    }
    box(lwd = 0.6, col = AXIS_COL)
    mtext(XLAB, side = 1, line = 1.5, cex = 0.78)
    mtext("Rejection rate", side = 2, line = 1.9, cex = 0.78, las = 0)
    for (k in c("min_p", "T_inf", "T_H", "T_H_hybrid", "proj")) draw_series(sm$tau, sm[[k]], STYLE[[k]])
    mtext("Smooth publication bias (Proposition 6)", side = 3, outer = TRUE,
          line = 0.25 + 1.1 * length(note), font = 2, cex = 0.95)
    for (j in seq_along(note)) {
      mtext(note[j], side = 3, outer = TRUE, line = 0.25 + 1.1 * (length(note) - j), cex = 0.68)
    }
    draw_legend_strip(LEG_ORDER, ncol)
    dev.off()
  }
  invisible(paste0(stem, c(".pdf", ".png")))
}

# LaTeX table: one row per test, one column per tau in CFG$table_taus
TEST_TEX <- c(proj = "Projection test", T_inf = "$T_\\infty$", T_H = "$T_H$",
              T_H_hybrid = "$T_H^{\\mathrm{hyb}}$", EKW_Binomial = "EKW: binomial",
              EKW_Discontinuity = "EKW: density discontinuity", EKW_Fisher = "EKW: Fisher",
              EKW_LCM = "EKW: LCM", EKW_CS1 = "EKW: CS1", EKW_CSUB = "EKW: CSUB", EKW_CS2B = "EKW: CS2B",
              stouffer = "$p$-curve (Stouffer)", min_p = "Any EKW test or $p$-curve (smallest $p$-value)")

write_tex <- function(sm, file) {
  rows <- vapply(CFG$table_taus, function(x) which(abs(sm$tau - x) < 1e-9)[1], 0L)
  rows <- rows[!is.na(rows)]
  if (length(rows) == 0) return(invisible(NULL))
  taus <- sm$tau[rows]
  reps <- min(sm$reps[rows])
  half <- qnorm(0.975) * sqrt(CFG$alpha * (1 - CFG$alpha) / reps)
  cell <- function(v, tau, sized) {
    if (!is.finite(v)) return("--")
    s <- sprintf("%.3f", v)
    if (sized && tau <= TAU_STAR + 1e-12 && v > CFG$alpha + half) s <- paste0("\\textbf{", s, "}")
    s
  }
  body <- character(0)
  for (tn in TESTS) {
    if (tn %in% c("T_inf", "EKW_Binomial", "min_p")) body <- c(body, "\\midrule")
    v    <- vapply(seq_along(rows), function(k) cell(sm[[tn]][rows[k]], taus[k], tn != "min_p"), "")
    body <- c(body, paste0(TEST_TEX[[tn]], " & ", paste(v, collapse = " & "), " \\\\"))
  }
  head  <- paste0("Test & ", paste(sprintf("$\\tau = %s$", vapply(taus, format, "")), collapse = " & "), " \\\\")
  notes <- paste0(
    sprintf("Rejection rates of %s\\%%-level tests over %d repetitions. ", format(100 * CFG$alpha), as.integer(reps)),
    sprintf(paste0("Each meta-sample has $n = %s$ $t$-statistics. A share $\\tau$ of the studies went through the ",
                   "publication filter of Proposition 6: true effects $H \\sim N(0, %s)$, $T = H + Z$, and a study is ",
                   "published with probability $1 - \\exp(-t^2/%s)$; the other studies are honest draws of $T$. "),
            gsub(",", "{,}", format(CFG$n, big.mark = ","), fixed = TRUE), format(CFG$var_h), format(CFG$filter_scale)),
    sprintf(paste0("For $\\tau \\le \\tau^* = %.3f$ the reported $t$-curve satisfies $H_0$ exactly, so the selection is ",
                   "undetectable there; at $\\tau = 1$ it is the $t$-curve of Proposition 6. Nothing is rounded. "), TAU_STAR),
    sprintf("Projection test: $J = %d$, $J_0 = %s$, centered at %s, %d bootstrap repetitions. ",
            as.integer(CFG$num_coeffs - 1), format(CFG$J0), format(CFG$shift), as.integer(CFG$proj_boots)),
    sprintf(paste0("Kudrin (2025): $T_\\infty$ with %d bootstrap repetitions, $T_H$ and $T_H^{\\mathrm{hyb}}$, with the ",
                   "sinc kernel and frequency cutoff 3. "), as.integer(CFG$kudrin_boots)),
    sprintf(paste0("EKW: the seven tests of Elliott, Kudrin and W\\\"uthrich (2022) on the $p$-values in $(0, 0.15]$ ",
                   "(Cox--Shi tests with %d bins); a test that fails counts as a non-rejection. "), as.integer(CFG$ekw_J)),
    "$p$-curve: the test for left skew of Simonsohn, Nelson and Simmons (2014), by Stouffer's method. ",
    sprintf("The last row rejects when the smallest of these eight $p$-values is at most %s; it is not a size-controlled test. ",
            format(CFG$alpha)),
    sprintf("Bold: above %.3f, the upper end of the 95\\%% Monte Carlo band around %s\\%%, at a $\\tau$ where $H_0$ holds.",
            CFG$alpha + half, format(100 * CFG$alpha)))
  writeLines(c("\\begin{table}[htbp]", "\\centering",
               "\\caption{Rejection rates under the DGP of Proposition 6}",
               sprintf("\\label{tab:%s}", TAG), "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
               sprintf("\\begin{tabular}{l%s}", strrep("c", length(rows))), "\\toprule", head, "\\midrule",
               body, "\\bottomrule", "\\end{tabular}", "\\par\\smallskip",
               sprintf("\\begin{minipage}{\\linewidth}\\footnotesize\\textit{Notes:} %s\\end{minipage}", notes),
               "\\end{table}"), file)
  invisible(file)
}


## ---- 6. Run --------------------------------------------------------------------------------

# Projection-test constants (basis and epsilon_U), cached as in cochrane_power_experiment.R
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
message(sprintf("Projection test: epsilon_U = %.7f per sqrt(n) (paper: 0.0006333757).", PROJ$epsilon_U_pre))

# Population distance d_w of the reported t-curve from the honest set at each tau, with the test's
# weights, shift and grid: theta_j = E[phi(T - shift) psi_j(T - shift)] by quadrature (the t-curve is
# symmetric, so symmetrizing changes nothing), then the projection of run_test_weight_fast(). The
# test can reject in large samples only where d_w exceeds epsilon_U (Proposition 2).
D_POP <- local({
  h  <- 0.001
  tt <- seq(-20, 20, by = h)
  Ph <- phi_psi_stable(tt - CFG$shift, CFG$num_coeffs) * h
  w  <- create_weights(CFG$num_coeffs, p = CFG$p_weight, J0 = CFG$J0)
  Uw <- PROJ$U * w
  vapply(CFG$taus, function(tau) {
    th <- colSums(f_tau(tt, tau) * Ph)
    project_simplex_nnls(w * th, rbind(Uw, 1e5), Uw, 1e5)$residual
  }, 0)
})
message(sprintf("Population distance d_w x 1000 by tau (epsilon_U x 1000 = %.2f): %s", 1e3 * PROJ$epsilon_U_pre,
                paste(sprintf("%s: %.2f", vapply(CFG$taus, format, ""), 1e3 * D_POP), collapse = ", ")))

RUN_KEY <- paste("prop6", CFG$n, paste(CFG$taus, collapse = ","), CFG$var_h, CFG$filter_scale,
                 CFG$num_coeffs, CFG$J0, CFG$p_weight, CFG$shift, CFG$L, CFG$numgrid, CFG$proj_boots,
                 PROJ$epsilon_U_pre, CFG$kudrin_boots, CFG$ekw_J, CFG$seed)

# Pilot: every test on one meta-sample at tau = 1, for the run time and a first look
set.seed(CFG$seed)
pilot_t <- draw_published(CFG$n)
invisible(kudrin_p(pilot_t))                       # warm-up: builds kudrin_tests()' cached bounds
secs <- system.time(pilot <- run_tests(pilot_t))[["elapsed"]]
message("Pilot meta-sample at tau = 1: ", paste(sprintf("%s %.3g", names(pilot), pilot), collapse = ", "))
message(sprintf("One battery of tests takes %.1f s; %d per repetition; %d repetitions on %d core(s): about %.1f hours.",
                secs, length(CFG$taus), as.integer(CFG$reps), as.integer(CFG$n_cores),
                secs * length(CFG$taus) * CFG$reps / CFG$n_cores / 3600))

cl <- NULL
if (CFG$n_cores > 1) {
  cl <- parallel::makeCluster(CFG$n_cores)
  parallel::clusterExport(cl, "CODE_FILES", envir = globalenv())
  invisible(parallel::clusterEvalQ(cl, { for (f in CODE_FILES) source(f); NULL }))
  parallel::clusterExport(cl, c("CFG", "PROJ", "RUN_KEY", "REP_DIR", "P_NOT", "KUDRIN_NAMES",
                                "draw_honest", "draw_published", "proj_test", "kudrin_p", "ekw_p",
                                "stouffer_p", "run_tests", "one_rep", "rep_file", "one_rep_cached"),
                          envir = globalenv())
}
RAW <- tryCatch({
  acc   <- NULL
  batch <- 4 * max(1L, CFG$n_cores)
  for (b0 in seq(1, CFG$reps, by = batch)) {
    rr  <- b0:min(CFG$reps, b0 + batch - 1)
    t0  <- Sys.time()
    out <- if (is.null(cl)) lapply(rr, one_rep_cached) else parallel::parLapply(cl, rr, one_rep_cached)
    acc <- rbind(acc, do.call(rbind, out))
    sm  <- summarise_rates(acc)
    tryCatch(plot_curves(sm, file.path(OUT_DIR, TAG)),
             error = function(e) warning("Figure not drawn: ", conditionMessage(e)))
    top <- sm[nrow(sm), ]
    message(sprintf(paste("Repetitions %d-%d of %d done (%.1f min). At tau = %s: projection %.3f, T_inf %.3f,",
                          "T_H %.3f, hybrid %.3f, EKW or p-curve %.3f."),
                    min(rr), max(rr), as.integer(CFG$reps), as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    format(top$tau), top$proj, top$T_inf, top$T_H, top$T_H_hybrid, top$min_p))
  }
  acc
}, finally = {
  if (!is.null(cl)) parallel::stopCluster(cl)
})

SM <- summarise_rates(RAW)
write.csv(SM, file.path(OUT_DIR, paste0(TAG, "_rates.csv")), row.names = FALSE)
plot_curves(SM, file.path(OUT_DIR, TAG))
write_tex(SM, file.path(OUT_DIR, paste0(TAG, ".tex")))
print(SM[, c("tau", "reps", "filtered", "d_w_pop", TESTS)], digits = 3, row.names = FALSE)

# Where H0 holds (tau <= tau*), a valid test should reject at most about alpha
half <- qnorm(0.975) * sqrt(CFG$alpha * (1 - CFG$alpha) / CFG$reps)
h0   <- SM[SM$tau <= TAU_STAR + 1e-12, , drop = FALSE]
over <- Filter(function(tn) any(h0[[tn]] > CFG$alpha + half, na.rm = TRUE), setdiff(TESTS, "min_p"))
message(if (length(over) == 0) {
  sprintf("At every tau <= tau* (H0 holds), every test rejects at most %.3f.", CFG$alpha + half)
} else {
  sprintf("Above %.3f at some tau <= tau* (H0 holds): %s.", CFG$alpha + half, paste(over, collapse = ", "))
})
message("Outputs in ", OUT_DIR)

## ---- 7. The same figure without the grey field for tau <= tau* ---------------------------------
SHADE_COL <- NA                                               # rect() with col = NA draws nothing
plot_curves(SM, file.path(OUT_DIR, paste0(TAG, "_noshade")))
SHADE_COL <- "grey95"