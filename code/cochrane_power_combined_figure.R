#### =====================================================================================================
#### cochrane_power_combined_figure.R
#### One figure from a finished run of cochrane_power_experiment.R, without rerunning it: four rows (the
#### projection test, Kudrin's T_H and T_inf, and the smallest p-value of the six EKW tests and the p-curve
#### test), each under a header naming its test, by two columns (switching to the most significant outcome;
#### threshold switching), with one legend and the same axes in all eight panels. The panels are drawn as
#### plot_power() draws them. The notes are not in the figure: they go to a text file, for the caption.
#### Reads, from the run's folder:
####   <name>_rates.csv     the rejection rates (projection test, T_H, T_inf)
####   reps/rep_*.rds       the p-values of every repetition, for the EKW-or-p-curve row; if they are missing
####                        and the p-curve test never rejects in the run, that row equals the EKW row exactly
####   <name>_sample.txt    the trial counts for the notes (optional)
#### Writes, in the same folder: <name>_combined.pdf, <name>_combined_caption.txt (the notes),
#### <name>_combined_rates.csv (the plotted rates) and <name>_combined_table.tex (every test at tau = 0 and
#### tau = 1, with full notes, for the appendix).
#### =====================================================================================================
rm(list = ls())

#### 0. Settings -------------------------------------------------------------------------------------------
ROOT     <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
RUN_NAME <- "cochrane_switching_power_desktop"          #the run's CFG$name: the prefix of its files
RUN_DIR  <- file.path(ROOT, "simulations", RUN_NAME)    #the run's OUT_DIR (a QUICK run adds "_quick")
ALPHA    <- 0.05                                        #the run's CFG$alpha
CV       <- 1.96                                        #the run's CFG$cv, for threshold switching
ALLOW_C  <- 0.13                                        #the run's rounding allowance constant (J0 = 20)
PANEL_W  <- 2.85                                        #width and height of each panel's plotting area, inches
PANEL_H  <- 1.45
HEADER_H <- 0.26                                        #height of the header above each row, inches
HEADER_FILL <- "grey85"                                 #the band behind each row's header; NA for none
OUT_STEM <- file.path(RUN_DIR, paste0(RUN_NAME, "_combined"))
TAB_LABEL     <- "tab:cochrane_power_all_tests"    #the appendix table's LaTeX label
TAB_MARK_SIZE <- TRUE                              #bold the sizes above the 95% Monte Carlo band
PROJ_BOOTS    <- 200                               #the run's CFG$proj_boots and CFG$kudrin_boots (not saved by the run)
KUDRIN_BOOTS  <- 199
RUN_FACTS <- c(n_min = "100", n_meta = "1,500", n_elig = "3,282", min_out = "2", max_out = "8", kappa_cap = "0.5")
#the run's sample, used only where <name>_sample.txt lacks it

# The rows, top to bottom: the column of rejection rates each shows, and the header above it
ROWS <- list(
  proj    = list(col = "proj",           label = "Projection test"),
  T_H     = list(col = "T_H",            label = expression(bold("Kudrin's" ~ T[H]))),
  T_inf   = list(col = "T_inf",          label = expression(bold("Kudrin's" ~ T[infinity]))),
  ekw_min = list(col = "EKW_pcurve_min", label = "Any EKW test or the p-curve test (smallest p-value)"))

METHOD_TITLE <- c(max       = "Outcome switching: most significant outcome",
                  threshold = sprintf("Threshold: another outcome if |z| <= %s", format(CV)))
OPT_STYLE <- list(                                      #as in cochrane_power_experiment.R
  none       = list(col = "#8C8C8C", lty = "longdash", lwd = 1.3, pch = 4,  label = "No de-rounding"),
  uniform    = list(col = "#E69F00", lty = "dashed",   lwd = 1.3, pch = 2,  label = "Uniform de-rounding"),
  smin       = list(col = "#0072B2", lty = "solid",    lwd = 1.3, pch = 1,  label = "Smallest consistent SE"),
  smin_noise = list(col = "#000000", lty = "solid",    lwd = 1.8, pch = 16, label = "Smallest consistent SE + noise"))
XLAB <- expression("Share of trials that switch (" * tau * "); at " * tau == 0 * ": size")

#### 1. The rejection rates --------------------------------------------------------------------------------
rates_file <- file.path(RUN_DIR, paste0(RUN_NAME, "_rates.csv"))
if (!file.exists(rates_file)) stop("No ", rates_file, ": set RUN_NAME and RUN_DIR to a finished run", call. = FALSE)
SM <- read.csv(rates_file, stringsAsFactors = FALSE)
need <- c("method", "tau", "option", "reps", "proj", "T_H", "T_inf", "EKW_min", "stouffer")
if (length(setdiff(need, names(SM)))) stop(basename(rates_file), " lacks the columns ", paste(setdiff(need, names(SM)), collapse = ", "),
                                           " (was the run made with run_kudrin, run_ekw and run_stouffer?)", call. = FALSE)
METHODS <- intersect(c("max", "threshold"), unique(SM$method))
OPTIONS <- intersect(names(OPT_STYLE), unique(SM$option))
n_reps  <- as.integer(max(SM$reps))
cell_of <- function(d) paste(d$method, round(d$tau, 6), d$option, sep = "|")

# The EKW-or-p-curve row: in each repetition, the smallest p-value of the six EKW tests (EKW_min, already the
# smallest of the six) and the p-curve test; a rejection when it is at most ALPHA. As summarise_rates(), the
# tau = 0 ("honest") repetitions count for both methods, and missing p-values are left out.
rep_files <- file.path(RUN_DIR, "reps", sprintf("rep_%05d.rds", seq_len(n_reps)))
have <- file.exists(rep_files)
RAW  <- NULL
if (any(have)) {
  reps <- lapply(rep_files[have], function(f) tryCatch(readRDS(f), error = function(e) NULL))
  ok   <- vapply(reps, function(z) is.list(z) && !is.null(z$key) && is.data.frame(z$res), TRUE)
  if (any(ok)) {
    keys <- vapply(reps[ok], function(z) paste(z$key, collapse = " "), "")
    main <- names(sort(table(keys), decreasing = TRUE))[1]   #the run's key: the one most repetitions carry
    use  <- ok
    use[ok] <- keys == main
    if (sum(use) < n_reps) warning("Only ", sum(use), " of the ", n_reps, " repetitions in ", basename(rates_file),
                                   " are in the reps folder with the run's key", call. = FALSE)
    RAW <- do.call(rbind, lapply(reps[use], `[[`, "res"))
  }
}
if (!is.null(RAW)) {
  h    <- RAW[RAW$method == "honest", , drop = FALSE]
  RAW2 <- rbind(RAW[RAW$method != "honest", , drop = FALSE],
                do.call(rbind, lapply(METHODS, function(m) { z <- h; z$method <- m; z })))
  cell_raw <- cell_of(RAW2)
  rate <- function(p) tapply(p <= ALPHA, cell_raw, mean, na.rm = TRUE)
  for (cn in c("EKW_min", "stouffer")) {                       #do the repetitions reproduce the rates file?
    r <- as.numeric(rate(RAW2[[cn]])[cell_of(SM)])
    if (!isTRUE(all.equal(r, as.numeric(SM[[cn]]), tolerance = 1e-9)))
      warning("The repetitions do not reproduce the ", cn, " rates of ", basename(rates_file), " (stale reps?)", call. = FALSE)
  }
  smallest <- suppressWarnings(pmin(RAW2$EKW_min, RAW2$stouffer, na.rm = TRUE))
  SM$EKW_pcurve_min <- as.numeric(rate(smallest)[cell_of(SM)])
  message("EKW-or-p-curve row computed from ", length(unique(RAW$rep)), " repetitions")
} else if (isTRUE(all(SM$stouffer == 0, na.rm = TRUE))) {
  SM$EKW_pcurve_min <- SM$EKW_min                               #exact: the p-curve test never rejects in this run
  message("No reps folder; the p-curve test never rejects in this run, so the EKW-or-p-curve row is the EKW row")
} else {
  stop("The EKW-or-p-curve row needs the repetitions in ", file.path(RUN_DIR, "reps"), call. = FALSE)
}

#### 2. The caption: the notes, written to a text file instead of into the figure ------------------------
smp_file <- file.path(RUN_DIR, paste0(RUN_NAME, "_sample.txt"))
smp  <- if (file.exists(smp_file)) readLines(smp_file, warn = FALSE) else character(0)
grab <- function(pattern, i = 1) {                              #a number from the run's sample summary
  l <- grep(pattern, smp, value = TRUE)
  if (length(l)) sub(pattern, paste0("\\", i), l[1]) else NA_character_
}
pat_trials <- "^Trials: ([0-9,]+) randomized trials with at least ([0-9]+) continuous efficacy outcomes of at least ([0-9]+) participants.*$"
facts <- c(n_elig    = grab(pat_trials, 1), min_out = grab(pat_trials, 2), n_min = grab(pat_trials, 3),
           n_meta    = grab("^Meta-samples: ([0-9,]+) trials each.*$"),
           max_out   = grab("^Trials: .*\\(their first ([0-9]+) in review order\\).*$"),
           kappa_cap = grab("^Rounding: kappa = 2a / SE_min; the cap ([0-9.]+) leaves out.*$"))
miss <- is.na(facts)
if (any(miss)) {
  facts[miss] <- RUN_FACTS[names(facts)[miss]]
  message("Not found in ", basename(smp_file), ": ", paste(names(facts)[miss], collapse = ", "), "; RUN_FACTS used for them")
}
note1 <- sprintf("Cochrane: continuous efficacy outcomes of randomized trials (n >= %s); meta-samples of %s of %s trials with >= %s outcomes; %d repetitions; shaded: 95%% band around %s%%.",
                 facts[["n_min"]], facts[["n_meta"]], facts[["n_elig"]], facts[["min_out"]], n_reps, format(100 * ALPHA))
notes <- c(
  note1,
  sprintf("A share tau of trials reports its most significant outcome (max) or switches when |z| <= %s (threshold); the rest report a random one.", format(CV)),
  sprintf("Under s_min and s_min + noise the projection test rejects only if B-hat > rounding allowance %s x mean(kappa^6) of the reported outcomes; otherwise if p <= %s.",
          format(ALLOW_C), format(ALPHA)),
  sprintf("Any EKW test or the p-curve test: a rejection when the smallest p-value of the six EKW tests and the p-curve test for left skew is at most %s.",
          format(ALPHA)))
writeLines(notes, paste0(OUT_STEM, "_caption.txt"))

#### 3. The figure -----------------------------------------------------------------------------------------
LINE_IN <- 8 / 72 * 1.2                                         #one margin line at pointsize 8, in inches
PMAR    <- c(1.3, 3.0, 0.4, 0.5)                                #panel margins, first column (lines)
PMAR2   <- c(1.3, 0.8, 0.4, 0.5)                                #other columns: no tick labels or axis title
nr      <- length(ROWS)
nc      <- length(METHODS)
mar_of  <- function(j) if (j == 1) PMAR else PMAR2
# Cell sizes in inches, and the cells, numbered in the order they are drawn: the column titles, then for each
# row its header (across all columns) and its panels, then the axis titles and the legend (across all columns)
W <- vapply(seq_len(nc), function(j) PANEL_W + (mar_of(j)[2] + mar_of(j)[4]) * LINE_IN, 0)
H <- c(0.24, rep(c(HEADER_H, PANEL_H + (PMAR[1] + PMAR[3]) * LINE_IN), nr), 0.30, 0.40)
M  <- matrix(seq_len(nc), nrow = 1)
id <- nc
for (r in seq_len(nr)) {
  M  <- rbind(M, rep(id + 1, nc), id + 1 + seq_len(nc))
  id <- id + 1 + nc
}
M <- rbind(M, id + seq_len(nc), rep(id + nc + 1, nc))

draw_combined <- function() {
  layout(M, widths = lcm(2.54 * W), heights = lcm(2.54 * H))
  par(mgp = c(1.9, 0.4, 0), tcl = -0.25, las = 1, cex = 1)
  cell  <- function(mar) { par(mar = mar, cex = 1); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i") }
  ticks <- seq(0, 1, 0.25)
  tlab  <- c("0", ".25", ".5", ".75", "1")
  for (j in seq_len(nc)) {                                      #the column titles
    cell(c(0, mar_of(j)[2], 0, mar_of(j)[4]))
    text(0.5, 0.4, METHOD_TITLE[[METHODS[j]]], font = 2, cex = 0.95)
  }
  for (r in seq_len(nr)) {
    cell(c(0, PMAR[2], 0, mar_of(nc)[4]))                       #the row's header, spanning its panels
    if (!is.na(HEADER_FILL)) rect(0, 0.08, 1, 0.92, col = HEADER_FILL, border = NA)
    text(0.5, 0.5, ROWS[[r]]$label, font = 2, cex = 1.05)
    for (j in seq_len(nc)) {                                    #the panels, drawn as plot_power() draws them
      d <- SM[SM$method == METHODS[j], , drop = FALSE]
      par(mar = mar_of(j), cex = 1)
      plot.new()
      plot.window(xlim = c(-0.03, 1.03), ylim = c(0, 1), xaxs = "i", yaxs = "i")
      half <- qnorm(0.975) * sqrt(ALPHA * (1 - ALPHA) / max(1, min(d$reps)))
      rect(-0.03, ALPHA - half, 1.03, ALPHA + half, col = "grey92", border = NA)
      abline(h = c(0.25, 0.5, 0.75), v = c(0.25, 0.5, 0.75), col = "grey88", lwd = 0.5)
      abline(h = ALPHA, lty = 3)
      for (o in OPTIONS) {
        z  <- d[d$option == o, , drop = FALSE]
        z  <- z[order(z$tau), , drop = FALSE]
        st <- OPT_STYLE[[o]]
        y  <- z[[ROWS[[r]]$col]]
        ok <- is.finite(y)
        lines(z$tau[ok], y[ok], col = st$col, lty = st$lty, lwd = st$lwd)
        points(z$tau[ok], y[ok], col = st$col, pch = st$pch, cex = 0.75, xpd = TRUE)
      }
      axis(1, at = ticks, labels = if (r == nr) tlab else FALSE)    #one scale: tick labels on the outer panels
      axis(2, at = ticks, labels = if (j == 1) tlab else FALSE)
      box()
      if (j == 1) mtext("Rejection rate", side = 2, line = 2.1, las = 0, cex = 0.8)
    }
  }
  for (j in seq_len(nc)) {                                      #the axis titles
    cell(c(0, mar_of(j)[2], 0, mar_of(j)[4]))
    text(0.5, 0.6, XLAB, cex = 0.85)
  }
  cell(c(0, 0, 0, 0))                                           #the legend
  sty <- OPT_STYLE[OPTIONS]
  legend("center", legend = vapply(sty, function(s) s$label, ""), ncol = 2, bty = "n",
         col = vapply(sty, function(s) s$col, ""), lty = vapply(sty, function(s) s$lty, ""),
         lwd = vapply(sty, function(s) s$lwd, 0), pch = vapply(sty, function(s) s$pch, 0),
         cex = 0.95, seg.len = 2.6)
}

grDevices::pdf(paste0(OUT_STEM, ".pdf"), width = sum(W) + 0.02, height = sum(H) + 0.02, pointsize = 8)  #a little slack
draw_combined()
grDevices::dev.off()

#### 4. The rates behind the figure ------------------------------------------------------------------------
shown <- SM[order(match(SM$method, METHODS), match(SM$option, OPTIONS), SM$tau),
            c("method", "option", "tau", "reps", vapply(ROWS, function(x) x$col, ""))]
write.csv(shown, paste0(OUT_STEM, "_rates.csv"), row.names = FALSE)
print(shown[shown$tau %in% c(0, 0.5, 1), ], digits = 3, row.names = FALSE)

#### 5. The appendix table: every test, with no switching (tau = 0) and when every trial switches (tau = 1) ---
# Rows: every test, in groups. Columns: size, then tau = 1 under each kind of switching, each under the four
# treatments of rounding. Written to <name>_combined_table.tex, with full notes.
alpha_txt  <- sub("^0\\.", ".", format(ALPHA))
TAB_GROUPS <- list(                                             #the rows: column of the rates file = row label
  c(proj = "Projection test", proj_pval = sprintf("Projection test, $p \\le %s$ only", alpha_txt)),
  c(T_inf = "Kudrin's $T_\\infty$", T_H = "Kudrin's $T_H$", T_H_hybrid = "Kudrin's $T_H^{\\mathrm{hyb}}$"),
  c(EKW_LCM = "EKW: LCM", EKW_Fisher = "EKW: Fisher", EKW_Discontinuity = "EKW: density discontinuity",
    EKW_Binomial = "EKW: binomial", EKW_CS1 = "EKW: CS1", EKW_CS2B = "EKW: CS2B", EKW_min = "EKW: smallest $p$-value"),
  c(stouffer = "$p$-curve (Stouffer)", EKW_pcurve_min = "EKW or $p$-curve: smallest $p$-value"))
OPT_TEX    <- c(none = "None", uniform = "Unif.", smin = "$s_{\\min}$", smin_noise = "$+$noise")   #as the run's tables
METHOD_TEX <- c(max = "Most significant", threshold = "Threshold")   #short, so each group is as wide as its numbers
upper <- ALPHA + qnorm(0.975) * sqrt(ALPHA * (1 - ALPHA) / n_reps)      #the top of the 95% Monte Carlo band
rate_at <- function(m, o, tau, col) {
  if (!col %in% names(SM)) return(NA_real_)
  z <- SM[SM$method == m & SM$option == o & abs(SM$tau - tau) < 1e-9, col]
  if (length(z)) as.numeric(z[1]) else NA_real_
}
fmt_rate <- function(v, size = FALSE) {
  if (!is.finite(v)) return("--")
  s <- sub("^0\\.", ".", sprintf("%.3f", v))
  if (size && TAB_MARK_SIZE && v > upper) paste0("\\textbf{", s, "}") else s
}
no <- length(OPTIONS)
ng <- 1 + length(METHODS)
head_tex <- c(
  paste0(" & ", paste(sprintf("\\multicolumn{%d}{c}{%s}", no, c("Size ($\\tau = 0$)", sprintf("%s ($\\tau = 1$)", METHOD_TEX[METHODS]))),
                      collapse = " & "), " \\\\"),
  paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + (seq_len(ng) - 1) * no, 1 + seq_len(ng) * no), collapse = " "),
  paste0("Test & ", paste(rep(OPT_TEX[OPTIONS], ng), collapse = " & "), " \\\\"))
body_tex <- unlist(lapply(seq_along(TAB_GROUPS), function(g) {
  rows <- vapply(names(TAB_GROUPS[[g]]), function(col) {
    v <- c(vapply(OPTIONS, function(o) fmt_rate(rate_at(METHODS[1], o, 0, col), size = TRUE), ""),   #tau = 0: the same for both
           unlist(lapply(METHODS, function(m) vapply(OPTIONS, function(o) fmt_rate(rate_at(m, o, 1, col)), ""))))
    paste0(TAB_GROUPS[[g]][[col]], " & ", paste(v, collapse = " & "), " \\\\")
  }, "")
  c(if (g > 1) "\\addlinespace", rows)
}))
tab_notes <- paste0(
  sprintf("Rejection rates of %s\\%%-level tests over %d repetitions, with no switching ($\\tau = 0$: size) and when every trial switches ($\\tau = 1$). ",
          format(100 * ALPHA), n_reps),
  sprintf(paste0("Each meta-sample draws %s of the %s randomized trials in the Cochrane Database of Systematic Reviews (Schwab's data set, as in ",
                 "van Zwet, Schwab and Senn 2021) with at least %s continuous efficacy outcomes of at least %s participants, without replacement, ",
                 "and reports one outcome per trial. "),
          facts[["n_meta"]], facts[["n_elig"]], facts[["min_out"]], facts[["n_min"]]),
  sprintf(paste0("A trial chooses among its first %s outcomes in review order. With probability $1 - \\tau$ it reports one chosen at random ",
                 "(its pre-specified outcome); with probability $\\tau$ it reports its most significant outcome (Most significant) or, when the ",
                 "pre-specified outcome has $|z| \\le %s$, the larger of it and one other outcome (Threshold). "),
          facts[["max_out"]], format(CV)),
  "Trials select on unrounded $z$: each repetition draws the means and standard deviations uniformly within their rounding cells. ",
  "Columns within each group: the analyst's treatment of rounding. None: the printed $z$. Unif.: means and standard deviations drawn ",
  "uniformly within their rounding cells. $s_{\\min}$: the mean difference over $SE_{\\min}$, the standard error from the smallest standard ",
  "deviations consistent with the printed ones. $+$noise: the same with uniform noise at each mean's last digit. ",
  sprintf("Outcomes whose rounding unit in $z$-units, $\\kappa = 2a/SE_{\\min}$ ($a$ the unit of the means), exceeds %s are left out. ",
          facts[["kappa_cap"]]),
  sprintf(paste0("Projection test: $J = 40$, $J_0 = 20$, centered at 1.96, %d bootstrap repetitions; under $s_{\\min}$ and $+$noise it rejects ",
                 "only when its breakdown statistic $\\hat B$ exceeds the rounding allowance %s times the mean of $\\kappa^6$ over the reported ",
                 "outcomes (Propositions 3 and 4), and otherwise when $p \\le %s$; the row ``$p \\le %s$ only'' ignores the allowance. "),
          as.integer(PROJ_BOOTS), format(ALLOW_C), alpha_txt, alpha_txt),
  sprintf("Kudrin (2025): $T_\\infty$, with %d bootstrap repetitions, and $T_H$ and $T_H^{\\mathrm{hyb}}$, with the sinc kernel and frequency cutoff 3. ",
          as.integer(KUDRIN_BOOTS)),
  paste0("EKW: the tests of Elliott, Kudrin and W\\\"uthrich (2022), as in the application, on the two-sided $p$-values below 0.15: LCM and ",
         "Fisher, both of a concave $p$-curve CDF, the density discontinuity test at .05, the binomial test on $[.04,.05]$, and the Cox--Shi ",
         "histogram tests of monotonicity (CS1) and 2-monotonicity with bounds (CS2B). "),
  "$p$-curve: the test for left skew of Simonsohn, Nelson and Simmons (2014), by Stouffer's method. ",
  sprintf("The two ``smallest $p$-value'' rows reject when the smallest of their tests' $p$-values is at most %s; they are not size-controlled tests. ",
          alpha_txt),
  "Every test treats trials as independent units.",
  if (TAB_MARK_SIZE) sprintf(" Bold: a size above %s, the upper end of the 95\\%% Monte Carlo band around %s\\%%.",
                             sub("^0\\.", ".", sprintf("%.3f", upper)), format(100 * ALPHA)))
tab_tex <- c("\\begin{table}[htbp]", "\\centering",
             "\\caption{Size and power of every test against outcome switching in Cochrane trials}",
             paste0("\\label{", TAB_LABEL, "}"),
             "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
             "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
             sprintf("\\begin{tabular}{l%s}", strrep("c", no * ng)), "\\toprule", head_tex, "\\midrule", body_tex,
             "\\bottomrule", "\\end{tabular}}", "\\par\\smallskip",
             paste0("\\begin{minipage}{\\linewidth}\\footnotesize\\textit{Notes:} ", tab_notes, "\\end{minipage}"),
             "\\end{table}")
writeLines(tab_tex, paste0(OUT_STEM, "_table.tex"))
cat("\nWritten to ", RUN_DIR, ":\n  ", basename(OUT_STEM), ".pdf (the figure)\n  ", basename(OUT_STEM),
    "_caption.txt (its notes)\n  ", basename(OUT_STEM), "_rates.csv (the plotted rates)\n  ", basename(OUT_STEM),
    "_table.tex (every test at tau = 0 and tau = 1, for the appendix)\n", sep = "")


#### 6. The figure without Kudrin's tests --------------------------------------------------------------------
# The same figure without its two Kudrin rows (T_H and T_inf), so the projection test and the EKW-or-p-curve row:
# draw_combined() as it is, on section 3's layout less those rows' lines, so the panels, headers, axes and legend
# are those of <name>_combined.pdf and only the height changes. The script's own ROWS, nr, H and M are untouched,
# and the caption file serves both figures, as its notes never mention Kudrin's tests. Writes <name>_combined_nokudrin.pdf.
nk_r   <- which(!names(ROWS) %in% c("T_H", "T_inf"))              #the rows kept
nk_lay <- c(1, rbind(2 * nk_r, 2 * nk_r + 1), nrow(M) - 1:0)      #their lines of the layout (header, panels), with the
#column titles, the axis titles and the legend
nk_M   <- M[nk_lay, , drop = FALSE]
nk_M[] <- match(nk_M, sort(unique(c(nk_M))))                      #the cells renumbered 1, 2, ... in the order they are drawn
nk_draw <- draw_combined                                          #draw_combined(), finding the kept rows and layout first
environment(nk_draw) <- list2env(list(ROWS = ROWS[nk_r], nr = length(nk_r), H = H[nk_lay], M = nk_M), parent = globalenv())
grDevices::pdf(paste0(OUT_STEM, "_nokudrin.pdf"), width = sum(W) + 0.02, height = sum(H[nk_lay]) + 0.02, pointsize = 8)
nk_draw()
grDevices::dev.off()
cat("  ", basename(OUT_STEM), "_nokudrin.pdf (the figure without Kudrin's tests)\n", sep = "")