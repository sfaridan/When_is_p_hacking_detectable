# Combine the per-option results of application_PAP_prereg.R into one table: every test, by pre-registration
# status, under each treatment of rounding, with p-values to three significant digits.
# It reads results/<tag>_Table_application_numbers_<option>, the full-precision numbers make_application_table()
# saves next to each .tex (the .tex files themselves round to two decimals), so run application_PAP_prereg.R once
# per option first. It writes results/table_rounding_<tag>.tex.

root     <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
res_dir  <- file.path(root, "results")
tag      <- "PAP"                       #"PAP_ekwinput" for the sample of BCHH's figure_a12
options  <- c("None" = "none", "Unif." = "uniform", "$s_{\\min}$" = "smin")   #column label = option key; add "$s_{\\min}$+noise" = "smin_noise" for option 4
J0       <- 20                          #the projection test's J0 to report
J        <- 40                          #the settings of application_PAP_prereg.R, for the notes
numgrid  <- 3000
L        <- 6.5
boots    <- 1000
eps_disc <- 0.00063                     #discretization error E(L, delta), as printed in the per-option tables
B_digits <- 1                           #decimals of B-hat x 1000
out_file <- file.path(res_dir, paste0("table_rounding_", tag, ".tex"))

grp <- c("NoPreReg", "PreReg")
lab <- c(NoPreReg = "Not pre-registered", PreReg = "Pre-registered")

#--- Read the numbers, one data frame per option, rows in the order of grp
nums <- lapply(options, function(opt) {
  f <- file.path(res_dir, paste0(tag, "_Table_application_numbers_", opt))
  if (!file.exists(f)) stop("No file ", f, ": run application_PAP_prereg.R with this de-rounding option first")
  d <- readRDS(f)
  if (!all(grp %in% rownames(d))) stop(basename(f), " has no rows for ", paste(setdiff(grp, rownames(d)), collapse = ", "))
  d[grp, ]
})
n_draws <- setNames(vapply(unname(options), function(opt) {    #de-roundings each option averages over
  f <- file.path(res_dir, paste0(tag, "_weights_results_draws_", opt))
  if (file.exists(f)) length(unique(readRDS(f)$draw)) else NA_integer_
}, integer(1)), unname(options))

same <- function(col) all(vapply(nums, function(d) isTRUE(all.equal(d[[col]], nums[[1]][[col]])), logical(1)))
if (!same("n") || !same("articles")) warning("n or article counts differ across options: were they all run on the same sample?")
pcol <- paste0("p_J0_", J0)
Bcol <- paste0("B_J0_", J0)
if (!pcol %in% names(nums[[1]])) stop("The numbers files hold no results for J0 = ", J0)
hyb_differs <- any(vapply(nums, function(d) any(d$hybrid_restricted %in% TRUE), logical(1)))   #show T_H^hybrid only if it ever differs from T_H

#--- Cell formats
fmt_p <- function(x) {                                       #three significant digits, no leading zero
  if (is.na(x)) return("$\\cdot$")
  if (x <= 0)   return("${<}\\,.001$")                       #a bootstrap p-value of zero: below 1/1,000
  if (x < 0.001) {                                           #e.g. 4.93e-06 -> $4.93{\times}10^{-6}$
    s <- sprintf("%.2e", x)
    return(sprintf("$%s{\\times}10^{%d}$", sub("e.*$", "", s), as.integer(sub("^.*e", "", s))))
  }
  sub("^0\\.", ".", sprintf("%#.3g", x))                     #e.g. .0500, .150, 1.00
}
fmt_B <- function(x) if (is.na(x)) "$\\cdot$" else sprintf(paste0("%.", B_digits, "f"), 1000 * x)

#--- Body: tests as rows; columns are the groups, each split by treatment of rounding
k     <- length(options)
cells <- function(col, f) unlist(lapply(grp, function(g) vapply(nums, function(d) f(as.numeric(d[g, col])[1]), character(1))))
row   <- function(label, x) paste0(paste(c(label, x), collapse = " & "), " \\\\")
block <- function(text) paste0("\\multicolumn{", 1 + 2 * k, "}{@{}l}{", text, "} \\\\")
sizes <- vapply(grp, function(g) paste0("\\multicolumn{", k, "}{c}{", formatC(as.integer(round(nums[[1]][g, "n"])), format = "d", big.mark = ","),
                                        " (", nums[[1]][g, "articles"], ")}"), character(1))
body <- c(
  row("Projection test", cells(pcol, fmt_p)),
  "\\addlinespace",
  block("Kudrin (2025)"),
  row("\\quad $T_\\infty$", cells("T_inf", fmt_p)),
  row("\\quad $T_H$", cells("T_H", fmt_p)),
  if (hyb_differs) row("\\quad $T_H^{\\mathrm{hybrid}}$", cells("T_H_hybrid", fmt_p)),
  "\\addlinespace",
  block("Elliott et al.\\ (2022)"),
  row("\\quad LCM", cells("LCM", fmt_p)),
  row("\\quad Fisher", cells("Fisher", fmt_p)),
  row("\\quad Disc.", cells("Disc", fmt_p)),
  row("\\quad Binomial", cells("Binomial", fmt_p)),
  row("\\quad CS1", cells("CS1", fmt_p)),
  row("\\quad CS2B", cells("CS2B", fmt_p)),
  "\\midrule",
  row("$\\widehat{B}$ ($\\times 10^{3}$)", cells(Bcol, fmt_B)),
  row("Tests (articles)", sizes))

#--- Notes
treat_text <- c(
  none       = "each $t$-statistic is the $z$-statistic of \\citet{brodeur2024} as stored: the ratio of the reported coefficient and standard error, or a conversion of a reported $t$-statistic or $p$-value.",
  uniform    = "where the coefficient and standard error are reported, each is replaced by an independent uniform draw from its rounding interval, and other $t$-statistics receive uniform noise at their last digit; %s.",
  smin       = "where the coefficient and standard error are reported, the coefficient is divided by the smallest standard error consistent with the reported one, and other $t$-statistics are used as reported (Section~\\ref{sec:rounding}); this treatment is deterministic.",
  smin_noise = "as $s_{\\min}$, with uniform noise at its last digit added to each coefficient and to every other $t$-statistic; %s.")
avg_txt <- function(key) paste("$p$-values and $\\widehat{B}$ are averages over",
                               if (is.na(n_draws[[key]])) "independent draws" else paste(n_draws[[key]], "draws"))
treat_notes <- vapply(seq_len(k), function(i) {
  txt <- treat_text[[options[[i]]]]
  if (grepl("%s", txt, fixed = TRUE)) txt <- sprintf(txt, avg_txt(options[[i]]))
  paste0("\\emph{", names(options)[i], "}: ", txt)
}, character(1))
sample_note <- if (grepl("ekwinput", tag)) {
  "The sample, the same in every column, is the one the replication files of \\citet{brodeur2024} pass to the tests of \\citet{elliott2022}: every test statistic with a recorded pre-registration status and a nonnegative $z$-statistic."
} else {
  "The sample is that of \\citet{brodeur2024} in every column: as in their replication files, articles in \\textit{Experimental Economics} and studies that are not experiments are dropped, and no test is dropped for rounding."
}
notes <- c(
  "The table reports $p$-values for the null hypothesis of no selective reporting in the data of \\citet{brodeur2024}, separately for studies that were and were not pre-registered, meaning registered before the trial end date.",
  "Columns differ in the treatment of rounding.",
  treat_notes,
  sample_note,
  sprintf("\\emph{Projection test}: $J_0=%d$, $J=%d$, and a grid of %s points on $[-%s,%s]$ plus the point 9999; the discretization error $E(L,\\delta)=%.5f$ is included in every $p$-value and in $\\widehat{B}$.",
          J0, J, formatC(numgrid, format = "d", big.mark = ","), L, L, eps_disc),
  paste0("\\emph{Kudrin (2025)}: the KS-type ($T_\\infty$) and histogram-based ($T_H$) tests of \\citet{kudrin2025}, using the sinc kernel with frequency cutoff 3 and, for $T_H$, the points $0.01, 0.02, \\ldots, 0.15$",
         if (hyb_differs) "; $T_H^{\\mathrm{hybrid}}$ is his hybrid test." else "; his hybrid test coincides with $T_H$ in every case."),
  paste0("\\emph{Elliott et al.\\ (2022)}: the tests of \\citet{elliott2022}, namely the least concave majorant test (LCM) and Fisher's test on $p$-values in $(0,0.15]$, ",
         "the density discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$, and histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) ",
         "and for 2-monotonicity with bounds (CS2B); as in the input files of \\citet{brodeur2024}, they use the $p$-values $2(1-\\Phi(|t|))$ in single precision."),
  sprintf("The projection test, $T_\\infty$, $T_H$, CS1 and CS2B cluster by article, and the projection test and $T_\\infty$ use %s article-level bootstrap repetitions.",
          formatC(boots, format = "d", big.mark = ",")),
  "$p$-values are reported to three significant digits, and a $p$-value of exactly zero as ${<}\\,.001$.",
  "$\\widehat{B}$: the breakdown statistic of the projection test (Section~\\ref{sec:breakdown}), a lower 95\\% confidence bound on the distance between the reported $t$-curve and the set of $t$-curves possible without selective reporting.",
  "Tests (articles): the number of test statistics and, in parentheses, of articles in each group.")

#--- LaTeX (needs booktabs and threeparttable)
tex <- c(
  "\\begin{table}[t]",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{threeparttable}",
  paste0("\\caption{Tests for selective reporting by pre-registration status under ", c("one", "two", "three", "four")[k], " treatments of rounding}"),
  paste0("\\label{tab:rounding_", tag, "}"),
  paste0("\\begin{tabular}{@{}l*{", k, "}{c}@{\\hspace{12pt}}*{", k, "}{c}@{}}"),
  "\\toprule",
  row("", paste0("\\multicolumn{", k, "}{c}{", lab[grp], "}")),
  paste0("\\cmidrule(lr){2-", 1 + k, "} \\cmidrule(l){", 2 + k, "-", 1 + 2 * k, "}"),
  row("", rep(names(options), 2)),
  "\\midrule",
  body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  paste("\\item \\textit{Notes:}", paste(notes, collapse = " ")),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}")
writeLines(tex, out_file)

#--- Console preview at full precision
tests <- c(Projection = pcol, T_inf = "T_inf", T_H = "T_H", if (hyb_differs) c(T_H_hybrid = "T_H_hybrid"),
           LCM = "LCM", Fisher = "Fisher", Disc = "Disc", Binomial = "Binomial", CS1 = "CS1", CS2B = "CS2B", B_x1000 = Bcol)
prev <- do.call(cbind, lapply(grp, function(g) sapply(nums, function(d) as.numeric(unlist(d[g, tests])))))
dimnames(prev) <- list(names(tests), paste(rep(grp, each = k), rep(unname(options), 2), sep = ":"))
prev["B_x1000", ] <- 1000 * prev["B_x1000", ]
print(signif(prev, 3))
cat("Draws averaged over:", paste(names(n_draws), n_draws, sep = " = ", collapse = ", "), "\n")
cat("LaTeX table written to", out_file, "\n")