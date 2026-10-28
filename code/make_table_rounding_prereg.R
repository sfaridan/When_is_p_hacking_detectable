# Combine the per-option results of application_PAP_prereg.R into one table: every test, by pre-registration
# status, under each treatment of rounding. The projection test's breakdown statistic and the rounding allowance
# sit right under its p-values, and significant results are in bold.
# It reads, from results/:
#   <tag>_Table_application_numbers_<option>  the full-precision numbers make_application_table() saves next to
#                                             each .tex (run application_PAP_prereg.R once per option)
#   <tag>_Sample                              the sample those runs share, with each t-statistic's kappa
#   <tag>_Sample_rows_<option>                each finished run's stamp: the merged.dta rows it used
# and writes results/table_rounding_<tag>.tex. It stops unless every column was run on the identical sample.

root     <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
res_dir  <- file.path(root, "results")
tag      <- "PAP_cs_kappa0p3"           #file_tag of the runs: "PAP" for BCHH's exact sample, "PAP_ekwinput..." for their figure_a12
options  <- c("None" = "none", "Unif." = "uniform", "$s_{\\min}$" = "smin")   #column label = option key; add "$s_{\\min}$+noise" = "smin_noise" for option 4
J0       <- 20                          #the projection test's J0 to report (the allowance constant exists for 20 and 30)
J        <- 40                          #the settings of application_PAP_prereg.R, for the notes
numgrid  <- 3000
L        <- 6.5
boots    <- 1000
eps_disc <- 0.00063                     #discretization error E(L, delta), as printed in the per-option tables
p_style  <- "decimals"                  #"decimals": .002, .046, <.001; "signif": three significant digits (.00200, .0459), <.001 below .001
alpha    <- 0.05                        #level for bold
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

#--- The sample: one per table. Every column's results must carry the stamp of exactly this sample.
f_sample <- file.path(res_dir, paste0(tag, "_Sample"))
if (!file.exists(f_sample)) stop("No file ", f_sample, ": run the current application_PAP_prereg.R, which saves it")
smp <- readRDS(f_sample)
for (opt in unname(options)) {
  f <- file.path(res_dir, paste0(tag, "_Sample_rows_", opt))
  if (!file.exists(f)) stop("No file ", basename(f), ": rerun application_PAP_prereg.R with this option so its results are stamped with their sample")
  if (!identical(readRDS(f), smp$rows)) stop("The ", opt, " results were run on a different sample from ", basename(f_sample), ": rerun that option")
}
n_smp <- as.numeric(table(factor(smp$group, grp)))
for (i in seq_along(nums)) if (!isTRUE(all.equal(as.numeric(nums[[i]][grp, "n"]), n_smp)))
  stop("The ", options[[i]], " numbers file has n = ", paste(nums[[i]][grp, "n"], collapse = "/"), ", but the sample has ", paste(n_smp, collapse = "/"))

#--- Rounding allowance: C * mean(kappa^6) over that sample (Propositions 3 and 4), valid under s_min and s_min plus noise only
C_tab <- c(`20` = 0.13, `30` = 0.31)                           #the constant C of Proposition 3 for J0 = 20 and 30
if (!as.character(J0) %in% names(C_tab)) stop("No rounding-allowance constant for J0 = ", J0)
C_allow   <- C_tab[[as.character(J0)]]
allow     <- setNames(C_allow * as.numeric(tapply(smp$kappa^6, factor(smp$group, grp), mean)), grp)
has_allow <- unname(options) %in% c("smin", "smin_noise")       #treatments the allowance covers

#--- Cell formats
k <- length(options)
fmt_p <- function(x, bold = FALSE) {                           #p-values; no leading zero
  if (is.na(x)) return("$\\cdot$")
  if (x < 0.001) return(if (bold) "\\textbf{{\\boldmath${<}$}\\,.001}" else "${<}\\,.001$")
  s <- sub("^0\\.", ".", if (p_style == "signif") sprintf("%#.3g", x) else sprintf("%.3f", x))
  if (bold) paste0("\\textbf{", s, "}") else s
}
fmt_d <- function(x, bold = FALSE) {                           #distances, x 1000, three significant digits
  if (is.na(x)) return("$\\cdot$")
  v <- 1000 * x
  if (v == 0) return("0")
  if (v < 0.001) return("${<}\\,0.001$")
  if (v >= 100) return(sprintf("%.0f", signif(v, 3)))
  sprintf("%#.3g", v)
}
sig_p     <- function(g, i, x) !is.na(x) && x < alpha            #bold rule of every test but the projection test
proj_rule <- function(g, i, x) {                                 #projection test: B-hat above the allowance where one exists
  if (has_allow[i]) {
    B <- as.numeric(nums[[i]][g, Bcol])[1]
    !is.na(B) && !is.na(allow[[g]]) && B > allow[[g]]
  } else !is.na(x) && x < alpha
}
cells <- function(col, f, bold = function(g, i, x) FALSE) {     #columns: groups outer, treatments inner
  unlist(lapply(grp, function(g) vapply(seq_len(k), function(i) {
    x <- as.numeric(nums[[i]][g, col])[1]
    f(x, isTRUE(bold(g, i, x)))
  }, character(1))))
}
allow_cells <- unlist(lapply(grp, function(g) vapply(seq_len(k), function(i) if (has_allow[i]) fmt_d(allow[[g]]) else "$\\cdot$", character(1))))
row   <- function(label, x) paste0(paste(c(label, x), collapse = " & "), " \\\\")
block <- function(text) paste0("\\multicolumn{", 1 + 2 * k, "}{@{}l}{", text, "} \\\\")
fmtn   <- function(x) formatC(as.integer(round(x)), format = "d", big.mark = ",")
pctf   <- function(a, b) sprintf("%.1f\\%%", 100 * a / b)
mc     <- function(x) paste0("\\multicolumn{", k, "}{c}{", x, "}")   #one cell spanning a group's columns
base_n <- if (!is.null(smp$n_bchh)) smp$n_bchh else smp$n_before   #BCHH's tests, before the two rules
kept_n <- as.numeric(nums[[1]][grp, "n"])
row_tests <- mc(paste0(fmtn(kept_n), " of ", fmtn(base_n), " (", pctf(kept_n, base_n), ")"))
row_win   <- if (is.null(smp$w_fig5)) rep("$\\cdot$", 2) else
  mc(paste0(fmtn(smp$w_kept), " of ", fmtn(smp$w_fig5), " (", pctf(smp$w_kept, smp$w_fig5), ")"))
row_art   <- mc(as.character(nums[[1]][grp, "articles"]))

#--- Body: tests as rows; columns are the groups, each split by treatment of rounding
body <- c(
  block("Projection test"),
  row("\\quad $p$-value", cells(pcol, fmt_p, proj_rule)),
  row("\\quad $\\widehat{B}$ ($\\times 10^{3}$)", cells(Bcol, fmt_d)),
  row("\\quad Rounding allowance ($\\times 10^{3}$)", allow_cells),
  "\\addlinespace",
  block("Kudrin (2025)"),
  row("\\quad $T_\\infty$", cells("T_inf", fmt_p, sig_p)),
  row("\\quad $T_H$", cells("T_H", fmt_p, sig_p)),
  if (hyb_differs) row("\\quad $T_H^{\\mathrm{hybrid}}$", cells("T_H_hybrid", fmt_p, sig_p)),
  "\\addlinespace",
  block("Elliott et al.\\ (2022)"),
  row("\\quad LCM", cells("LCM", fmt_p, sig_p)),
  row("\\quad Fisher", cells("Fisher", fmt_p, sig_p)),
  row("\\quad Disc.", cells("Disc", fmt_p, sig_p)),
  row("\\quad Binomial", cells("Binomial", fmt_p, sig_p)),
  row("\\quad CS1", cells("CS1", fmt_p, sig_p)),
  row("\\quad CS2B", cells("CS2B", fmt_p, sig_p)),
  "\\midrule",
  row("Tests kept, of BCHH's", row_tests),
  row("Tests with $p<.15$, of BCHH's Figure 5", row_win),
  row("Articles", row_art))

#--- Notes
cs <- isTRUE(smp$settings$require_coef_se)                   #only tests with a reported coefficient and positive SE?
treat_text <- if (cs) {
  c(none       = "each $t$-statistic is the ratio of the reported coefficient and standard error.",
    uniform    = "the coefficient and the standard error are each replaced by an independent uniform draw from its rounding interval; %s.",
    smin       = "each coefficient is divided by the smallest standard error consistent with the reported one (Section~\\ref{sec:rounding}); this treatment is deterministic.",
    smin_noise = "as $s_{\\min}$, with uniform noise at its last digit added to each coefficient; %s.")
} else {
  c(none       = "each $t$-statistic is the $z$-statistic of \\citet{brodeur2024} as stored: the ratio of the reported coefficient and standard error, or a conversion of a reported $t$-statistic or $p$-value.",
    uniform    = "where the coefficient and standard error are reported, each is replaced by an independent uniform draw from its rounding interval, and other $t$-statistics receive uniform noise at their last digit; %s.",
    smin       = "where the coefficient and standard error are reported, the coefficient is divided by the smallest standard error consistent with the reported one, and other $t$-statistics are used as reported (Section~\\ref{sec:rounding}); this treatment is deterministic.",
    smin_noise = "as $s_{\\min}$, with uniform noise at its last digit added to each coefficient and to every other $t$-statistic; %s.")
}
avg_txt <- function(key) paste("$p$-values and $\\widehat{B}$ are averages over",
                               if (is.na(n_draws[[key]])) "independent draws" else paste(n_draws[[key]], "draws"))
treat_notes <- vapply(seq_len(k), function(i) {
  txt <- treat_text[[options[[i]]]]
  if (grepl("%s", txt, fixed = TRUE)) txt <- sprintf(txt, avg_txt(options[[i]]))
  paste0("\\emph{", names(options)[i], "}: ", txt)
}, character(1))
cap    <- smp$settings$kappa_cap
n_bchh <- if (!is.null(smp$n_bchh)) smp$n_bchh else smp$n_before   #sample files from before require_coef_se have n_before
n_cs   <- if (!is.null(smp$n_rule1)) smp$n_rule1 else if (!is.null(smp$n_cs)) smp$n_cs else n_bchh   #after rule 1
pct    <- function(a, b) sprintf("%.1f and %.1f percent", 100 * (1 - a[1] / b[1]), 100 * (1 - a[2] / b[2]))
drop_txt <- paste0(
  if (cs) paste0("we keep only tests reported with both a coefficient and a standard error, which drops ", pct(n_cs, n_bchh),
                 " of the tests without and with pre-registration") else "",
  if (cs && !is.null(cap)) "; we then " else if (!is.null(cap)) "we " else "",
  if (!is.null(cap)) paste0("drop every $t$-statistic whose rounding unit in $t$-units, $\\kappa$, exceeds ", format(cap),
                            if (cs) ", counting a standard error printed as zero as infinite $\\kappa$" else "", " (Section~\\ref{sec:rounding}): ",
                            if (cs) paste("a further", pct(smp$n_after, n_cs)) else paste(pct(smp$n_after, n_cs), "of the tests without and with pre-registration")) else "",
  if (cs && !is.null(cap)) ". These rules depend only on how results are printed, so they drop the same tests under every treatment."
  else if (cs || !is.null(cap)) ". The rule depends only on how results are printed, so it drops the same tests under every treatment."
  else "no test is dropped for rounding.")
sample_note <- paste("The sample is the same in every column.", if (identical(smp$settings$bchh_sample, "ekw_input")) {
  "It starts from the sample the replication files of \\citet{brodeur2024} pass to the tests of \\citet{elliott2022}, every test statistic with a recorded pre-registration status and a nonnegative $z$-statistic;"
} else {
  "It starts from that of \\citet{brodeur2024}, as in their replication files, without articles in \\textit{Experimental Economics} or studies that are not experiments;"
}, drop_txt)
allow_cols <- paste(names(options)[has_allow], collapse = " and ")
ex <- smp$extra_fig5
fig5_extra <- if (!is.null(ex) && ex[["tests"]] > 0) sprintf(paste0(
  ", which come from the input file of their EKW tests and so also include %s tests (%s with $p<.15$) from \\textit{Experimental Economics} ",
  "and from studies that are not experiments, which their other analyses drop"), fmtn(ex[["tests"]]), fmtn(ex[["window"]])) else ""
notes <- c(
  "The table reports tests of the null hypothesis of no selective reporting in the data of \\citet{brodeur2024}, separately for studies that were and were not pre-registered, meaning registered before the trial end date.",
  "Columns differ in the treatment of rounding.",
  treat_notes,
  sample_note,
  sprintf("\\emph{Projection test}: $J_0=%d$, $J=%d$, and a grid of %s points on $[-%s,%s]$ plus the point 9999; the discretization error $E(L,\\delta)=%.5f$ is included in every $p$-value and in $\\widehat{B}$.",
          J0, J, formatC(numgrid, format = "d", big.mark = ","), L, L, eps_disc),
  "$\\widehat{B}$ is its breakdown statistic (Section~\\ref{sec:breakdown}), a lower 95\\% confidence bound on the distance between the reported $t$-curve and the set of $t$-curves possible without selective reporting.",
  sprintf(paste0("The rounding allowance, %.2f times the average of $\\kappa^{6}$ over the $t$-statistics used (Propositions~\\ref{prop:rounding} and~\\ref{prop:smin}), bounds the distance that rounding can create under %s, ",
                 "where $\\kappa$ is a $t$-statistic's rounding unit in $t$-units; ",
                 if (cs) "it covers every $t$-statistic in the sample" else "it does not cover $t$-statistics converted from rounded $p$-values",
                 ", and no such bound is available under the other treatments."),
          C_allow, allow_cols),
  paste0("\\emph{Kudrin (2025)}: the KS-type ($T_\\infty$) and histogram-based ($T_H$) tests of \\citet{kudrin2025}, using the sinc kernel with frequency cutoff 3 and, for $T_H$, the points $0.01, 0.02, \\ldots, 0.15$",
         if (hyb_differs) "; $T_H^{\\mathrm{hybrid}}$ is his hybrid test." else "; his hybrid test coincides with $T_H$ in every case."),
  paste0("\\emph{Elliott et al.\\ (2022)}: the tests of \\citet{elliott2022}, namely the least concave majorant test (LCM) and Fisher's test on $p$-values in $(0,0.15]$, ",
         "the density discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$, and histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) ",
         "and for 2-monotonicity with bounds (CS2B); as in the input files of \\citet{brodeur2024}, they use the $p$-values $2(1-\\Phi(|t|))$ in single precision."),
  sprintf("The projection test, $T_\\infty$, $T_H$, CS1 and CS2B cluster by article, and the projection test and $T_\\infty$ use %s article-level bootstrap repetitions.",
          formatC(boots, format = "d", big.mark = ",")),
  if (p_style == "signif") "$p$-values are reported to three significant digits, and those below .001 as ${<}\\,.001$."
  else "$p$-values are rounded to three decimals, and those below .001 are shown as ${<}\\,.001$.",
  sprintf(paste0("Bold: $p<%s$, except that the projection test under %s is bold only when $\\widehat{B}$ exceeds the rounding allowance, ",
                 "so that its rejection cannot be attributed to rounding."), sub("^0\\.", ".", format(alpha)), allow_cols),
  paste0("Tests kept, of BCHH's: the tests in our sample, out of the ", fmtn(sum(base_n)), " in theirs. ",
         "Tests with $p<.15$, of BCHH's Figure 5: the tests in the window of the EKW tests, counted with BCHH's $z$-statistics in both samples ",
         "and set against the counts in their Figure 5", fig5_extra, ". Articles: the number of articles in each group."))

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
tests <- c(Projection = pcol, B_x1000 = Bcol, T_inf = "T_inf", T_H = "T_H", if (hyb_differs) c(T_H_hybrid = "T_H_hybrid"),
           LCM = "LCM", Fisher = "Fisher", Disc = "Disc", Binomial = "Binomial", CS1 = "CS1", CS2B = "CS2B")
prev <- do.call(cbind, lapply(grp, function(g) sapply(nums, function(d) as.numeric(unlist(d[g, tests])))))
dimnames(prev) <- list(names(tests), paste(rep(grp, each = k), rep(unname(options), 2), sep = ":"))
prev["B_x1000", ] <- 1000 * prev["B_x1000", ]
prev <- rbind(prev, Allow_x1000 = unlist(lapply(grp, function(g) ifelse(has_allow, 1000 * allow[[g]], NA))))
print(signif(prev, 3))
cat("Draws averaged over:", paste(names(n_draws), n_draws, sep = " = ", collapse = ", "), "\n")
cat("LaTeX table written to", out_file, "\n")