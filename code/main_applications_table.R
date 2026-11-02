#### =====================================================================================================
#### main_applications_table.R
#### The main-text table of the two applications: Brodeur, Cook and Heyes (2020), by research design, and the
#### Cochrane Database of Systematic Reviews. For each sample, every test on the printed statistics and on
#### de-rounded ones, and the projection test's certified result under s_min. It reads only what
#### application_MM_analysis.R and application_cochrane.R saved: run them first, clear R, then run this.
#### Writes table_main_applications.tex in out_dir (needs \usepackage{booktabs, threeparttable, multirow}).
#### The full results (every rounding treatment, standard deviations across draws) are in the appendix tables.
#### =====================================================================================================
rm(list = ls())

#### 0. Settings -------------------------------------------------------------------------------------------
root      <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
out_dir   <- file.path(root, "results")
mm_dir    <- file.path(root, "results")                       #application_MM_analysis.R's res_dir
mm_tag    <- "_keepmu0"                                      #its run_tag
co_dir    <- file.path(root, "results", "cochrane")           #application_cochrane.R's res_dir
co_tag    <- "_n100_all_trial_kappa0.5_zcounts_cont"          #its tag
co_groups <- "all"                                           #Cochrane subsamples shown: any of "all", "early", "late"
deround   <- "uniform"                                       #the de-rounded rows: "uniform" or "smin_noise"
certify   <- "smin"                                          #the certified column: "smin" or "smin_noise"
J0        <- 20                                              #the projection test's main specification
ref_fig_tinf <- "Figure~\\ref{fig:cochrane_tinf_gap}"         #your labels for the appendix figure and tables
ref_appendix <- "Appendix Tables~\\ref{tab:brodeur_rounding} and~\\ref{tab:cochrane}"

stopifnot(deround %in% c("uniform", "smin_noise"), certify %in% c("smin", "smin_noise"))
meth <- c("RCT", "IV", "DID", "RDD")

#### 1. Brodeur, Cook and Heyes (2020): application_MM_analysis.R's results ---------------------------------
ms_file <- file.path(mm_dir, paste0("Settings_application", mm_tag))
if (!file.exists(ms_file)) stop("No ", ms_file, ": run application_MM_analysis.R first, or set mm_tag")
MS <- readRDS(ms_file)
if (!identical(MS$status, "complete")) stop("The Brodeur run with run tag '", mm_tag, "' did not finish")
for (k in c("none", deround, certify)) {
  if (!k %in% MS$keys_run) stop("The Brodeur run did not include the de-rounding option '", k, "'")
}
if (!J0 %in% MS$J0_list) stop("The Brodeur run did not include J0 = ", J0)

mm_get <- function(key) {                                     #one row per design: every test under option key
  f <- c(proj = paste0("MM_weights_results_", J0, "_J0_", key, mm_tag), kud = paste0("Kudrin_results_", key, mm_tag),
         ekw = paste0("EWK_results_", key, mm_tag), st = paste0("Stouffer_results_", key, mm_tag))
  f <- setNames(file.path(mm_dir, f), names(f))
  miss <- f[!file.exists(f)]
  if (length(miss)) stop("Missing Brodeur results: ", paste(basename(miss), collapse = ", "))
  pr <- readRDS(f[["proj"]]); kd <- readRDS(f[["kud"]]); ek <- readRDS(f[["ekw"]]); st <- readRDS(f[["st"]])
  num <- function(v) as.numeric(unlist(v))
  ekw_min <- num(ek["min", paste0("ewk_", tolower(meth))])     #the smallest of the six averaged EKW p-values
  data.frame(sample = meth,
             proj = num(pr$pval)[match(meth, pr$methods)], B = num(pr$breakdown)[match(meth, pr$methods)],
             T_inf = num(kd$T_inf)[match(meth, kd$methods)], T_H = num(kd$T_H)[match(meth, kd$methods)],
             T_H_hybrid = num(kd$T_H_hybrid)[match(meth, kd$methods)],
             ekw_st = pmin(ekw_min, num(st$pval)[match(meth, st$methods)], na.rm = TRUE),
             stringsAsFactors = FALSE)
}
mm_printed <- mm_get("none")
mm_derounded <- mm_get(deround)
mm_cert <- mm_get(certify)
mm_allow <- if (!is.null(MS$allowance)) MS$allowance[[as.character(J0)]][meth] else MS$C_round[[as.character(J0)]] * MS$kappa6[meth]
mm_size <- sprintf("%s tests, %s articles", format(MS$n_by_design[meth], big.mark = ","),
                   format(MS$articles_by_design[meth], big.mark = ","))

#### 2. Cochrane: application_cochrane.R's results ----------------------------------------------------------
co_file <- file.path(co_dir, paste0("Cochrane_results", co_tag, ".rds"))
if (!file.exists(co_file)) stop("No ", co_file, ": run application_cochrane.R first, or set co_tag")
RES <- readRDS(co_file)
cs  <- RES$settings
for (k in c("none", deround, certify)) {
  if (!k %in% cs$treatments) stop("The Cochrane run did not include the treatment '", k, "'")
}
if (!J0 %in% cs$J0_list) stop("The Cochrane run did not include J0 = ", J0)
co_groups <- intersect(co_groups, unique(RES$results$group))
if (!length(co_groups)) stop("None of co_groups is in the Cochrane results")

co_get <- function(tr) {                                     #one row per subsample: every test under treatment tr
  r  <- RES$results[RES$results$treatment == tr, , drop = FALSE]
  r  <- r[match(co_groups, r$group), , drop = FALSE]
  ek <- intersect(c("LCM", "Fisher", "Disc", "Binomial", "CS1", "CS2B", "Stouffer"), names(r))
  data.frame(sample = co_groups,
             proj = r[[paste0("proj_p_J0_", J0)]], B = r[[paste0("proj_B_J0_", J0)]],
             T_inf = r$T_inf, T_H = r$T_H, T_H_hybrid = r$T_H_hybrid,
             ekw_st = apply(r[, ek, drop = FALSE], 1, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}
co_printed   <- co_get("none")
co_derounded <- co_get(deround)
co_cert      <- co_get(certify)
co_rows      <- RES$results[RES$results$treatment == certify, , drop = FALSE]
co_rows      <- co_rows[match(co_groups, co_rows$group), , drop = FALSE]
co_allow     <- co_rows[[paste0("allowance_J0_", J0)]]
co_label     <- c(all = "All trials", early = sprintf("Published by %d", cs$year_split),
                  late = sprintf("Published after %d", cs$year_split))[co_groups]
co_size      <- sprintf("%s outcomes, %s trials", format(co_rows$n, big.mark = ","), format(co_rows$clusters, big.mark = ","))

#### 3. Do the two runs match? -----------------------------------------------------------------------------
same <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b)))
if (!same(MS$kappa_cap, cs$kappa_cap)) warning("The kappa caps differ: Brodeur ", MS$kappa_cap, ", Cochrane ", cs$kappa_cap, call. = FALSE)
if (!same(MS$n_draws_noise, cs$n_draws)) warning("The numbers of de-rounding draws differ: Brodeur ", MS$n_draws_noise,
                                                 ", Cochrane ", cs$n_draws, call. = FALSE)
if (!same(MS$boots, cs$boots)) warning("The bootstrap repetitions differ: Brodeur ", MS$boots, ", Cochrane ", cs$boots, call. = FALSE)
if (!isFALSE(MS$boot_seed_by_draw)) warning("The Brodeur run predates the common bootstrap seeds (boot_seed_by_draw = FALSE): ",
                                            "rerun application_MM_analysis.R so both applications use the same convention", call. = FALSE)
n_draws <- as.integer(cs$n_draws)

#### 4. The table ------------------------------------------------------------------------------------------
fp <- function(x) {                                           #a p-value: .974, <.001; bold when at most .05
  if (!is.finite(x)) return("--")
  s <- if (x < 0.0005) "$<$.001" else sub("^0\\.", ".", sprintf("%.3f", x))
  if (x <= 0.05) paste0("\\textbf{", s, "}") else s
}
fB <- function(B, a) {                                        #10^3 B-hat (10^3 allowance); bold when B-hat exceeds it
  if (!is.finite(B) || !is.finite(a)) return("--")
  fa <- if (1e3 * a < 0.05) "$<$0.1" else sprintf("%.1f", 1e3 * a)
  s  <- sprintf("%.1f (%s)", 1e3 * B, fa)
  if (B > a) paste0("\\textbf{", s, "}") else s
}
pair <- function(label, size, pr, de, cert) {                 #the two rows of one sample
  c(sprintf("%s & Printed & %s & %s & %s & %s & %s & \\multirow{2}{*}{%s} \\\\", label,
            fp(pr$proj), fp(pr$T_inf), fp(pr$T_H), fp(pr$T_H_hybrid), fp(pr$ekw_st), cert),
    sprintf("{\\footnotesize %s} & De-rounded & %s & %s & %s & %s & %s & \\\\", size,
            fp(de$proj), fp(de$T_inf), fp(de$T_H), fp(de$T_H_hybrid), fp(de$ekw_st)))
}
body_mm <- unlist(lapply(seq_along(meth), function(i) {
  c(if (i > 1) "\\addlinespace[2pt]",
    pair(meth[i], mm_size[i], mm_printed[i, ], mm_derounded[i, ], fB(mm_cert$B[i], mm_allow[[i]])))
}))
body_co <- unlist(lapply(seq_along(co_groups), function(i) {
  c(if (i > 1) "\\addlinespace[2pt]",
    pair(co_label[[i]], co_size[i], co_printed[i, ], co_derounded[i, ], fB(co_cert$B[i], co_allow[i])))
}))

deround_text <- if (deround == "uniform") {
  "De-rounded: every printed number drawn uniformly within its rounding interval"
} else {
  "De-rounded: the smallest standard error consistent with the printed one, with uniform noise at the last printed digit of the estimate"
}
certify_text <- if (certify == "smin") "the smallest standard error consistent with the printed one" else
  "the smallest standard error consistent with the printed one, with uniform noise at the last printed digit of the estimate"
notes <- paste0(
  "Entries are $p$-values; bold marks $p \\le .05$. Printed: the test statistics as computed from the printed coefficients and standard errors (Panel A) or the printed arm means and standard deviations (Panel B). ",
  deround_text, sprintf("; $p$-values averaged over %d draws. ", n_draws),
  sprintf("Projection test: $J = 40$, $J_0 = %d$, centered at 1.96. Certified: the projection test under %s (Section 5.2); its rejection is certified when $10^3 \\hat B$ exceeds the rounding allowance $10^3 \\times %s\\,\\overline{\\kappa^6}$ (Propositions 3 and 4), shown in parentheses; bold when it does. ",
          J0, certify_text, format(if (J0 == 20) 0.13 else 0.31)),
  "EKW and $p$-curve: the smallest of the $p$-values of the six tests of Elliott, Kudrin and W\\\"uthrich (2022) and the $p$-curve test for left skew (Simonsohn, Nelson and Simmons 2014), each averaged over draws; it is not itself a valid $p$-value and overstates their evidence. ",
  "Every test clusters by article (Panel A) or trial (Panel B). Statistics whose rounding interval exceeds half a standard error are left out ($\\kappa \\le 0.5$); estimates printed as zero are kept. ",
  sprintf("Panel B: continuous efficacy outcomes of randomized trials with at least %d participants. ", as.integer(cs$n_min)),
  sprintf("$T_\\infty$'s rejections in Panel B come from $p$-values near 1: outcomes whose two arms have equal printed means (%s). ", ref_fig_tinf),
  sprintf("Every rounding treatment, and standard deviations across draws: %s.", ref_appendix))

tex <- c(sprintf("%% Written by main_applications_table.R on %s from Settings_application%s and Cochrane_results%s.rds",
                 format(Sys.time(), "%Y-%m-%d %H:%M"), mm_tag, co_tag),
         "% Needs \\usepackage{booktabs, threeparttable, multirow}",
         "\\begin{table}[!htbp]\\centering",
         "\\begin{threeparttable}",
         "\\caption{Tests for selective reporting in economics and medicine, on printed and de-rounded statistics}",
         "\\label{tab:main_applications}",
         "\\begin{tabular}{llcccccc}",
         "\\toprule",
         " & & Projection & \\multicolumn{3}{c}{Kudrin (2025)} & EKW and & Projection test \\\\",
         "\\cmidrule(lr){4-6}",
         "Sample & Statistics & test & $T_\\infty$ & $T_H$ & $T_H^{\\mathrm{hyb}}$ & $p$-curve & certified: $10^3 \\hat B$ \\\\",
         " & & $p$ & $p$ & $p$ & $p$ & smallest $p$ & (allowance) \\\\",
         "\\midrule",
         "\\multicolumn{8}{l}{\\textit{A. Economics: Brodeur, Cook and Heyes (2020), by research design}} \\\\",
         body_mm,
         "\\midrule",
         "\\multicolumn{8}{l}{\\textit{B. Medicine: Cochrane Database of Systematic Reviews, continuous efficacy outcomes}} \\\\",
         body_co,
         "\\bottomrule",
         "\\end{tabular}",
         "\\begin{tablenotes}[flushleft]\\footnotesize",
         paste0("\\item ", notes),
         "\\end{tablenotes}",
         "\\end{threeparttable}",
         "\\end{table}")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, "table_main_applications.tex")
writeLines(tex, out_file)

#### 5. The same numbers on the console, for a quick look -------------------------------------------------
show <- rbind(
  data.frame(panel = "Brodeur", sample = meth, statistics = "printed", mm_printed[, -1], allowance = NA),
  data.frame(panel = "Brodeur", sample = meth, statistics = deround, mm_derounded[, -1], allowance = NA),
  data.frame(panel = "Brodeur", sample = meth, statistics = certify, mm_cert[, -1], allowance = as.numeric(mm_allow)),
  data.frame(panel = "Cochrane", sample = co_groups, statistics = "printed", co_printed[, -1], allowance = NA),
  data.frame(panel = "Cochrane", sample = co_groups, statistics = deround, co_derounded[, -1], allowance = NA),
  data.frame(panel = "Cochrane", sample = co_groups, statistics = certify, co_cert[, -1], allowance = co_allow))
show <- show[order(show$panel, match(show$sample, c(meth, co_groups))), ]
show$B <- 1e3 * show$B; show$allowance <- 1e3 * show$allowance
names(show)[names(show) %in% c("B", "allowance")] <- c("B_x1000", "allowance_x1000")
print(format(show, digits = 3), row.names = FALSE)
cat("\nTable written to ", out_file, "\n", sep = "")