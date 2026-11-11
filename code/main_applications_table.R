#### =====================================================================================================
#### main_applications_table.R
#### The main-text table of the two applications: Brodeur, Cook and Heyes (2020), by research design, and the
#### Cochrane Database of Systematic Reviews. For each sample, every test on the printed statistics and on
#### de-rounded ones, and the projection test's certified result under s_min. It reads only what
#### application_MM_analysis.R and application_cochrane.R saved: run them first, clear R, then run this.
#### Writes table_main_applications.tex in out_dir (needs \usepackage{booktabs, threeparttable, multirow}).
#### The full results (every rounding treatment, standard deviations across draws) are in the appendix tables.
#### Section 6 also writes application_MM_output.R's two rounding tables with the Cochrane results below RDD:
#### table_rounding_ekw<mm_tag><rt_suffix>.tex and table_rounding_appendix<mm_tag><rt_suffix>.tex, from the same
#### saved results (no analysis is rerun).
#### =====================================================================================================
rm(list = ls())

#### 0. Settings -------------------------------------------------------------------------------------------
root      <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
out_dir   <- file.path(root, "results")
mm_dir    <- file.path(root, "results")                       #application_MM_analysis.R's res_dir
mm_tag    <- "_keepmu0_drop_se_mi"                                      #its run_tag
co_dir    <- file.path(root, "results", "cochrane")           #application_cochrane.R's res_dir
co_tag    <- "_n100_all_trial_kappa0.5_zcounts_cont"          #its tag
co_groups <- c("all", "early", "late")                       #Cochrane subsamples in the main table: any of "all", "early", "late"
deround   <- "uniform"                                       #the de-rounded rows: "uniform" or "smin_noise"
certify   <- "smin"                                          #the certified column: "smin" or "smin_noise"
J0        <- 20                                              #the projection test's main specification
ref_fig_tinf <- "Figure~\\ref{fig:cochrane_tinf_gap}"         #your labels for the appendix figure and tables
ref_appendix <- NULL   #NULL: the rounding table of section 6 (application_MM_output.R's if rt_write is FALSE) and tab:cochrane

#--- Section 6: application_MM_output.R's two rounding tables, with the Cochrane results below RDD
rt_write  <- TRUE                                            #FALSE: skip them
rt_groups <- "all"                                           #Cochrane subsamples there: any of "all", "early", "late"
rt_suffix <- "_cochrane"                                     #added to the file names and labels of application_MM_output.R's
#tables, so both versions can coexist
rt_bold   <- TRUE                                            #as there: bold p < .05; under s_min only certified projection tests

stopifnot(deround %in% c("uniform", "smin_noise"), certify %in% c("smin", "smin_noise"))
meth <- c("RCT", "IV", "DID", "RDD")
rt_file  <- c(main     = paste0("table_rounding_ekw", mm_tag, rt_suffix, ".tex"),
              appendix = paste0("table_rounding_appendix", mm_tag, rt_suffix, ".tex"))
rt_label <- c(main     = paste0("tab:rounding_ekw", mm_tag, rt_suffix),
              appendix = paste0("tab:rounding_appendix", mm_tag, rt_suffix))
if (is.null(ref_appendix)) ref_appendix <- sprintf("Appendix Tables~\\ref{%s} and~\\ref{tab:cochrane}",
                                                   if (rt_write) rt_label[["main"]] else paste0("tab:rounding_ekw", mm_tag))

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
co_by        <- if (is.null(cs$cluster_by)) "trial" else cs$cluster_by
co_unit        <- switch(co_by, trial = "trial", trial_in_review = "trial within review", review = "review", "cluster")
co_unit_plural <- switch(co_by, trial = "trials", trial_in_review = "trials within reviews", review = "reviews", "clusters")
year_note <- function(groups) {                              #how the period subsamples are defined, when any is shown
  if (!any(c("early", "late") %in% groups)) return("")
  n_of <- function(g) {
    r <- RES$results[RES$results$group == g, , drop = FALSE]
    if (nrow(r)) as.numeric(r$n[1]) else NA_real_
  }
  no_year <- n_of("all") - n_of("early") - n_of("late")
  paste0(sprintf("Published by %d and after %d: trials by year of publication (the study's year, else the year in the trial's label)",
                 as.integer(cs$year_split), as.integer(cs$year_split)),
         if (isTRUE(no_year > 0)) sprintf("; the %s outcomes without a year are in neither", format(no_year, big.mark = ",")) else "",
         ". ")
}

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
  sprintf("The projection test, Kudrin's tests, CS1 and CS2B cluster by article (Panel A) or %s (Panel B); the other EKW tests and the $p$-curve test treat the statistics as independent. ", co_unit),
  "Statistics whose rounding interval exceeds half a standard error are left out ($\\kappa \\le 0.5$); estimates printed as zero are kept. ",
  sprintf("Panel B: continuous efficacy outcomes of randomized trials with at least %d participants. ", as.integer(cs$n_min)),
  year_note(co_groups),
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

#### 6. application_MM_output.R's two rounding tables, with the Cochrane results below RDD ------------------
# The same tables as its table_rounding_ekw<mm_tag>.tex and table_rounding_appendix<mm_tag>.tex, built by the same
# rules from the same Brodeur results files (so the Brodeur blocks are identical to its), plus one block for each
# Cochrane subsample in rt_groups, from Cochrane_results<co_tag>.rds. Saved results only: no analysis is rerun.
if (rt_write) {
  rt_J0    <- 20                                                #the main specification; the appendix adds J0 = 30
  rt_treat <- c(none = "None", uniform = "Unif.", smin = "$s_{\\min}$", smin_noise = "$s_{\\min}$+noise")
  rt_treat <- rt_treat[intersect(names(rt_treat), MS$keys_run)] #the treatments the Brodeur run has
  rt_keys  <- names(rt_treat)
  rt_Brows <- intersect(c("smin", "smin_noise"), rt_keys)       #rows with B-hat: the treatments the allowance covers
  rt_tests <- list(main     = c("proj", "T_inf", "T_H", "Disc", "Binomial", "CS1", "CS2B", "LCM", "Fisher", "Stouffer"),
                   appendix = c("proj30", "proj_noshift"))                   #the projection test's variants
  rt_cat <- data.frame(
    family  = c("Projection test", "Projection, $J_0=30$", "Projection, no shift", "Kudrin (2025)", "Kudrin (2025)",
                "EKW", "EKW", "EKW", "EKW", "$p$-curve", "EKW", "EKW"),
    colhead = c("$p$", "$p$", "$p$", "$T_\\infty$", "$T_H$", "Disc.", "Binom.", "CS1", "CS2B", "Stouffer", "LCM", "Fisher"),
    short   = c("the projection test", "the projection test with $J_0=30$", "the unshifted projection test", "$T_\\infty$",
                "$T_H$", "Disc.\\@", "Binom.", "CS1", "CS2B", "the $p$-curve test", "LCM", "Fisher"),
    clust   = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE),
    has_B   = c(TRUE, TRUE, rep(FALSE, 10)),
    row.names = c("proj", "proj30", "proj_noshift", "T_inf", "T_H", "Disc", "Binomial", "CS1", "CS2B", "Stouffer", "LCM", "Fisher"),
    stringsAsFactors = FALSE)
  rt_cite <- c(brodeur = "Brodeur et al.\\ (2020)", ekw = "Elliott et al.\\ (2022)", kudrin = "Kudrin (2025)",
               coxshi = "Cox and Shi (2023)", simonsohn = "Simonsohn et al.\\ (2014)")
  if (!all(c(rt_J0, 30) %in% MS$J0_list) || !all(c(rt_J0, 30) %in% cs$J0_list))
    stop("Section 6 needs J0 = 20 and J0 = 30 in both runs", call. = FALSE)
  if (!same(MS$epsilon_U_pre, cs$epsilon_U_pre))
    warning("E(L, delta) differs: Brodeur ", MS$epsilon_U_pre, ", Cochrane ", cs$epsilon_U_pre, call. = FALSE)
  
  #--- Brodeur: every test under every treatment, from the files application_MM_output.R reads
  rt_num  <- function(x) vapply(x, function(v) as.numeric(v)[1], numeric(1), USE.NAMES = FALSE)   #results hold list-columns
  rt_read <- function(stem, key) {
    f <- file.path(mm_dir, paste0(stem, "_", key, mm_tag))
    if (!file.exists(f)) stop("Missing Brodeur results: ", basename(f), call. = FALSE)
    readRDS(f)
  }
  rt_pick <- function(d, what) {                                #rows in the order of meth
    i <- match(meth, d$methods)
    if (anyNA(i)) stop(what, " has no results for ", paste(meth[is.na(i)], collapse = ", "), call. = FALSE)
    d[i, ]
  }
  rt_mm <- lapply(setNames(rt_keys, rt_keys), function(key) {
    pr  <- rt_pick(rt_read(paste0("MM_weights_results_", rt_J0, "_J0"), key), "MM_weights_results_20_J0")
    p30 <- rt_pick(rt_read("MM_weights_results_30_J0", key), "MM_weights_results_30_J0")
    pns <- rt_pick(rt_read(paste0("MM_weights_results_", rt_J0, "_J0_noshift"), key), "MM_weights_results_20_J0_noshift")
    kd  <- rt_pick(rt_read("Kudrin_results", key), "Kudrin_results")
    st  <- rt_pick(rt_read("Stouffer_results", key), "Stouffer_results")
    ek  <- rt_read("EWK_results", key)                          #matrix: tests x designs
    ek  <- matrix(rt_num(ek), nrow(ek), dimnames = dimnames(ek))[, paste0("ewk_", tolower(meth)), drop = FALSE]
    er  <- function(r) if (r %in% rownames(ek)) unname(ek[r, ]) else stop("EWK_results_", key, " has no ", r, call. = FALSE)
    list(p = list(proj = rt_num(pr$pval), proj30 = rt_num(p30$pval), proj_noshift = rt_num(pns$pval),
                  T_inf = rt_num(kd$T_inf), T_H = rt_num(kd$T_H), Stouffer = rt_num(st$pval),
                  LCM = er("lcms_EWK"), Fisher = er("Fisher_EWK"), Disc = er("disconts_EWK"),
                  Binomial = er("binomial_EWK"), CS1 = er("CS1_EWK_s"), CS2B = er("CS2_EWK_s")),
         B = rt_num(pr$breakdown), B30 = rt_num(p30$breakdown),
         p_hyb = if (is.null(kd$T_H_hybrid)) rep(NA_real_, length(meth)) else rt_num(kd$T_H_hybrid),
         restricted = if (is.null(kd$restricted)) rep(NA, length(meth)) else as.logical(rt_num(kd$restricted)),
         n = rt_num(pr$n), articles = rt_num(pr$articles))
  })
  if (!isTRUE(all(rt_mm[[1]]$n == MS$n_by_design[meth])))
    warning("The Brodeur results files have n = ", paste(rt_mm[[1]]$n, collapse = "/"), " but the run's settings say ",
            paste(MS$n_by_design[meth], collapse = "/"), " (stale files?)", call. = FALSE)
  rt_allow_mm <- function(j0) {                                 #the run's rounding allowance by design, as in section 1
    if (!is.null(MS$allowance)) MS$allowance[[as.character(j0)]][meth] else MS$C_round[[as.character(j0)]] * MS$kappa6[meth]
  }
  rt_A20 <- as.numeric(rt_allow_mm(rt_J0))
  rt_A30 <- as.numeric(rt_allow_mm(30))
  rt_blocks <- lapply(seq_along(meth), function(i) list(
    title = meth[i], n = rt_mm[[1]]$n[i], clusters = rt_mm[[1]]$articles[i], what = "$t$-statistics", unit = "articles",
    allow = c("20" = rt_A20[i], "30" = rt_A30[i]),
    res = lapply(rt_mm, function(r) list(p = lapply(r$p, `[`, i), B = r$B[i], B30 = r$B30[i],
                                         p_hyb = r$p_hyb[i], restricted = r$restricted[i]))))
  
  #--- Cochrane: the same tests under the same treatments, one block per subsample in rt_groups
  rt_co_title <- c(all   = "Cochrane, all trials",
                   early = sprintf("Cochrane, published by %d", as.integer(cs$year_split)),
                   late  = sprintf("Cochrane, published after %d", as.integer(cs$year_split)))
  rt_groups   <- intersect(intersect(rt_groups, names(rt_co_title)), unique(RES$results$group))
  rt_first_co <- if (length(rt_groups)) length(rt_blocks) + 1 else NA_integer_
  rt_col <- function(r, nm) if (nrow(r) && nm %in% names(r)) as.numeric(r[[nm]][1]) else NA_real_
  for (g in rt_groups) {
    rows <- RES$results[RES$results$group == g, , drop = FALSE]
    rt_blocks[[length(rt_blocks) + 1]] <- list(
      title = rt_co_title[[g]], n = rt_col(rows, "n"), clusters = rt_col(rows, "clusters"), what = "outcomes",
      unit = co_unit_plural, allow = c("20" = rt_col(rows, paste0("allowance_J0_", rt_J0)), "30" = rt_col(rows, "allowance_J0_30")),
      res = lapply(setNames(rt_keys, rt_keys), function(key) {
        r <- rows[rows$treatment == key, , drop = FALSE]
        if (!nrow(r)) return(NULL)                              #a treatment the Cochrane run lacks: its cells are dots
        list(p = list(proj = rt_col(r, paste0("proj_p_J0_", rt_J0)), proj30 = rt_col(r, "proj_p_J0_30"),
                      proj_noshift = rt_col(r, "proj_p_noshift"), T_inf = rt_col(r, "T_inf"), T_H = rt_col(r, "T_H"),
                      Stouffer = rt_col(r, "Stouffer"), LCM = rt_col(r, "LCM"), Fisher = rt_col(r, "Fisher"),
                      Disc = rt_col(r, "Disc"), Binomial = rt_col(r, "Binomial"), CS1 = rt_col(r, "CS1"), CS2B = rt_col(r, "CS2B")),
             B = rt_col(r, paste0("proj_B_J0_", rt_J0)), B30 = rt_col(r, "proj_B_J0_30"), p_hyb = rt_col(r, "T_H_hybrid"),
             restricted = as.logical(rt_col(r, "hybrid_restricted")))   #the share of draws in which it bound
      }))
  }
  
  #--- Cells, exactly as application_MM_output.R formats them
  rt_fmt_p <- function(x) ifelse(is.na(x), "$\\cdot$", sub("^0\\.", ".", sprintf("%.2f", x)))
  rt_fmt_d <- function(x) ifelse(is.na(x), "$\\cdot$",
                                 ifelse(x == 0, "0", ifelse(1e3 * x < 0.1, "$<0.1$", sprintf("%.1f", 1e3 * x))))
  rt_big   <- function(x) formatC(x, big.mark = ",", format = "d")
  rt_Bof   <- function(tk, rr) if (tk == "proj30") rr$B30 else rr$B                   #B-hat of a projection test
  rt_Aof   <- function(tk, blk) if (tk == "proj30") blk$allow[["30"]] else blk$allow[["20"]]   #and its allowance
  rt_hyb_same <- function(rr) {                                 #does Kudrin's hybrid test coincide with T_H here?
    isFALSE(rr$restricted) || (!is.na(rr$p_hyb) && !is.na(rr$p$T_H) && abs(rr$p_hyb - rr$p$T_H) < 1e-12)
  }
  rt_cell_p <- function(tk, key, blk) {
    rr <- blk$res[[key]]
    if (is.null(rr)) return("$\\cdot$")
    p <- rr$p[[tk]]
    s <- rt_fmt_p(p)
    certified <- !rt_cat[tk, "has_B"] || !(key %in% rt_Brows) || isTRUE(rt_Bof(tk, rr) > rt_Aof(tk, blk))
    if (rt_bold && !is.na(p) && p < 0.05 && certified) s <- paste0("\\textbf{", s, "}")
    if (tk == "T_H") {
      s <- if (rt_hyb_same(rr)) paste0(s, "$^{\\dagger}$") else if (!is.na(rr$p_hyb)) paste0(s, " [", rt_fmt_p(rr$p_hyb), "]") else s
    }
    s
  }
  rt_cell_B <- function(tk, key, blk) {
    rr <- blk$res[[key]]
    if (!(key %in% rt_Brows)) "" else if (is.null(rr)) "$\\cdot$" else rt_fmt_d(rt_Bof(tk, rr))
  }
  
  #--- LaTeX: one block per design or subsample, one row per treatment; a rule between the two applications
  rt_build <- function(tests, caption, label, notes) {
    fam    <- unlist(lapply(tests, function(tk) rep(rt_cat[tk, "family"], 1 + rt_cat[tk, "has_B"])))
    heads  <- unlist(lapply(tests, function(tk) c(rt_cat[tk, "colhead"], if (rt_cat[tk, "has_B"]) "$10^3\\widehat{B}$")))
    ncols  <- 1 + length(heads)
    runs   <- rle(fam)
    ends   <- 1 + cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1
    b_tk   <- tests[rt_cat[tests, "has_B"]][1]                  #the projection test whose allowance the block line gives
    body <- unlist(lapply(seq_along(rt_blocks), function(b) {
      blk  <- rt_blocks[[b]]
      line <- sprintf("\\textit{%s}: %s %s in %s %s%s", blk$title, rt_big(blk$n), blk$what, rt_big(blk$clusters), blk$unit,
                      if (!is.na(b_tk)) paste0("; rounding allowance ", rt_fmt_d(rt_Aof(b_tk, blk))) else "")
      rows <- vapply(rt_keys, function(key) paste0(rt_treat[[key]], " & ", paste(unlist(lapply(tests, function(tk)
        c(rt_cell_p(tk, key, blk), if (rt_cat[tk, "has_B"]) rt_cell_B(tk, key, blk)))), collapse = " & "), " \\\\"), "")
      c(if (isTRUE(b == rt_first_co)) "\\midrule" else if (b > 1) "\\addlinespace",
        sprintf("\\multicolumn{%d}{l}{%s} \\\\", ncols, line), rows)
    }))
    c("\\begin{table}[htbp]", "\\centering", paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
      "\\small",
      "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
      sprintf("\\begin{tabular}{l%s}", strrep("c", ncols - 1)), "\\toprule",
      paste0(" & ", paste(sprintf("\\multicolumn{%d}{c}{%s}", runs$lengths, runs$values), collapse = " & "), " \\\\"),
      paste(sprintf("\\cmidrule(lr){%d-%d}", starts, ends), collapse = " "),
      paste0("Treatment & ", paste(heads, collapse = " & "), " \\\\"),
      "\\midrule", body, "\\bottomrule", "\\end{tabular}}", "", "\\medskip",
      "\\begin{minipage}{\\linewidth}", "\\footnotesize", notes, "\\end{minipage}", "\\end{table}")
  }
  
  #--- Notes: short; the appendix table's define only what is new there
  rt_and <- function(x) if (length(x) < 2) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
  rt_gl  <- function(lbl) if (substr(lbl, 1, 1) == "$") lbl else paste0("\\textit{", lbl, "}")
  rt_def <- c(none       = "printed statistics",
              uniform    = "every printed number drawn uniformly within its rounding interval",
              smin       = "the smallest standard error consistent with the printed numbers (Section 5.2)",
              smin_noise = "the same plus uniform noise at the estimate's last digit")
  rt_lab   <- vapply(rt_keys, function(k) rt_gl(rt_treat[[k]]), "")
  rt_noisy <- intersect(c("uniform", "smin_noise"), rt_keys)
  rt_open  <- paste0("\\textit{Notes:} $p$-values. ", rt_cite[["brodeur"]], ": $t$-statistics by research design. Cochrane: ",
                     sprintf("$z$-statistics of continuous efficacy outcomes of randomized trials with at least %d participants.", as.integer(cs$n_min)),
                     if (any(c("early", "late") %in% rt_groups))
                       sprintf(" Its periods split trials by year of publication at %d.", as.integer(cs$year_split)))
  rt_rows  <- paste0("Rows: ", paste(paste0(rt_lab, ", ", rt_def[rt_keys]), collapse = "; "), ".",
                     if (length(rt_noisy)) paste0(" ", rt_and(rt_lab[rt_noisy]), " average over ",
                                                  paste(unique(as.integer(MS$n_draws[rt_noisy])), collapse = " and "), " draws",
                                                  if (!same(cs$n_draws, MS$n_draws_noise)) sprintf(" (Cochrane: %d)", as.integer(cs$n_draws)), "."))
  rt_excl  <- sprintf("Statistics with $\\kappa>%s$ are left out (Section 5.2).", format(MS$kappa_cap))
  rt_projn <- sprintf(paste0("Projection test: $J_0=%d$, $J=%d$, centered at 1.96. $10^3\\widehat{B}$: a lower ",
                             "95\\%% confidence bound on the distance to the nearest honest $t$-curve; the rounding allowance (first line ",
                             "of each block), $%s\\,\\overline{\\kappa^6}$, bounds what rounding can create (Propositions 3 and 4)."),
                      as.integer(rt_J0), as.integer(MS$J), format(MS$C_round[[as.character(rt_J0)]]))
  rt_any_hyb <- any(vapply(rt_blocks, function(blk) any(vapply(blk$res, function(rr) !is.null(rr) && !rt_hyb_same(rr), TRUE)), TRUE))
  rt_kudn  <- paste0("$T_\\infty$, $T_H$: ", rt_cite[["kudrin"]], "; $^\\dagger$: his hybrid test gives the same $p$-value",
                     if (rt_any_hyb) " (otherwise it follows in brackets)", ".")
  rt_ekwn  <- paste0("EKW, as in ", rt_cite[["ekw"]], ": density discontinuity at .05; binomial on $[.04,.05]$; on $(0,.15]$, Cox--Shi ",
                     "monotonicity (CS1) and 2-monotonicity with bounds (CS2B), and LCM and Fisher tests of a concave CDF.")
  rt_pcn   <- paste0("Stouffer: the $p$-curve test for left skew of ", rt_cite[["simonsohn"]], ".")
  rt_cln   <- paste0("The projection test, $T_\\infty$, $T_H$, CS1 and CS2B cluster by article or ", co_unit, "; the others do not.")
  rt_boldn <- if (rt_bold) "Bold: $p<.05$; in the $s_{\\min}$ rows, the projection test only when $\\widehat{B}$ exceeds the allowance."
  rt_ns <- if (is.null(MS$noshift_n_draws)) NULL else MS$noshift_n_draws[rt_keys]   #draws behind the unshifted test (Brodeur)
  rt_ns_fewer <- if (is.null(rt_ns)) logical(0) else (!is.na(rt_ns) & rt_ns != MS$n_draws[rt_keys])
  rt_ns_note <- if (any(rt_ns_fewer)) paste0("; under ", rt_and(unname(rt_treat[rt_keys[rt_ns_fewer]])), " its Brodeur $p$-values average over ",
                                             rt_and(as.character(rt_ns[rt_ns_fewer])), " draws") else ""
  rt_time_file <- file.path(mm_dir, paste0("Run_times", mm_tag))
  rt_runtime <- if (file.exists(rt_time_file)) {
    rt_rt <- readRDS(rt_time_file)
    rt_rr <- function(x) if (diff(range(x)) < 1) sprintf("%.1f", max(x)) else sprintf("%.0f--%.0f", min(x), max(x))
    sprintf("Run time per design of %s, on one core: %s s for the projection test, %s s for Kudrin's tests, %s s for the EKW tests.",
            rt_cite[["brodeur"]], rt_rr(rt_rt[, "projection"]), rt_rr(rt_rt[, "kudrin"]), rt_rr(rt_rt[, "ekw"]))
  } else NULL
  rt_notes <- list(
    main = c(rt_open, rt_rows, rt_excl, rt_projn, rt_kudn, rt_ekwn, rt_pcn, rt_cln, rt_boldn,
             paste0("Projection tests with $J_0=30$ and without the shift: Table~\\ref{", rt_label[["appendix"]], "}.")),
    appendix = c(
      paste0("\\textit{Notes:} $p$-values for the samples and treatments of Table~\\ref{", rt_label[["main"]], "}."),
      sprintf(paste0("$J_0=30$: weights $w_j=\\min\\{1,2^{30-j}\\}$; $10^3\\widehat{B}$ and the rounding allowance (first line of each ",
                     "block, $%s\\,\\overline{\\kappa^6}$) are for $J_0=30$."), format(MS$C_round[["30"]])),
      paste0("No shift: $t$-statistics not shifted by 1.96 (Remark 13); no allowance is computed for it, so its rejections are not ",
             "certified against rounding", rt_ns_note, "."),
      sprintf("Both cluster by article or %s and use %s bootstrap repetitions.", co_unit, rt_big(MS$boots)),
      if (rt_bold) "Bold: $p<.05$; in the $s_{\\min}$ rows, the $J_0=30$ test only when $\\widehat{B}$ exceeds the allowance.",
      rt_runtime))
  rt_count   <- c("one", "two", "three", "four")[length(rt_keys)]
  rt_zero    <- if (isTRUE(MS$drop_zero_mu)) "" else " (keeping coefficients reported as zero)"
  rt_caption <- c(main     = paste0("Tests for selective reporting under ", rt_count, " treatments of rounding, in economics", rt_zero, " and medicine"),
                  appendix = paste0("The projection test with $J_0=30$ and without the shift, under ", rt_count,
                                    " treatments of rounding, in economics", rt_zero, " and medicine"))
  for (nm in names(rt_tests)) {
    writeLines(rt_build(rt_tests[[nm]], rt_caption[[nm]], rt_label[[nm]], rt_notes[[nm]]), file.path(out_dir, rt_file[[nm]]))
    cat("Rounding table with Cochrane written to ", file.path(out_dir, rt_file[[nm]]), "\n", sep = "")
  }
  
  #--- The Cochrane blocks on the console (the Brodeur blocks are application_MM_output.R's)
  rt_all_tests <- unlist(rt_tests, use.names = FALSE)
  for (blk in rt_blocks[if (is.na(rt_first_co)) integer(0) else seq(rt_first_co, length(rt_blocks))]) {
    M <- t(vapply(rt_keys, function(key) vapply(rt_all_tests, function(tk) {
      rr <- blk$res[[key]]
      if (is.null(rr)) NA_real_ else rr$p[[tk]]
    }, 0), numeric(length(rt_all_tests))))
    cat("\n", blk$title, ": p-values by treatment\n", sep = "")
    print(round(M, 3))
  }
}

#### 7. The main table without Kudrin's tests ---------------------------------------------------------------
# Section 4's table without its three Kudrin (2025) columns (T_inf, T_H, T_H^hyb), and without the two places its
# notes mention them: the sentence on T_inf's rejections in Panel B, and Kudrin's tests among the tests that cluster.
# Everything else is section 4's: the same rows, bold rules, certified column and notes (whose citation of Elliott,
# Kudrin and Wuthrich stays, as it is the EKW tests'). Saved results only. Writes the table to nk_file.
nk_file  <- file.path(out_dir, "table_main_applications_nokudrin.tex")
nk_label <- "tab:main_applications_nokudrin"                 #"tab:main_applications" if it replaces section 4's table in the paper

nk_pair <- function(label, size, pr, de, cert) {              #pair() without Kudrin's three columns
  c(sprintf("%s & Printed & %s & %s & \\multirow{2}{*}{%s} \\\\", label, fp(pr$proj), fp(pr$ekw_st), cert),
    sprintf("{\\footnotesize %s} & De-rounded & %s & %s & \\\\", size, fp(de$proj), fp(de$ekw_st)))
}
nk_body_mm <- unlist(lapply(seq_along(meth), function(i) {
  c(if (i > 1) "\\addlinespace[2pt]",
    nk_pair(meth[i], mm_size[i], mm_printed[i, ], mm_derounded[i, ], fB(mm_cert$B[i], mm_allow[[i]])))
}))
nk_body_co <- unlist(lapply(seq_along(co_groups), function(i) {
  c(if (i > 1) "\\addlinespace[2pt]",
    nk_pair(co_label[[i]], co_size[i], co_printed[i, ], co_derounded[i, ], fB(co_cert$B[i], co_allow[i])))
}))

nk_notes <- sub(sprintf("$T_\\infty$'s rejections in Panel B come from $p$-values near 1: outcomes whose two arms have equal printed means (%s). ",
                        ref_fig_tinf), "", notes, fixed = TRUE)     #section 4's notes, less the sentence on T_inf
nk_notes <- sub("Kudrin's tests, ", "", nk_notes, fixed = TRUE)    #and less Kudrin's tests among those that cluster
if (grepl("T_\\\\infty|T_H|Kudrin's|Kudrin \\(2025\\)", nk_notes)) {
  stop("Section 4's notes still mention Kudrin's tests after the cuts: their wording has changed, so update the two sub() calls",
       call. = FALSE)
}

nk_tex <- c(sprintf("%% Written by main_applications_table.R, section 7 (no Kudrin columns), on %s from Settings_application%s and Cochrane_results%s.rds",
                    format(Sys.time(), "%Y-%m-%d %H:%M"), mm_tag, co_tag),
            "% Needs \\usepackage{booktabs, threeparttable, multirow}",
            "\\begin{table}[!htbp]\\centering",
            "\\begin{threeparttable}",
            "\\caption{Tests for selective reporting in economics and medicine, on printed and de-rounded statistics}",
            sprintf("\\label{%s}", nk_label),
            "\\begin{tabular}{llccc}",
            "\\toprule",
            " & & Projection & EKW and & Projection test \\\\",
            "Sample & Statistics & test & $p$-curve & certified: $10^3 \\hat B$ \\\\",
            " & & $p$ & smallest $p$ & (allowance) \\\\",
            "\\midrule",
            "\\multicolumn{5}{l}{\\textit{A. Economics: Brodeur, Cook and Heyes (2020), by research design}} \\\\",
            nk_body_mm,
            "\\midrule",
            "\\multicolumn{5}{l}{\\textit{B. Medicine: Cochrane Database of Systematic Reviews, continuous efficacy outcomes}} \\\\",
            nk_body_co,
            "\\bottomrule",
            "\\end{tabular}",
            "\\begin{tablenotes}[flushleft]\\footnotesize",
            paste0("\\item ", nk_notes),
            "\\end{tablenotes}",
            "\\end{threeparttable}",
            "\\end{table}")
writeLines(nk_tex, nk_file)
cat("Table without Kudrin's tests written to ", nk_file, "\n", sep = "")