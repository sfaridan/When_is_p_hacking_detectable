#### =====================================================================================================
#### application_MM_output.R: every table and figure of the application, from application_MM_analysis.R
#### Runs no test and takes seconds. It reads the settings and results that application_MM_analysis.R saved in
#### <root>/results for one run, and writes, tagged like that run:
####   table_application_<option><tag>.tex   per de-rounding option: every test, J0 = 20 and 30
####   table_rounding_ekw<tag>.tex             Table 3: one block per design, one row per treatment of rounding
####   table_rounding_appendix<tag>.tex        its companion: J0 = 30, no shift, Stouffer, LCM and Fisher
####   Figure_residual_tcurve_<option><tag>    the part of the reported t-curve no honest t-curve explains
####   Figure_tcurve_shape_<option><tag>       reported |t| with the closest honest t-curve and the curvature check
#### Every setting in the notes (sample, kappa cap, draws, bootstrap repetitions, grid, E(L, delta), allowances)
#### comes from the run's settings file, so the notes always describe the run behind the numbers. Below are
#### only the choice of run and choices of presentation.
#### =====================================================================================================
rm(list = ls())
root    <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
res_dir <- file.path(root, "results")
run_tag <- "_keepmu0"   #the run to report: "_keepmu0" keeps zero coefficients (application_MM_analysis.R's default),
#"" is a run with Kudrin's (2025) drop of them; plus any run_label the run had

#--- Presentation
bold               <- TRUE   #Table 3: bold p < .05; under smin rows a projection test only when its B-hat exceeds the allowance
p_digits           <- 2
show_runtime       <- TRUE   #run-time sentence in the appendix table's notes, if the run timed the tests
fig_show_pointwise <- FALSE  #reported-versus-honest figure: also shade where the curvature check fails pointwise
fig_show_honest    <- TRUE   #and draw the closest honest t-curve
cite <- c(brodeur = "Brodeur et al.\\ (2020)", ekw = "Elliott et al.\\ (2022)", kudrin = "Kudrin (2025)",
          coxshi  = "Cox and Shi (2023)", simonsohn = "Simonsohn et al.\\ (2014)")   #Table 3's; or "\\citet{key}"

#### The run's settings ---------------------------------------------------------------------------------
settings_file <- file.path(res_dir, paste0("Settings_application", run_tag))
if (!file.exists(settings_file)) stop("No ", basename(settings_file), ": run application_MM_analysis.R with this run_tag first")
S <- readRDS(settings_file)
stopifnot(identical(S$run_tag, run_tag))
if (!identical(S$status, "complete")) warning("The run with run_tag '", run_tag, "' did not finish: its settings say '", S$status, "'")
if (!identical(S$tests_from, "this run")) message("Test results from ", S$tests_from)
cat(sprintf("Reporting the run saved %s: %s zero coefficients, kappa cap %s, %d de-roundings under options 2 and 4, %d bootstrap repetitions\n",
            S$saved, if (S$drop_zero_mu) "drops" else "keeps", S$kappa_cap, S$n_draws_noise, S$boots))
meth          <- c("RCT", "IV", "DID", "RDD")
drop_zero_mu  <- S$drop_zero_mu    #the names the table code below uses, all from the run
kappa_cap     <- S$kappa_cap
se_rounded_up <- S$se_rounded_up
J             <- S$J
coeffs        <- S$coeffs
L             <- S$L
numgrid       <- S$numgrid
boots         <- S$boots
J0_list       <- S$J0_list
J0            <- 20                #the main specification
epsilon_U_pre <- S$epsilon_U_pre   #E(L, delta)
E_LD          <- S$epsilon_U_pre
C_round       <- S$C_round

#### 1. One table per de-rounding option: every test, J0 = 20 and 30 -----------------------------------
# Reads the results files saved below for one de-rounding option (so after a full run you can tweak and
# re-run this block alone), prints the table and writes results/table_application<suffix>.tex.
# LaTeX needs booktabs, graphicx (shrink-to-fit wrapper) and hyperref or url (for \url in the notes).

make_application_table <- function(J0_values = c(20),            # one (p, B-hat) group per J0, in this order; each must be in J0_list above
                                   n_draws   = 1,                # number of de-roundings the results files average over
                                   spec,                         # the de-rounding option: list(key, label)
                                   suffix    = paste0("_", spec$key, run_tag),
                                   res_dir   = paste0(root, "/results"),
                                   out_file  = paste0("table_application", suffix, ".tex")) {
  meth <- c("RCT", "IV", "DID", "RDD")
  
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
  tab <- data.frame(n = proj[[1]]$n, articles = proj[[1]]$articles, row.names = meth)
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
  cells <- cbind(meth, sprintf("%.0f", tab$n), sprintf("%.0f", tab$articles))
  for (j0 in J0_values) cells <- cbind(cells, fmt(tab[[paste0("p_J0_", j0)]], 2), fmt(tab[[paste0("B_J0_", j0)]], 4))
  for (k in c(names(ekw_rows), names(kud_cols))) cells <- cbind(cells, fmt(tab[[k]], 2))
  unr <- !tab$hybrid_restricted                                  #dagger: hybrid restriction did not bind
  if (any(unr)) cells[unr, ncol(cells) - 1] <- paste0(cells[unr, ncol(cells) - 1], "$^{\\dagger}$")
  
  #--- Notes (one sentence per line in the .tex; edit freely)
  deround_notes <- c(
    none       = "De-rounding: none. Each $t$-statistic is the ratio of the reported coefficient and standard error, and reported $t$-statistics are used as printed.",
    uniform    = paste0("De-rounding: uniform, following ", cite[["brodeur16"]], ". The coefficient and the standard error are each replaced by an independent uniform draw from its rounding interval, ",
                        "and a reported $t$-statistic gets uniform noise at its last digit."),
    smin       = "De-rounding: smallest consistent standard error (Section 5.2). Each coefficient is divided by the smallest standard error consistent with the reported one, and reported $t$-statistics are used as printed.",
    smin_noise = paste0("De-rounding: smallest consistent standard error with uniform noise (Section 5.2). Each coefficient, plus uniform noise at its last digit, is divided by the smallest standard error ",
                        "consistent with the reported one, and a reported $t$-statistic gets uniform noise at its last digit."))
  notes <- c(
    paste0("Results for the ", cite[["brodeur"]], " data. ",
           if (drop_zero_mu) paste0("Following ", cite[["kudrin"]], ", we drop tests whose reported coefficient is missing or zero.") else
             paste0("We drop tests whose reported coefficient is missing; unlike ", cite[["kudrin"]], ", we keep those reported as zero.")),
    deround_notes[[spec$key]],
    sprintf("The sample is the same for every de-rounding choice: we also drop $t$-statistics whose standard error is reported as zero or whose coefficient's rounding unit exceeds %s times the smallest standard error consistent with the reported one.", kappa_cap),
    if (n_draws > 1) sprintf("All $p$-values and breakdown statistics are averages over %d independent de-roundings.", n_draws),
    paste0("Columns headed $J_0$ report $p$-values and breakdown statistics $\\widehat{B}$ of the projection test with $J=", J, "$",
           if (length(J0_values) > 1) "; several values of $J_0$ are included to show low sensitivity to the smoothing parameter", "."),
    sprintf("We form $\\mathcal{X}$ using a grid of %d elements evenly spaced between $\\pm %s$ plus the scalar 9999.", numgrid, L),
    sprintf("Discretization incurs error $\\mathcal{E}(L,\\delta)=%.5f$, which is already accounted for in $\\widehat{B}$ and the $p$-value.", epsilon_U_pre),
    "The discretization was chosen to be made of round numbers and to give an error that is small relative to the critical values and test statistics.",
    paste0("EKW $p$-values are for the tests of ", cite[["ekw"]], ": the least concave majorant test (LCM) and Fisher's test (Fisher) on $p$-values in $(0,0.15]$, ",
           "the density discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$ (Binomial), and the ", cite[["coxshi"]],
           " histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) and for 2-monotonicity with bounds (CS2B)."),
    paste0(cite[["kudrin"]], " $p$-values are for his KS-type test ($T_\\infty$), histogram-based test ($T_H$) and hybrid histogram-based test ($T_H^{\\mathrm{hybrid}}$), ",
           "using the sinc kernel with frequency cutoff 3 and, for $T_H$ and $T_H^{\\mathrm{hybrid}}$, the points $0.01,0.02,\\ldots,0.15$."),
    "All three cluster by article: $T_H$ and $T_H^{\\mathrm{hybrid}}$ use cluster-robust covariance matrices and $T_\\infty$ uses an article-level bootstrap.",
    if (any(unr)) "$^\\dagger$The deconvolution estimate already satisfies the hybrid test's shape restrictions, so $T_H^{\\mathrm{hybrid}}=T_H$.",
    if (n_draws > 1) "Min is the smallest of the averaged $p$-values in each group of tests." else "Min is the smallest $p$-value in each group of tests.",
    sprintf("The projection test and $T_\\infty$ use %d bootstrap repetitions.", boots),
    paste0("Computation used osqp package in R by ", cite[["osqp"]], ". Data: ", cite[["brodeur"]], "."),
    "Transparency R package available at: \\url{https://github.com/sfaridan/When_is_p_hacking_detectable}")
  
  #--- LaTeX
  width <- c(rep(2, length(J0_values)), 7, 4)                    #column groups after the 3 leading columns
  last  <- 3 + cumsum(width)
  grp   <- c(sprintf("$J_0=%s$", J0_values), "EKW $p$-values", paste(cite[["kudrin"]], "$p$-values"))
  head2 <- c("", "$n$", "Articles", rep(c("$p$", "$\\widehat{B}$"), length(J0_values)),
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
           tex_row(c(rep("", 3), paste0("\\multicolumn{", width, "}{c}{", grp, "}"))),
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
  dimnames(txt) <- list(meth, c("n", "Articles", paste0(c("p_J0=", "B_J0="), rep(J0_values, each = 2)),
                                "LCM", "Fisher", "Disc.", "Binomial", "CS1", "CS2B", "EKW_min", "T_inf", "T_H", "T_H_hyb", "Kud_min"))
  cat("\nDe-rounding: ", spec$label, "\n", sep = "")
  print(noquote(txt), right = TRUE)
  if (any(unr)) cat("+ hybrid restriction did not bind, so T_H_hybrid = T_H\n")
  cat("LaTeX table written to", file.path(res_dir, out_file), "\n")
  invisible(list(numbers = tab, latex = tex))
}

for (key in S$keys_run) {
  make_application_table(J0_values = J0_list, n_draws = S$n_draws[[key]],
                         spec = list(key = key, label = S$deround_labels[match(key, S$deround_keys)]))
}

#### 2. Table 3 and its appendix companion: one block per design, one row per treatment of rounding -----
#--- Layout
treatments   <- c(none = "None", uniform = "Unif.", smin = "$s_{\\min}$", smin_noise = "$s_{\\min}$+noise")
B_rows       <- c("smin", "smin_noise")   #rows that show B-hat: the treatments the rounding allowance covers
tables_spec  <- list(
  main     = list(tests   = c("proj", "T_inf", "T_H", "Disc", "Binomial", "CS1", "CS2B"),
                  file    = paste0("table_rounding_ekw", run_tag, ".tex"), label = paste0("tab:rounding_ekw", run_tag),
                  caption = paste0("Tests for selective reporting under four treatments of rounding",
                                   if (drop_zero_mu) "" else ", keeping coefficients reported as zero")),
  appendix = list(tests   = c("proj30", "proj_noshift", "Stouffer", "LCM", "Fisher"),
                  file    = paste0("table_rounding_appendix", run_tag, ".tex"), label = paste0("tab:rounding_appendix", run_tag),
                  caption = paste0("Further tests for selective reporting under four treatments of rounding",
                                   if (drop_zero_mu) "" else ", keeping coefficients reported as zero")))
treatments <- treatments[intersect(names(treatments), S$keys_run)]   #the treatments this run has

#### Setup
meth      <- c("RCT", "IV", "DID", "RDD")
test_keys <- c("proj", "proj30", "proj_noshift", "T_inf", "T_H", "T_H_hyb",
               "LCM", "Fisher", "Disc", "Binomial", "CS1", "CS2B", "EKW_min", "Stouffer")
catalog <- data.frame(
  family  = c("Projection test", "Projection, $J_0=30$", "Projection, no shift",
              rep(cite[["kudrin"]], 3), rep("EKW", 7), "$p$-curve"),
  colhead = c("$p$", "$p$", "$p$", "$T_\\infty$", "$T_H$", "$T_H^{\\mathrm{hyb}}$",
              "LCM", "Fisher", "Disc.", "Binom.", "CS1", "CS2B", "Min", "Stouffer"),
  short   = c("the projection test", "the projection test with $J_0=30$", "the unshifted projection test",
              "$T_\\infty$", "$T_H$", "$T_H^{\\mathrm{hyb}}$", "LCM", "Fisher", "Disc.\\@", "Binom.", "CS1", "CS2B",
              "EKW Min", "the $p$-curve test"),
  plain   = c("Proj", "Proj30", "ProjNS", "Tinf", "TH", "THhyb", "LCM", "Fisher", "Disc", "Binom", "CS1", "CS2B",
              "EKWmin", "Stouffer"),
  src     = c("proj", "proj30", "projns", rep("kud", 3), rep("ewk", 7), "stouffer"),
  field   = c("pval", "pval", "pval", "T_inf", "T_H", "T_H_hybrid", "lcms_EWK", "Fisher_EWK", "disconts_EWK",
              "binomial_EWK", "CS1_EWK_s", "CS2_EWK_s", "min", "pval"),
  clust   = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, NA, FALSE),
  has_B   = c(TRUE, TRUE, rep(FALSE, 12)),     #projection tests with a certified B-hat column
  row.names = test_keys, stringsAsFactors = FALSE)
opt_of      <- c(none = 1, uniform = 2, smin = 3, smin_noise = 4)
averaged    <- c("uniform", "smin_noise")
treat_plain <- c(none = "None", uniform = "Unif", smin = "smin", smin_noise = "smin+U")
treat_keys  <- names(treatments)
B_rows      <- intersect(B_rows, treat_keys)
shown       <- unique(unlist(lapply(tables_spec, `[[`, "tests")))
stopifnot(all(treat_keys %in% names(opt_of)), all(shown %in% test_keys), as.character(J0) %in% names(C_round),
          all(vapply(tables_spec, function(s) !anyDuplicated(s$tests), TRUE)))

num  <- function(x) vapply(x, function(v) as.numeric(v)[1], numeric(1), USE.NAMES = FALSE)
pick <- function(d, what) {
  i <- match(meth, d$methods)
  if (anyNA(i)) stop(what, " has no results for ", paste(meth[is.na(i)], collapse = ", "), call. = FALSE)
  d[i, ]
}
res_file <- function(stem, key) {
  f <- file.path(res_dir, paste0(stem, "_", key, run_tag))
  if (!file.exists(f)) stop("No results file ", basename(f), ": run application_MM_analysis.R with option ",
                            opt_of[[key]], " first", call. = FALSE)
  f
}

#### Read results
read_results <- function(key) {
  proj_name <- paste0("MM_weights_results_", J0, "_J0")
  proj <- pick(readRDS(res_file(proj_name, key)), paste0(proj_name, "_", key))
  kud  <- pick(readRDS(res_file("Kudrin_results", key)), paste0("Kudrin_results_", key))
  ewk  <- readRDS(res_file("EWK_results", key))
  ewk  <- matrix(num(ewk), nrow(ewk), dimnames = dimnames(ewk))[, paste0("ewk_", tolower(meth)), drop = FALSE]
  extra  <- function(stem) pick(readRDS(res_file(stem, key)), paste0(stem, "_", key, run_tag))
  proj30 <- if ("proj30" %in% shown) extra("MM_weights_results_30_J0") else NULL
  projns <- if ("proj_noshift" %in% shown) extra(paste0("MM_weights_results_", J0, "_J0_noshift")) else NULL
  stou   <- if ("Stouffer" %in% shown) extra("Stouffer_results") else NULL
  if (!is.null(proj$numcoeffs) && any(num(proj$numcoeffs) != J + 1))
    warning(proj_name, "_", key, " was run with numcoeffs != ", J + 1, call. = FALSE)
  p <- lapply(setNames(shown, shown), function(k) {
    f <- catalog[k, "field"]
    switch(catalog[k, "src"],
           proj     = num(proj$pval),
           proj30   = num(proj30$pval),
           projns   = num(projns$pval),
           stouffer = num(stou$pval),
           kud      = if (f %in% names(kud)) num(kud[[f]]) else stop("Kudrin_results_", key, " has no ", f, call. = FALSE),
           ewk      = if (f %in% rownames(ewk)) unname(ewk[f, ]) else stop("EWK_results_", key, " has no ", f, call. = FALSE))
  })
  list(p = p, B = num(proj$breakdown), B30 = if (is.null(proj30)) NULL else num(proj30$breakdown),
       p_hyb = if ("T_H_hybrid" %in% names(kud)) num(kud$T_H_hybrid) else rep(NA_real_, length(meth)),
       restricted = if (is.null(kud$restricted)) rep(NA, length(meth)) else as.logical(num(kud$restricted)),
       n_proj = num(proj$n), n_kud = num(kud$n), articles = num(proj$articles), articles_kud = num(kud$articles))
}
res <- setNames(lapply(treat_keys, read_results), treat_keys)

count_draws <- function(key) {
  from_file <- function(stem, f) {
    path <- file.path(res_dir, paste0(stem, "_", key, run_tag))
    if (file.exists(path)) as.integer(f(readRDS(path))) else NA_integer_
  }
  n <- c(projection = from_file("MM_weights_results_draws", function(d) {
    if (!is.null(d$J0)) d <- d[d$J0 == J0, ]
    length(unique(d$draw)) }),
    Kudrin = from_file("Kudrin_results_draws", function(d) length(unique(d$draw))),
    EKW    = from_file("EWK_results_draws", function(a) if (length(dim(a)) == 3) dim(a)[3] else 1L))
  if (length(unique(n[!is.na(n)])) > 1)
    warning("'", key, "' files average over different numbers of draws (",
            paste(names(n), n, sep = " = ", collapse = ", "), ")", call. = FALSE)
  if (all(is.na(n))) NA_integer_ else unname(n[!is.na(n)][1])
}
n_draws <- vapply(treat_keys, count_draws, integer(1))
if (!isTRUE(all(n_draws == S$n_draws[names(n_draws)])))
  warning("The draws files hold ", paste(n_draws, collapse = "/"), " de-roundings, but the settings say ",
          paste(S$n_draws[names(n_draws)], collapse = "/"), " (stale files?)", call. = FALSE)

#### Rounding allowance: C times the mean of kappa^6, from the run itself (its kappa, on the sample it tested)
n_kept      <- as.vector(S$n_by_design[meth])
allowance   <- as.vector(S$allowance[[as.character(J0)]][meth])
allowance30 <- as.vector(S$allowance[["30"]][meth])                   #for the J0 = 30 test
n_des <- res[[1]]$n_proj                                   #t-statistics and articles per design, as tested
a_des <- res[[1]]$articles
for (key in treat_keys) {
  if (!isTRUE(all(res[[key]]$articles == a_des)) || !isTRUE(all(res[[key]]$articles_kud == a_des)))
    warning("Article counts in the '", key, "' results differ from the '", treat_keys[1], "' results", call. = FALSE)
}
if (!isTRUE(all(n_des == n_kept)))
  warning("The results files have n = ", paste(n_des, collapse = "/"), " but the run's settings say ",
          paste(n_kept, collapse = "/"), " (stale files?): the allowances would be off", call. = FALSE)

#### Cells
fmt_p <- function(x) ifelse(is.na(x), "$\\cdot$", sub("^0\\.", ".", sprintf(paste0("%.", p_digits, "f"), x)))
fmt_d <- function(x) ifelse(is.na(x), "$\\cdot$",
                            ifelse(x == 0, "0", ifelse(1e3 * x < 0.1, "$<0.1$", sprintf("%.1f", 1e3 * x))))
big   <- function(x) formatC(x, big.mark = ",", format = "d")
B_of  <- function(tk, key) if (tk == "proj30") res[[key]]$B30 else res[[key]]$B   #B-hat of a projection test
A_of  <- function(tk) if (tk == "proj30") allowance30 else allowance              #and its rounding allowance
hyb_same <- function(key, i) {                             #does Kudrin's hybrid test coincide with T_H here?
  hy <- res[[key]]$p_hyb[i]
  p  <- res[[key]]$p[["T_H"]][i]
  isFALSE(res[[key]]$restricted[i]) || (!is.na(hy) && !is.na(p) && abs(hy - p) < 1e-12)
}
cell_p <- function(tk, key, i) {
  p <- res[[key]]$p[[tk]][i]
  s <- fmt_p(p)
  certified <- !catalog[tk, "has_B"] || !(key %in% B_rows) || isTRUE(B_of(tk, key)[i] > A_of(tk)[i])
  if (bold && !is.na(p) && p < 0.05 && certified) s <- paste0("\\textbf{", s, "}")
  if (tk == "T_H") {
    hy <- res[[key]]$p_hyb[i]
    s  <- if (hyb_same(key, i)) paste0(s, "$^{\\dagger}$") else if (!is.na(hy)) paste0(s, " [", fmt_p(hy), "]") else s
  }
  s
}
cell_B <- function(tk, key, i) if (key %in% B_rows) fmt_d(B_of(tk, key)[i]) else ""

#### Notes
gl       <- function(lbl) if (startsWith(lbl, "$")) lbl else paste0("\\textit{", lbl, "}")
and_list <- function(x) if (length(x) < 2) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
firstup  <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))
treat_note <- function(key) {
  avg <- if (is.na(n_draws[[key]])) "$p$-values are averages over independent draws" else
    sprintf("$p$-values are averages over %d draws", n_draws[[key]])
  paste0(gl(treatments[[key]]), ": ", switch(key,
                                             none       = "each $t$-statistic is the ratio of the reported coefficient and standard error, and reported $t$-statistics are used as printed.",
                                             uniform    = paste0("the coefficient and the standard error are each replaced by an independent uniform draw from its rounding interval, and reported $t$-statistics receive uniform noise at their last digit; ", avg, "."),
                                             smin       = "each coefficient is divided by the smallest standard error consistent with the reported one, and reported $t$-statistics are used as printed (Section 5.2); this treatment is deterministic.",
                                             smin_noise = paste0("each coefficient, plus uniform noise at its last digit, is divided by the smallest standard error consistent with the reported one, and reported $t$-statistics receive uniform noise at their last digit (Section 5.2); ", avg, ".")))
}
clust_note_for <- function(tests) {
  cl_yes  <- tests[catalog[tests, "clust"] %in% TRUE]
  cl_no   <- tests[catalog[tests, "clust"] %in% FALSE]
  by_boot <- intersect(c("proj", "proj30", "proj_noshift", "T_inf"), tests)
  if (!length(c(cl_yes, cl_no))) return(NULL)
  paste0(
    if (length(cl_yes)) paste0(firstup(and_list(catalog[cl_yes, "short"])), if (length(cl_yes) > 1) " cluster" else " clusters", " by article"),
    if (length(cl_yes) && length(cl_no)) ", and ",
    if (length(cl_no)) paste0(if (length(cl_yes)) and_list(catalog[cl_no, "short"]) else firstup(and_list(catalog[cl_no, "short"])),
                              if (length(cl_no) > 1) " do not" else " does not", if (!length(cl_yes)) " cluster by article"),
    if (length(by_boot)) paste0("; ", and_list(catalog[by_boot, "short"]), if (length(by_boot) > 1) " use " else " uses ",
                                big(boots), " article-level bootstrap repetitions"), ".")
}
bold_note_for <- function(tests) if (bold) paste0("Bold: $p<.05$",
                                                  if (any(catalog[tests, "has_B"]) && length(B_rows))
                                                    paste0(", except that under ", and_list(unname(treatments[B_rows])), " a projection test is bold only when its ",
                                                           "$\\widehat{B}$ exceeds the rounding allowance"), ".")
any_hyb_differs <- any(vapply(treat_keys, function(key) any(!vapply(seq_along(meth), function(i) hyb_same(key, i), TRUE)), TRUE))
ns_draws <- vapply(intersect(averaged, treat_keys), function(key) {   #draws behind the unshifted test's averages
  path <- file.path(res_dir, paste0("MM_weights_results_draws_", key, run_tag))
  if (!file.exists(path)) return(NA_integer_)
  d <- readRDS(path)
  if (is.null(d$shift)) return(NA_integer_)
  as.integer(length(unique(d$draw[d$J0 == 20 & d$shift == 0])))
}, integer(1))
ns_fewer <- !is.na(ns_draws) & ns_draws != n_draws[names(ns_draws)]

notes_main <- c(
  paste0("\\textit{Notes:} $p$-values for the null hypothesis of no selective reporting in the data of ", cite[["brodeur"]],
         ", by research design (RCT: randomized controlled trials; IV: instrumental variables; DID: difference-in-differences; ",
         "RDD: regression discontinuity), with one row for each treatment of rounding."),
  vapply(treat_keys, treat_note, ""),
  paste0("The sample is the same in every row: it excludes tests whose reported coefficient is missing",
         if (drop_zero_mu) paste0(" or zero, following ", cite[["kudrin"]]) else
           paste0(" (unlike ", cite[["kudrin"]], ", it keeps those reported as zero)"),
         ", and tests whose reported standard error is zero or whose rounding unit in $t$-units, $\\kappa$, exceeds ",
         kappa_cap, ". The first line of each block gives the design's numbers of $t$-statistics and articles."),
  sprintf(paste0("%s: $J_0=%s$, $J=%s$, and a grid of %s points on $[-%s,%s]$ plus the point 9999; the discretization error ",
                 "$E(L,\\delta)=%s$ is included in every $p$-value and in $\\widehat{B}$."),
          gl("Projection test"), J0, J, big(numgrid), L, L, format(E_LD, scientific = FALSE)),
  paste0("$10^3\\widehat{B}$: the breakdown statistic (Section 5.1), a lower 95\\% confidence bound on the distance between the ",
         "reported $t$-curve and the set of $t$-curves possible without selective reporting, shown under ",
         and_list(unname(treatments[B_rows])), ". ", gl("Rounding allowance"), " ($\\times 10^3$, first line of each block): ",
         "an upper bound on the distance that rounding can create under these treatments, ", C_round[[as.character(J0)]],
         " times the average of $\\kappa^6$ (Propositions 3 and 4). A rejection with $\\widehat{B}$ above it cannot be attributed to rounding."),
  paste0("$T_\\infty$, $T_H$: the KS-type and histogram-based tests of ", cite[["kudrin"]], ", using the sinc kernel with frequency ",
         "cutoff 3 and, for $T_H$, the points $0.01,0.02,\\ldots,0.15$. His hybrid test $T_H^{\\mathrm{hyb}}$ was run as well; ",
         "$^\\dagger$ marks the cells where it coincides with $T_H$ because its shape restrictions do not bind",
         if (any_hyb_differs) "; elsewhere its $p$-value follows in brackets", "."),
  paste0(gl("EKW"), ": the tests of ", cite[["ekw"]], ", run with their Tests.R on the two-sided $p$-values $2(1-\\Phi(|t|))$: ",
         "the density discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$ (Binom.), and the ", cite[["coxshi"]],
         " histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) and for 2-monotonicity with bounds (CS2B)."),
  clust_note_for(tables_spec$main$tests),
  bold_note_for(tables_spec$main$tests),
  paste0("Appendix Table~\\ref{", tables_spec$appendix$label, "} reports the projection test with $J_0=30$ and without the shift, ",
         "the $p$-curve test of ", cite[["simonsohn"]], ", and the EKW tests LCM and Fisher."))

runtime_note <- if (isTRUE(show_runtime) && file.exists(file.path(res_dir, paste0("Run_times", run_tag)))) {
  rt <- readRDS(file.path(res_dir, paste0("Run_times", run_tag)))
  rr <- function(x) if (diff(range(x)) < 1) sprintf("%.1f", max(x)) else sprintf("%.0f--%.0f", min(x), max(x))
  sprintf("Run times on one core, per design: %s seconds for the projection test with %s bootstrap repetitions, %s seconds for the three tests of %s and %s seconds for the six EKW tests.",
          rr(rt[, "projection"]), big(boots), rr(rt[, "kudrin"]), cite[["kudrin"]], rr(rt[, "ekw"]))
}
notes_appendix <- c(
  paste0("\\textit{Notes:} $p$-values for the null hypothesis of no selective reporting, by research design and treatment of rounding, ",
         "for the sample and the treatments of Table~\\ref{", tables_spec$main$label, "}, which defines them."),
  paste0("Projection, $J_0=30$: the projection test with weights $w_j=\\min\\{1,2^{30-j}\\}$; its $10^3\\widehat{B}$ and the rounding ",
         "allowance in the first line of each block are for $J_0=30$, with ", C_round[["30"]], " in place of ",
         C_round[[as.character(J0)]], "."),
  paste0("Projection, no shift: the projection test on $t$-statistics that are not shifted by 1.96 (Remark 13); no rounding ",
         "allowance is computed for it, so its rejections under $s_{\\min}$ are not certified against rounding",
         if (any(ns_fewer)) paste0("; under ", and_list(unname(treatments[names(ns_draws)[ns_fewer]])), " its $p$-values average over ",
                                   and_list(as.character(ns_draws[ns_fewer])), " draws"), "."),
  paste0("$p$-curve: the test for left skew of ", cite[["simonsohn"]], ", which multiplies each $p<.05$ by 20, transforms the results ",
         "by the normal quantile function and combines them by Stouffer's method; like $p$-curve, it treats the $p$-values as independent."),
  paste0(gl("EKW"), ": the least concave majorant test (LCM) and Fisher's test of ", cite[["ekw"]], " on $p$-values in $(0,0.15]$, ",
         "both of a concave $p$-curve CDF."),
  clust_note_for(tables_spec$appendix$tests),
  bold_note_for(tables_spec$appendix$tests),
  runtime_note)

#### LaTeX: one block of rows per design, one row per treatment, one column per test (two for a projection test)
build_table <- function(spec, notes) {
  tests  <- spec$tests
  fam    <- unlist(lapply(tests, function(tk) rep(catalog[tk, "family"], 1 + catalog[tk, "has_B"])))
  heads  <- unlist(lapply(tests, function(tk) c(catalog[tk, "colhead"], if (catalog[tk, "has_B"]) "$10^3\\widehat{B}$")))
  ncols  <- 1 + length(heads)
  runs   <- rle(fam)
  ends   <- 1 + cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  b_tk   <- tests[catalog[tests, "has_B"]][1]                 #the projection test whose allowance the block line gives
  body <- unlist(lapply(seq_along(meth), function(i) {
    line <- sprintf("\\textit{%s}: %s $t$-statistics in %s articles%s", meth[i], big(n_des[i]), big(a_des[i]),
                    if (!is.na(b_tk)) paste0("; rounding allowance ", fmt_d(A_of(b_tk)[i])) else "")
    rows <- vapply(treat_keys, function(key) paste0(treatments[[key]], " & ", paste(unlist(lapply(tests, function(tk)
      c(cell_p(tk, key, i), if (catalog[tk, "has_B"]) cell_B(tk, key, i)))), collapse = " & "), " \\\\"), "")
    c(if (i > 1) "\\addlinespace", sprintf("\\multicolumn{%d}{l}{%s} \\\\", ncols, line), rows)
  }))
  c("\\begin{table}[htbp]", "\\centering", paste0("\\caption{", spec$caption, "}"), paste0("\\label{", spec$label, "}"),
    "\\small",
    "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
    sprintf("\\begin{tabular}{l%s}", strrep("c", ncols - 1)), "\\toprule",
    paste0(" & ", paste(sprintf("\\multicolumn{%d}{c}{%s}", runs$lengths, runs$values), collapse = " & "), " \\\\"),
    paste(sprintf("\\cmidrule(lr){%d-%d}", starts, ends), collapse = " "),
    paste0("Treatment & ", paste(heads, collapse = " & "), " \\\\"),
    "\\midrule", body, "\\bottomrule", "\\end{tabular}}", "", "\\medskip",
    "\\begin{minipage}{\\linewidth}", "\\footnotesize", notes, "\\end{minipage}", "\\end{table}")
}
notes_of <- list(main = notes_main, appendix = notes_appendix)
for (nm in names(tables_spec)) {
  writeLines(build_table(tables_spec[[nm]], notes_of[[nm]]), file.path(res_dir, tables_spec[[nm]]$file))
  cat("LaTeX table written to", file.path(res_dir, tables_spec[[nm]]$file), "\n")
}

#### Console version and numbers
plain <- function(s) gsub("$^{\\dagger}$", "+", gsub("$\\cdot$", "-", gsub("$<0.1$", "<0.1", gsub("^\\\\textbf\\{([^}]*)\\}", "\\1*", s),
                                                                           fixed = TRUE), fixed = TRUE), fixed = TRUE)
for (nm in names(tables_spec)) {
  tests <- tables_spec[[nm]]$tests
  cat("\n==========", tables_spec[[nm]]$caption, "(* p < .05 and certified; + T_H^hyb = T_H) ==========\n")
  for (i in seq_along(meth)) {
    M <- t(vapply(treat_keys, function(key) unlist(lapply(tests, function(tk)
      c(cell_p(tk, key, i), if (catalog[tk, "has_B"]) cell_B(tk, key, i)))), character(length(tests) + sum(catalog[tests, "has_B"]))))
    colnames(M) <- unlist(lapply(tests, function(tk) c(catalog[tk, "plain"], if (catalog[tk, "has_B"]) paste0(catalog[tk, "plain"], ":Bhat"))))
    rownames(M) <- treat_plain[treat_keys]
    b_tk <- tests[catalog[tests, "has_B"]][1]
    cat(sprintf("\n%s: n = %s, %s articles%s\n", meth[i], big(n_des[i]), big(a_des[i]),
                if (!is.na(b_tk)) paste0(", rounding allowance (x1e3) ", plain(fmt_d(A_of(b_tk)[i]))) else ""))
    print(noquote(plain(M)), right = TRUE)
  }
}
numbers <- do.call(rbind, lapply(shown, function(tk) do.call(rbind, lapply(treat_keys, function(key)
  data.frame(design = meth, test = tk, treatment = key, p = res[[key]]$p[[tk]], stringsAsFactors = FALSE)))))
saveRDS(list(p = numbers,
             breakdown = sapply(treat_keys, function(key) setNames(res[[key]]$B, meth)),
             breakdown_J0_30 = if ("proj30" %in% shown) sapply(treat_keys, function(key) setNames(res[[key]]$B30, meth)),
             rounding_allowance = setNames(allowance, meth), rounding_allowance_J0_30 = setNames(allowance30, meth),
             hybrid_equals_T_H = sapply(treat_keys, function(key) vapply(seq_along(meth), function(i) hyb_same(key, i), TRUE)),
             n = setNames(n_des, meth), articles = setNames(a_des, meth), n_draws = n_draws, kappa_cap = kappa_cap, J0 = J0),
        file.path(res_dir, paste0("Table_rounding_ekw_numbers", run_tag)))

#### 3. Figure: the part of the reported t-curve no honest t-curve explains (R1 17b) --------------------
fr      <- readRDS(file.path(res_dir, paste0("Figdata_residual_tcurve_", S$deround_keys[S$fig_resid_option], run_tag)))
fig_key <- fr$key
u_grid  <- fr$u
fig_res <- fr$results
runs_text <- function(mask) {                                       #"1.42-1.71, 2.03-2.39": stretches of |t| where mask holds
  if (!any(mask)) return("")
  k <- which(mask); brk <- c(0, which(diff(k) > 1), length(k))
  paste(vapply(seq_len(length(brk) - 1), function(i) sprintf("%.2f-%.2f", u_grid[k[brk[i] + 1]], u_grid[k[brk[i + 1]]]), ""),
        collapse = ", ")
}
B_fig <- tryCatch({                                                 #B-hat from the run's results
  d <- readRDS(file.path(res_dir, paste0("MM_weights_results_20_J0_", fig_key, run_tag)))
  setNames(num(d$breakdown), d$methods)
}, error = function(e) NULL)
fig_summary <- do.call(rbind, lapply(meth, function(m) {
  z <- fig_res[[m]]
  data.frame(design = m, distance_x1000 = round(1e3 * z$d, 2),
             peak_excess_at = u_grid[which.max(z$R)], peak_deficit_at = u_grid[which.min(z$R)],
             excess_where  = if (is.null(z$band)) NA else runs_text(z$band[1, ] > 0),
             deficit_where = if (is.null(z$band)) NA else runs_text(z$band[2, ] < 0), stringsAsFactors = FALSE)
}))
cat("\nWhere each design's distortion lies (", S$deround_labels[S$fig_resid_option], "); excess_where and deficit_where: ",
    "stretches of |t| whose 90% band excludes zero:\n", sep = "")
print(fig_summary, row.names = FALSE)
fig_file <- file.path(res_dir, paste0("Figure_residual_tcurve_", fig_key, run_tag))
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
cat("Figure written to", paste0(fig_file, ".pdf/.png"), "\n")

#### 4. Figure: reported |t|, the closest honest t-curve, and where the shape rules out honest reporting --
fs     <- readRDS(file.path(res_dir, paste0("Figdata_tcurve_shape_", S$deround_keys[S$fig_shape_option], run_tag)))
key    <- fs$key
res    <- fs$results
breaks <- fs$breaks
u_c    <- fs$u_c
u_h    <- fs$u_h
max_t  <- max(breaks)
level  <- fs$level
window <- fs$window
show_pointwise <- fig_show_pointwise
show_honest    <- fig_show_honest
as_text <- function(st) if (!nrow(st)) "none" else paste(sprintf("%.2f-%.2f", st[, 1], st[, 2]), collapse = ", ")
tests <- tryCatch({                                                 #the projection test's p-value and B-hat under this option
  r <- readRDS(file.path(res_dir, paste0("MM_weights_results_20_J0_", key, run_tag)))
  data.frame(p = num(r$pval), B = num(r$breakdown), row.names = as.character(r$methods))
}, error = function(e) NULL)
if (is.null(tests)) message("No results file MM_weights_results_20_J0_", key, run_tag, ": the panels go without p-values")
dark  <- grDevices::adjustcolor("#D55E00", alpha.f = 0.45)
light <- grDevices::adjustcolor("#E69F00", alpha.f = 0.25)
file  <- file.path(res_dir, paste0("Figure_tcurve_shape_", key, run_tag))
ymax  <- 1.08 * max(unlist(lapply(res, function(z) c(z$density, if (show_honest) z$honest))))
any_sim <- any(vapply(res, function(z) nrow(z$sim) > 0, TRUE))                    #legend: only what some panel shows
any_pw  <- show_pointwise && any(vapply(res, function(z) nrow(z$pw) > 0, TRUE))
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
      keep <- c(TRUE, any_sim, any_pw, show_honest)
      legend("topright", bty = "n", cex = 0.8,
             legend = c("Reported |t|", sprintf("Too peaked for honest reporting (%g%%, %g-%g)", 100 * level, window[1], window[2]),
                        sprintf("Too peaked for honest reporting (%g%% pointwise)", 100 * level), "Closest honest t-curve")[keep],
             fill = c("grey75", dark, light, NA)[keep], border = c("white", NA, NA, NA)[keep],
             lty = c(NA, NA, NA, 1)[keep], lwd = c(NA, NA, NA, 1.8)[keep], col = c(NA, NA, NA, "black")[keep])
    }
  }
  par(op)
  grDevices::dev.off()
}
shape_summary <- data.frame(design = meth, n = vapply(res, `[[`, 0, "n"),
                            too_peaked_simultaneous = vapply(res, function(z) as_text(z$sim), ""),
                            too_peaked_pointwise = vapply(res, function(z) as_text(z$pw), ""),
                            min_curvature = round(vapply(res, function(z) min(z$D), 0), 3),
                            at_t = vapply(res, function(z) u_c[which.min(z$D)], 0), row.names = NULL)
cat("\nStretches of |t| where the reported t-curve is too peaked for honest reporting (", S$deround_labels[S$fig_shape_option],
    ", h = ", fs$bandwidth, "; simultaneous over ", window[1], "-", window[2], "):\n", sep = "")
print(shape_summary, row.names = FALSE)
write.csv(shape_summary, paste0(file, "_check.csv"), row.names = FALSE)
cat("Figure written to ", file, ".pdf/.png\n", sep = "")