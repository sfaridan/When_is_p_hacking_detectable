#### Table 3 and its appendix companion: tests for selective reporting under four treatments of rounding
# Main table (Table 3): the projection test (J0 = 20) with its B-hat, Kudrin's T_inf and T_H, and the EKW tests
# Disc., Binomial, CS1 and CS2B. Kudrin's hybrid test is run as well: a dagger on T_H marks every cell where
# the two coincide, and where they do not, the hybrid's p-value follows in brackets. Appendix table: the
# projection test with J0 = 30 (with its own B-hat and allowance) and without the 1.96 shift, Stouffer's p-curve
# test, and the EKW tests LCM and Fisher. Both tables have one block of rows per design, whose first line gives
# its sample and rounding allowance, and one row per treatment of rounding, so each column shows how a test's
# verdict moves with the treatment.
# Runs no tests: reads what application_new_deround_all.R saves in <root>/results for every treatment shown
# (1 "none", 2 "uniform", 3 "smin", 4 "smin_noise"). The rounding allowance is recomputed from the printed values
# with that script's rules.

root          <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
data_location <- "C:/Users/sfaridani6/Documents/Data/120246-V1"
res_dir       <- file.path(root, "results")

#--- Kudrin's (2025) drop of coefficients printed as zero: must match application_MM_all_four_rounding_withstouffer.R.
#    When that script sources this one, its value is used; run on its own, set it here. FALSE reads and writes the
#    files tagged "_keepmu0", so the tables of the run without the drop never overwrite the others.
if (!exists("drop_zero_mu")) drop_zero_mu <- TRUE
run_tag <- if (drop_zero_mu) "" else "_keepmu0"
cat("Tables for the sample that", if (drop_zero_mu) "drops" else "keeps", "coefficients reported as zero",
    if (nzchar(run_tag)) paste0("(files tagged ", run_tag, ")") else "(untagged files)", "\n")

#--- Must match the settings application_new_deround_all.R ran with
kappa_cap     <- 0.5      #its kappa_cap
se_rounded_up <- FALSE
J0            <- 20
J             <- 40
L             <- 6.5
numgrid       <- 3000
boots         <- 1000
E_LD          <- 0.00063  #compute_epsilons(L, nx = numgrid)$epsilon_U
C_round       <- c("20" = 0.13, "30" = 0.31)   #C of Propositions 3 and 4, by J0

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
bold         <- TRUE      #bold p < .05; under B_rows a projection test only when its B-hat exceeds the allowance
p_digits     <- 2
show_runtime <- TRUE      #run-time sentence in the appendix notes, when application_new_deround_all.R saved run times
cite <- c(brodeur = "Brodeur et al.\\ (2020)", ekw = "Elliott et al.\\ (2022)", kudrin = "Kudrin (2025)",
          coxshi  = "Cox and Shi (2023)", simonsohn = "Simonsohn et al.\\ (2014)")   #or "\\citet{key}"

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
  if (!file.exists(f)) stop("No results file ", basename(f), ": run application_new_deround_all.R with option ",
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

#### Rounding allowance (rules copied from application_new_deround_all.R)
n_dec <- function(x, max_dec = 12) {
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
rounding_units <- function(d) {
  dec_mu <- n_dec(d$mu); dec_sd <- n_dec(d$sd); dec_t <- n_dec(d$t)
  mu_p   <- round(d$mu, dec_mu); sd_p <- round(d$sd, dec_sd); t_p <- round(d$t, dec_t)
  tab_id <- paste(d$title, d$table, sep = " | ")
  b_sd   <- 10^(-dec_sd)
  a_tab  <- 10^(-ave(dec_mu, tab_id, FUN = function(x) if (all(is.na(x))) NA_integer_ else max(x, na.rm = TRUE)))
  a_mu   <- pmax(a_tab, b_sd)
  s_min  <- sd_p - if (se_rounded_up) b_sd else b_sd / 2
  a_t    <- 10^(-dec_t)
  has_ratio  <- !is.na(mu_p) & !is.na(sd_p)
  drop_se0   <- has_ratio & sd_p <= 0
  drop_kappa <- has_ratio & !drop_se0 & a_mu / s_min > kappa_cap
  use_ratio  <- has_ratio & !drop_se0 & !drop_kappa
  direct     <- !has_ratio & !is.na(t_p)
  kappa      <- ifelse(use_ratio, a_mu / s_min, ifelse(direct, a_t, NA))
  keep       <- !(drop_se0 | drop_kappa)
  data.frame(method = as.character(d$method)[keep], kappa = kappa[keep], stringsAsFactors = FALSE)
}
MM_raw      <- haven::read_dta(file.path(data_location, "MM data.dta"))
MM          <- MM_raw[!(is.na(MM_raw$mu) | (drop_zero_mu & MM_raw$mu %in% 0)), ]
kap         <- rounding_units(MM)
by_design   <- factor(kap$method, levels = meth)
n_kept      <- as.vector(table(by_design))
kappa6      <- as.vector(tapply(kap$kappa^6, by_design, mean, na.rm = TRUE))
allowance   <- C_round[[as.character(J0)]] * kappa6
allowance30 <- C_round[["30"]] * kappa6                    #for the J0 = 30 test
n_des <- res[[1]]$n_proj                                   #t-statistics and articles per design, as tested
a_des <- res[[1]]$articles
for (key in treat_keys) {
  if (!isTRUE(all(res[[key]]$articles == a_des)) || !isTRUE(all(res[[key]]$articles_kud == a_des)))
    warning("Article counts in the '", key, "' results differ from the '", treat_keys[1], "' results", call. = FALSE)
}
if (!isTRUE(all(n_des == n_kept)))
  warning("The results have n = ", paste(n_des, collapse = "/"), " but kappa_cap = ", kappa_cap, " keeps ",
          paste(n_kept, collapse = "/"), ": match application_new_deround_all.R, or the allowances are off", call. = FALSE)

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