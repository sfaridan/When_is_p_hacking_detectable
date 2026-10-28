#### Table 4 with the EKW tests added
# Same layout as Table 4: one row per design; for each test, one p-value column per treatment of rounding; then
# B-hat and the rounding allowance (x 10^3). The six EKW tests are two more blocks of rows with the same columns
# ("EKW LCM", ... as Table 4 has "Kudrin T_inf"), keeping Table 4's width and font sizes.
# Runs no tests: reads what application_new_deround.R saves in <root>/results, one run per treatment shown
# (deround_option = 1 "none", 2 "uniform", 3 "smin"). The rounding allowance is recomputed from the printed values
# with that script's rules; this script stops unless the implied sample matches n in every results file.

root          <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
data_location <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data"
res_dir       <- file.path(root, "results")
out_file      <- "table_rounding_ekw.tex"

#--- Must match the settings application_new_deround.R ran with
kappa_cap     <- 1        #its current value; Section 5.2.1 and Table 4's notes use 0.3
se_rounded_up <- FALSE
J0            <- 20
J             <- 40
L             <- 6.5
numgrid       <- 3000
boots         <- 1000
E_LD          <- 0.00063  #compute_epsilons(L, nx = numgrid)$epsilon_U
C_round       <- c("20" = 0.13, "30" = 0.31)   #C of Propositions 3 and 4, by J0

#--- Layout
treatments  <- c(none = "None", uniform = "Unif.", smin = "$s_{\\min}$")  #add smin_noise = "$s_{\\min}$+noise"
B_treatment <- "smin"     #the rounding allowance covers "smin" and "smin_noise" only
cite <- c(brodeur = "Brodeur et al.\\ (2020)", ekw = "Elliott et al.\\ (2022)", kudrin = "Kudrin (2025)",
          coxshi  = "Cox and Shi (2023)")                       #or "\\citet{key}"
panels <- list(           #blocks of up to three tests; a non-NULL title adds a "Panel A." row
  list(title = NULL, tests = c("proj", "T_inf", "T_H")),
  list(title = NULL, tests = c("LCM", "Fisher", "Disc")),
  list(title = NULL, tests = c("Binomial", "CS1", "CS2B")))
#Also available: "T_H_hyb", "EKW_min"
bold        <- FALSE
p_digits    <- 2
caption     <- sprintf("Tests for selective reporting under %s treatments of rounding",
                       c("one", "two", "three", "four")[length(treatments)])
label       <- "tab:rounding_ekw"
extra_notes <- if (!"smin_noise" %in% names(treatments))       #set the label of your appendix table
  "A fourth treatment that adds uniform noise to $s_{\\min}$ is reported in Appendix Table~\\ref{tab:application}."

#### Setup
meth      <- c("RCT", "IV", "DID", "RDD")
test_keys <- c("proj", "T_inf", "T_H", "T_H_hyb", "LCM", "Fisher", "Disc", "Binomial", "CS1", "CS2B", "EKW_min")
catalog   <- data.frame(
  head  = c("Projection test", "Kudrin $T_\\infty$", "Kudrin $T_H$", "Kudrin $T_H^{\\mathrm{hyb}}$",
            "EKW LCM", "EKW Fisher", "EKW Disc.", "EKW Binomial", "EKW CS1", "EKW CS2B", "EKW Min"),
  short = c("the projection test", "$T_\\infty$", "$T_H$", "$T_H^{\\mathrm{hyb}}$",
            "LCM", "Fisher", "Disc.\\@", "Binomial", "CS1", "CS2B", "EKW Min"),
  plain = c("Proj", "Tinf", "TH", "THhyb", "LCM", "Fisher", "Disc", "Binom", "CS1", "CS2B", "EKWmin"),
  src   = c("proj", "kud", "kud", "kud", rep("ewk", 7)),
  field = c("pval", "T_inf", "T_H", "T_H_hybrid", "lcms_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK",
            "CS1_EWK_s", "CS2_EWK_s", "min"),
  clust = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, NA),
  row.names = test_keys, stringsAsFactors = FALSE)
opt_of      <- c(none = 1, uniform = 2, smin = 3, smin_noise = 4)
averaged    <- c("uniform", "smin_noise")
treat_plain <- c(none = "None", uniform = "Unif", smin = "smin", smin_noise = "smin+U")
treat_keys  <- names(treatments)
shown       <- unlist(lapply(panels, `[[`, "tests"))
show_dist   <- "proj" %in% shown
keys_read   <- unique(c(treat_keys, if (show_dist) B_treatment))
stopifnot(all(treat_keys %in% names(opt_of)), all(shown %in% test_keys), !anyDuplicated(shown),
          B_treatment %in% c("smin", "smin_noise"), as.character(J0) %in% names(C_round))
treat_lab <- c(none = "None", uniform = "Unif.", smin = "$s_{\\min}$", smin_noise = "$s_{\\min}$+noise")
treat_lab[treat_keys] <- treatments

num  <- function(x) vapply(x, function(v) as.numeric(v)[1], numeric(1), USE.NAMES = FALSE)
pick <- function(d, what) {
  i <- match(meth, d$methods)
  if (anyNA(i)) stop(what, " has no results for ", paste(meth[is.na(i)], collapse = ", "), call. = FALSE)
  d[i, ]
}
res_file <- function(stem, key) {
  f <- file.path(res_dir, paste0(stem, "_", key))
  if (!file.exists(f)) stop("No results file ", basename(f), ": run application_new_deround.R with deround_option = ",
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
  if (!is.null(proj$numcoeffs) && any(num(proj$numcoeffs) != J + 1))
    warning(proj_name, "_", key, " was run with numcoeffs != ", J + 1, call. = FALSE)
  p <- lapply(setNames(shown, shown), function(k) {
    f <- catalog[k, "field"]
    switch(catalog[k, "src"],
           proj = num(proj$pval),
           kud  = if (f %in% names(kud)) num(kud[[f]]) else stop("Kudrin_results_", key, " has no ", f, call. = FALSE),
           ewk  = if (f %in% rownames(ewk)) unname(ewk[f, ]) else stop("EWK_results_", key, " has no ", f, call. = FALSE))
  })
  list(p = p, B = num(proj$breakdown), n_proj = num(proj$n), n_kud = num(kud$n),
       restricted = if (is.null(kud$restricted)) rep(NA, length(meth)) else as.logical(num(kud$restricted)))
}
res <- setNames(lapply(keys_read, read_results), keys_read)

count_draws <- function(key) {
  from_file <- function(stem, f) {
    path <- file.path(res_dir, paste0(stem, "_", key))
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

#### Rounding allowance (rules copied from application_new_deround.R)
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
MM_raw    <- haven::read_dta(file.path(data_location, "MM data.dta"))
MM        <- MM_raw[!(is.na(MM_raw$mu) | MM_raw$mu == 0), ]
kap       <- rounding_units(MM)
by_design <- factor(kap$method, levels = meth)
n_kept    <- as.vector(table(by_design))
allowance <- C_round[[as.character(J0)]] * as.vector(tapply(kap$kappa^6, by_design, mean, na.rm = TRUE))
check_n <- function(n, what) if (!isTRUE(all(n == n_kept)))
  stop(what, " has n = ", paste(n, collapse = "/"), " but kappa_cap = ", kappa_cap, " and se_rounded_up = ",
       se_rounded_up, " keep ", paste(n_kept, collapse = "/"), ". Match application_new_deround.R or re-run it.",
       call. = FALSE)
for (key in keys_read) {
  check_n(res[[key]]$n_proj, paste0("MM_weights_results_", J0, "_J0_", key))
  check_n(res[[key]]$n_kud,  paste0("Kudrin_results_", key))
}

#### Cells
fmt_p <- function(x) ifelse(is.na(x), "$\\cdot$", sub("^0\\.", ".", sprintf(paste0("%.", p_digits, "f"), x)))
fmt_d <- function(x) ifelse(is.na(x), "$\\cdot$",
                            ifelse(x == 0, "0", ifelse(1e3 * x < 0.1, "$<0.1$", sprintf("%.1f", 1e3 * x))))
cells <- lapply(setNames(shown, shown), function(tk) lapply(setNames(treat_keys, treat_keys), function(key) {
  p <- res[[key]]$p[[tk]]
  b <- bold & !is.na(p) & p < 0.05
  if (tk == "proj" && key %in% c("smin", "smin_noise"))
    b <- b & !is.na(res[[key]]$B) & !is.na(allowance) & res[[key]]$B > allowance
  ifelse(b, paste0("\\textbf{", fmt_p(p), "}"), fmt_p(p))
}))
dist_cells <- if (show_dist) cbind(fmt_d(res[[B_treatment]]$B), fmt_d(allowance))

#### Notes
gl       <- function(lbl) if (startsWith(lbl, "$")) lbl else paste0("\\textit{", lbl, "}")
and_list <- function(x) if (length(x) < 2) x else paste(paste(x[-length(x)], collapse = ", "), "and", x[length(x)])
big      <- function(x) formatC(x, big.mark = ",", format = "d")
firstup  <- function(s) paste0(toupper(substring(s, 1, 1)), substring(s, 2))
nT <- length(treat_keys); G <- max(lengths(lapply(panels, `[[`, "tests"))); P <- G * nT
titled <- !vapply(panels, function(pn) is.null(pn$title), logical(1))
intro <- c(
  paste0("\\textit{Notes:} ",
         if (length(panels) > 1 && all(titled)) sprintf("Panels A--%s report", LETTERS[length(panels)]) else
           if (length(panels) > 1) "The table reports" else paste("The first", P, "columns report"),
         " $p$-values for the null hypothesis of no selective reporting in the data of ", cite[["brodeur"]],
         ", by research design (RCT: randomized controlled trials; IV: instrumental variables; ",
         "DID: difference-in-differences; RDD: regression discontinuity)."),
  "Columns differ in the treatment of rounding.")
treat_note <- function(key) {
  avg <- if (is.na(n_draws[[key]])) "$p$-values are averages over independent draws" else
    sprintf("$p$-values are averages over %d draws", n_draws[[key]])
  paste0(gl(treat_lab[[key]]), ": ", switch(key,
                                            none       = "each $t$-statistic is the ratio of the reported coefficient and standard error, and reported $t$-statistics are used as printed.",
                                            uniform    = paste0("the coefficient and the standard error are each replaced by an independent uniform draw from its rounding interval, and reported $t$-statistics receive uniform noise at their last digit; ", avg, "."),
                                            smin       = "each coefficient is divided by the smallest standard error consistent with the reported one, and reported $t$-statistics are used as printed (Section 5.2); this treatment is deterministic.",
                                            smin_noise = paste0("each coefficient, plus uniform noise at its last digit, is divided by the smallest standard error consistent with the reported one, and reported $t$-statistics receive uniform noise at their last digit (Section 5.2); ", avg, ".")))
}
sample_note <- c("The sample is the same in every column.",
                 paste0("It excludes tests whose reported coefficient is missing or zero, following ", cite[["kudrin"]],
                        ", and tests whose reported standard error is zero or whose rounding unit in $t$-units, $\\kappa$, exceeds ",
                        kappa_cap, "."))
proj_note <- if ("proj" %in% shown) sprintf(paste0(
  "%s: $J_0=%s$, $J=%s$, and a grid of %s points on $[-%s,%s]$ plus the point 9999; ",
  "the discretization error $E(L,\\delta)=%s$ is included in every $p$-value and in $\\widehat{B}$."),
  gl("Projection test"), J0, J, big(numgrid), L, L, format(E_LD, scientific = FALSE))
kud_shown  <- intersect(c("T_inf", "T_H", "T_H_hyb"), shown)
kud_hist   <- intersect(c("T_H", "T_H_hyb"), shown)
restricted <- unlist(lapply(res[treat_keys], `[[`, "restricted"))
hyb_same   <- isTRUE(!any(restricted))
if (!"T_H_hyb" %in% shown && any(restricted, na.rm = TRUE))
  warning("T_H_hybrid differs from T_H in some cells: consider showing \"T_H_hyb\"", call. = FALSE)
kud_note <- if (length(kud_shown)) paste0(
  paste(catalog[kud_shown, "short"], collapse = ", "), ": the ",
  and_list(c(T_inf = "KS-type", T_H = "histogram-based", T_H_hyb = "hybrid histogram-based")[kud_shown]),
  if (length(kud_shown) > 1) " tests" else " test", " of ", cite[["kudrin"]], ", using the sinc kernel with frequency cutoff 3",
  if (length(kud_hist)) paste0(" and, for ", and_list(catalog[kud_hist, "short"]), ", the points $0.01,0.02,\\ldots,0.15$"),
  if (!"T_H_hyb" %in% shown && hyb_same) "; his hybrid test coincides with $T_H$ in every case", ".")
ekw_note <- if (length(intersect(c("LCM", "Fisher", "Disc", "Binomial", "CS1", "CS2B"), shown))) paste0(
  gl("EKW"), ": the tests of ", cite[["ekw"]], ", run with their Tests.R on the two-sided $p$-values ",
  "$2(1-\\Phi(|t|))$, namely the least concave majorant test (LCM) and Fisher's test on $(0,0.15]$, the density ",
  "discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$, and the ", cite[["coxshi"]],
  " histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) and for 2-monotonicity with bounds (CS2B).")
min_note <- if ("EKW_min" %in% shown) paste0(gl("EKW Min"), ": the smallest of the six EKW $p$-values",
                                             if (any(treat_keys %in% averaged))
                                               paste0(" (of the averaged $p$-values under ", and_list(treat_lab[intersect(averaged, treat_keys)]), ")"), ".")
cl_yes <- shown[catalog[shown, "clust"] %in% TRUE]
cl_no  <- shown[catalog[shown, "clust"] %in% FALSE]
by_boot <- intersect(c("proj", "T_inf"), shown)
clust_note <- if (length(c(cl_yes, cl_no))) paste0(
  if (length(cl_yes)) paste0(firstup(and_list(catalog[cl_yes, "short"])),
                             if (length(cl_yes) > 1) " cluster" else " clusters", " by article"),
  if (length(cl_yes) && length(cl_no)) ", and ",
  if (length(cl_no)) paste0(if (length(cl_yes)) and_list(catalog[cl_no, "short"]) else firstup(and_list(catalog[cl_no, "short"])),
                            if (length(cl_no) > 1) " do not" else " does not", if (!length(cl_yes)) " cluster by article"),
  if (length(by_boot)) paste0("; ", and_list(catalog[by_boot, "short"]), if (length(by_boot) > 1) " use " else " uses ",
                              big(boots), " article-level bootstrap repetitions"), ".")
dist_note <- if (show_dist) c(
  paste0("$\\widehat{B}$: the breakdown statistic under ", treat_lab[[B_treatment]], " (Section 5.1), a lower 95\\% ",
         "confidence bound on the distance between the reported $t$-curve and the set of $t$-curves possible without selective reporting."),
  paste0(gl("Rounding"), ": an upper bound on the distance that rounding can create under ", treat_lab[[B_treatment]],
         ", equal to ", C_round[[as.character(J0)]], " times the average of $\\kappa^6$ (Propositions 3 and 4)."),
  "A rejection with $\\widehat{B}$ above this bound cannot be attributed to rounding.")
bold_note <- if (bold) paste0("Bold: $p<.05$",
                              if ("proj" %in% shown && any(treat_keys %in% c("smin", "smin_noise")))
                                paste0(", except that the projection test under ", and_list(treat_lab[intersect(c("smin", "smin_noise"), treat_keys)]),
                                       " is bold only when its $\\widehat{B}$ exceeds the rounding allowance"), ".")
notes <- c(intro, vapply(treat_keys, treat_note, ""), sample_note, proj_note, kud_note, ekw_note, min_note,
           clust_note, dist_note, bold_note, extra_notes)

#### LaTeX
ncols <- 1 + P + 2 * show_dist
mc    <- function(k, text, align = "c") sprintf("\\multicolumn{%d}{%s}{%s}", as.integer(k), align, text)
cmid  <- function(from, to) paste(sprintf("\\cmidrule(lr){%d-%d}", as.integer(from), as.integer(to)), collapse = " ")
tex_row <- function(x) {
  span  <- rep(1L, length(x))
  is_mc <- startsWith(x, "\\multicolumn{")
  span[is_mc] <- as.integer(sub("^\\\\multicolumn\\{([0-9]+)\\}.*$", "\\1", x[is_mc]))
  stopifnot(sum(span) == ncols)
  paste0(paste(x, collapse = " & "), " \\\\")
}
panel_tex <- function(k) {
  tests <- panels[[k]]$tests; g <- length(tests)
  pad   <- rep("", nT * (G - g))
  has_B <- show_dist && "proj" %in% tests
  rows  <- vapply(seq_along(meth), function(i) tex_row(c(
    meth[i],
    unlist(lapply(tests, function(tk) vapply(treat_keys, function(key) cells[[tk]][[key]][i], ""))),
    pad,
    if (show_dist) (if (has_B) dist_cells[i, ] else c("", "")))), "")
  c(if (k > 1) "\\midrule",
    if (titled[k]) tex_row(mc(ncols, sprintf("\\textit{Panel %s. %s}", LETTERS[k], panels[[k]]$title), "l")),
    tex_row(c("", mc(nT, catalog[tests, "head"]), pad, if (show_dist) c("", ""))),
    cmid(2 + (seq_len(g) - 1) * nT, 1 + seq_len(g) * nT),
    tex_row(c("Design", rep(unname(treatments), g), pad,
              if (show_dist) (if (has_B) c("$\\widehat{B}$", "Rounding") else c("", "")))),
    "\\midrule",
    rows)
}
head_top <- c(tex_row(c("", mc(P, "$p$-value"), if (show_dist) mc(2, "Distance ($\\times 10^3$)"))),
              paste(c(cmid(2, P + 1), if (show_dist) cmid(P + 2, P + 3)), collapse = " "))
tex <- c("\\begin{table}[htbp]", "\\centering",
         paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
         "\\small",
         "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
         paste0("\\begin{tabular}{l", strrep("c", ncols - 1), "}"),
         "\\toprule", head_top, unlist(lapply(seq_along(panels), panel_tex)), "\\bottomrule",
         "\\end{tabular}}", "", "\\medskip",
         "\\begin{minipage}{\\linewidth}", "\\footnotesize", notes, "\\end{minipage}", "\\end{table}")
writeLines(tex, file.path(res_dir, out_file))

#### Console version and numbers
plain <- function(s) gsub("$\\cdot$", "-", gsub("$<0.1$", "<0.1", sub("^\\\\textbf\\{(.*)\\}$", "\\1*", s),
                                                fixed = TRUE), fixed = TRUE)
for (k in seq_along(panels)) {
  tests <- panels[[k]]$tests
  M <- do.call(cbind, lapply(tests, function(tk) sapply(treat_keys, function(key) cells[[tk]][[key]])))
  colnames(M) <- paste(rep(catalog[tests, "plain"], each = nT), treat_plain[treat_keys], sep = ":")
  if (show_dist && "proj" %in% tests)
    M <- cbind(M, "Bhat(x1e3)" = dist_cells[, 1], "Rounding(x1e3)" = dist_cells[, 2])
  rownames(M) <- meth
  cat("\n", paste(catalog[tests, "plain"], collapse = ", "), "\n", sep = "")
  print(noquote(plain(M)), right = TRUE)
}
cat(sprintf("\nSample check passed: n = %s (%s), kappa_cap = %s.\n",
            paste(n_kept, collapse = "/"), paste(meth, collapse = "/"), kappa_cap))
cat("LaTeX table written to", file.path(res_dir, out_file), "\n")
numbers <- do.call(rbind, lapply(shown, function(tk) do.call(rbind, lapply(treat_keys, function(key)
  data.frame(design = meth, test = tk, treatment = key, p = res[[key]]$p[[tk]], stringsAsFactors = FALSE)))))
saveRDS(list(p = numbers, breakdown = if (show_dist) setNames(res[[B_treatment]]$B, meth),
             rounding_allowance = setNames(allowance, meth), n = setNames(n_kept, meth),
             n_draws = n_draws, kappa_cap = kappa_cap, J0 = J0),
        file.path(res_dir, "Table_rounding_ekw_numbers"))