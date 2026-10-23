################################################################################
#  make_sim_table.R
#
#  Turns the summary file saved by run_sims() (functions_all_weighted.R) into a
#  table like Table 1 of the paper, with three more columns for the Kudrin
#  tests T_inf, T_H and T_H_hybrid (kudrin_tests.R).
#
#  The summary file is the parms data frame: one row per parameterization and
#  one rejection-rate column per test. The script prints the table and writes it
#  as a LaTeX table, with a note built from the saved simulation settings, next
#  to the results file. Base R only. Set root below and source the script after
#  run_sims_weighted.R has finished.
################################################################################


# ---------------------------------------------------------------------------- #
#  Settings
# ---------------------------------------------------------------------------- #

root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
simulation_results <- paste0(root, "/simulations")

#First combine EKW and Kudrin results from the J0=30 run with the J0=20 run
J0_30 <- readRDS(paste0(simulation_results,"/to keep/to keep revision/sims_weighted_500014_nsims500_2026-10-23_01-01-05.RDS"))
J0_20 <- readRDS(paste0(simulation_results,"/to keep/to keep revision/sims_weighted_500014_nsims500_2026-10-23_09-59-45.RDS"))
J0_20[,c("rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK", "CS1_EWK",
            "CS2B_EWK", "Kudrin_T_inf", "Kudrin_T_H", "Kudrin_T_H_hybrid")] <- J0_30[,c("rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK", "CS1_EWK",
                                                                                        "CS2B_EWK", "Kudrin_T_inf", "Kudrin_T_H", "Kudrin_T_H_hybrid")] 
saveRDS(J0_20,file=paste0(simulation_results,"/sims_weighted_5000_J0_20_combined.RDS"))

#make the tablea
results_prefix <- "sims_weighted_5000"   # the sim_file_prefix given to run_sims()
results_file   <- paste0(simulation_results,"/sims_weighted_5000_J0_20_combined.RDS")                   # NULL = newest file in simulation_results starting with results_prefix
table_caption  <- NULL                   # NULL = named after the p-hacking, e.g. "Threshold $p$-hacking"
shrink_to_fit  <- TRUE                   # scale the table down to the text width if wider (needs graphicx)
digits         <- 2
na_string      <- "--"                   # shown for tests that were skipped or never returned a p-value

ekw_cite    <- "Elliott et al. (2022)"   # plain text; e.g. "\\citet{key}" for a linked citation
kudrin_cite <- "Kudrin (2024)"           # year as in the paper's references (kudrin_tests.R says 2025)

# Test columns, left to right: column in the results file, LaTeX header, console header, and the
# header spanning its block of columns ("" for none). Blocks are separated by vertical rules.
tests <- data.frame(
  col   = c("reject_rate", "rej_lcm_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK", "CS1_EWK",
            "CS2B_EWK", "Kudrin_T_inf", "Kudrin_T_H", "Kudrin_T_H_hybrid"),
  tex   = c("Proj.", "LCM", "Fisher", "Disc.", "Binom.", "CS1", "CS2B",
            "$T_{\\infty}$", "$T_{H}$", "$T_{H,\\mathrm{hyb}}$"),
  text  = c("Proj.", "LCM", "Fisher", "Disc.", "Binom.", "CS1", "CS2B", "T_inf", "T_H", "T_H_hyb"),
  group = c("", rep(ekw_cite, 6), rep(kudrin_cite, 3)),
  stringsAsFactors = FALSE
)

# Pi_0 column: labels for pi0_shape. Shapes not listed appear under their own names.
pi0_labels <- c(
  "poisson(2)"      = "Poisson(2)",
  "chi2(2)"         = "$\\chi^2(2)$",
  "null"            = "$H = 0$",
  "onehalf"         = "$H = 0.5$",
  "double_normal"   = "Mix Normals",        # 50-50 mixture of N(h_center - 2, sigma_h^2) and N(h_center, sigma_h^2)
  "fitted"          = "Fitted (z-curve)",   # discrete on {0, ..., 6}, fitted by z-curve
  "unif_half"       = "Unif(-0.5,0.5)",
  "standard_normal" = "N(0,1)",
  "hard"            = "$N(2, 0.7^2)$"
)

# p-Hacking column (rows with prob_hack > 0), and the names used in the caption and note
hack_labels <- c(threshold = "Threshold", max = "Max of 8 $t$-stats")
hack_names  <- c(threshold = "threshold $p$-hacking", max = "maximization $p$-hacking")


# ---------------------------------------------------------------------------- #
#  Helpers
# ---------------------------------------------------------------------------- #

escape_tex  <- function(x) gsub("([_%&#$])", "\\\\\\1", x)                         # raw names into LaTeX
tex_to_text <- function(x) gsub("\\$|\\\\(?=[A-Za-z_%&#])", "", x, perl = TRUE)    # LaTeX labels for the console
fmt_rate    <- function(x) ifelse(is.na(x), na_string, formatC(as.numeric(x), format = "f", digits = digits))
fmt_value   <- function(x) if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE) else escape_tex(as.character(x))
uniq        <- function(x) paste(fmt_value(unique(x)), collapse = ", ")            # settings quoted in the note
tex_row     <- function(cells) paste0(paste(cells, collapse = " & "), " \\\\")


# ---------------------------------------------------------------------------- #
#  Read the results
# ---------------------------------------------------------------------------- #

if (is.null(results_file)) {
  files <- list.files(simulation_results, pattern = "\\.rds$", full.names = TRUE)
  files <- files[startsWith(basename(files), results_prefix)]
  if (length(files) == 0) stop("No .rds file starting with '", results_prefix, "' in ", simulation_results)
  # run_sims() ends each file name with a %Y-%m-%d_%H-%M-%S time stamp, which sorts chronologically
  stamp <- sub("^.*_([0-9]{4}-[0-9]{2}-[0-9]{2}_[0-9]{2}-[0-9]{2}-[0-9]{2})\\.rds$", "\\1", basename(files))
  stamp[stamp == basename(files)] <- ""
  results_file <- files[order(stamp, decreasing = TRUE)[1]]
}
message("Reading ", results_file)
res <- readRDS(results_file)

table_file     <- paste0(root,"/tables/", "table_", sub("\\.RDS$", ".tex", basename(results_file))) # NULL = table_<results file name>.tex next to the results file


needed <- c("prob_hack", "hack_type", "theta", "pi0_shape", "finished_at", "n", "nsims", "nboots", "cv",
            "num_coeffs", "numgrid", "L", "J0", "p")
if (length(setdiff(needed, names(res)))) {
  stop("Not a run_sims() results file; it lacks: ", paste(setdiff(needed, names(res)), collapse = ", "))
}
for (v in setdiff(tests$col, names(res))) {
  warning("The results file has no column ", v, "; it is shown as '", na_string, "'.")
  res[[v]] <- NA_real_
}

# run_sims() saves the file after every parameterization and creates each result column with
# parms$col[1] <- x, which R recycles to all rows. Rows that had not run yet when the file was
# saved therefore carry row 1's results, finish time included.



# ---------------------------------------------------------------------------- #
#  Build the table
# ---------------------------------------------------------------------------- #

shape     <- as.character(res$pi0_shape)
type      <- as.character(res$hack_type)
pi0_order <- if (is.factor(res$pi0_shape)) levels(res$pi0_shape) else unique(shape)
hacked    <- unique(type[res$prob_hack > 0])
what      <- if (length(hacked) == 1 && hacked %in% names(hack_names)) hack_names[[hacked]] else NA

# p-Hacking column: the kind of selective reporting in a row
reporting <- function(prob_hack, hack_type, theta, shape) {
  if (shape == "smooth_pub") return("Smooth pub. bias")        # this DGP replaces hacking and truncation
  parts <- character(0)
  if (prob_hack > 0) {
    lab <- if (hack_type %in% names(hack_labels)) hack_labels[[hack_type]] else escape_tex(hack_type)
    if (prob_hack < 1) lab <- sprintf("%s (%g\\%%)", lab, 100 * prob_hack)   # share of researchers who hack
    parts <- c(parts, lab)
  }
  if (theta < 1) parts <- c(parts, sprintf("Pub. bias ($\\theta = %g$)", theta))
  if (length(parts) == 0) "None" else paste(parts, collapse = " + ")
}

tab <- data.frame(
  panel = ifelse(res$prob_hack == 0 & res$theta == 1 & shape != "smooth_pub", "Size", "Power"),
  hack  = mapply(reporting, res$prob_hack, type, res$theta, shape, USE.NAMES = FALSE),
  pi0   = ifelse(shape %in% names(pi0_labels), pi0_labels[shape], escape_tex(shape)),
  stringsAsFactors = FALSE
)

# Settings that vary across rows and are not in the labels (e.g. n) get a column of their own
not_settings <- c("epsilon_U", "nu_resid", "reject_rate_no_adj", "reject_rate_no_nu", "mean_resid",
                  "breakdown", "boot95", "mean_cvs_no_nu", "mean_cvs_full", "finished_at", "time_taken",
                  "script_path", tests$col, grep("^Kudrin_", names(res), value = TRUE),
                  "prob_hack", "hack_type", "theta", "pi0_shape", "omit_proj", "omit_EWK", "kudrin_tests")
settings   <- setdiff(names(res), not_settings)
extra_cols <- settings[vapply(res[settings], function(x) length(unique(x)) > 1, logical(1))]
extra_hdr  <- ifelse(extra_cols == "n", "$n$", ifelse(extra_cols == "nu", "$\\nu$", escape_tex(extra_cols)))
for (v in extra_cols) tab[[v]] <- fmt_value(res[[v]])
for (i in seq_len(nrow(tests))) tab[[tests$text[i]]] <- fmt_rate(res[[tests$col[i]]])

# Size panel first; within a panel, rows by p-hacking, then Pi_0 in the order of the parms grid
ord <- do.call(order, unname(c(list(tab$panel != "Size", match(tab$hack, unique(tab$hack)),
                                    match(shape, pi0_order)), as.list(res[extra_cols]))))
tab <- tab[ord, , drop = FALSE]

dup <- duplicated(tab[c("panel", "hack", "pi0", extra_cols)])
if (any(dup)) {
  warning("Some rows have the same labels, so the table cannot tell them apart: ",
          paste(unique(tex_to_text(paste(tab$hack[dup], tab$pi0[dup], sep = " / "))), collapse = "; "))
}

if (is.null(table_caption)) {
  table_caption <- if (is.na(what)) "Simulated rejection rates" else paste0(toupper(substr(what, 1, 1)), substring(what, 2))
}


# ---------------------------------------------------------------------------- #
#  Print
# ---------------------------------------------------------------------------- #

txt        <- tab
txt[]      <- lapply(txt, tex_to_text)
names(txt) <- c("panel", "p-Hacking", "Pi_0", tex_to_text(extra_hdr), tests$text)
old_opts   <- options(width = max(getOption("width"), 160))
cat("\n", tex_to_text(table_caption), "   (", basename(results_file), ")\n", sep = "")
for (pn in unique(txt$panel)) {
  cat("\n", pn, "\n", sep = "")
  print(txt[txt$panel == pn, -1, drop = FALSE], row.names = FALSE, right = FALSE)
}
options(old_opts)


# ---------------------------------------------------------------------------- #
#  LaTeX
# ---------------------------------------------------------------------------- #

n_lab <- 2 + length(extra_cols)                  # p-Hacking, Pi_0 and any extra settings
n_all <- n_lab + nrow(tests)
runs  <- rle(tests$group)                        # blocks of test columns
spec  <- paste0("ll", strrep("c", length(extra_cols)), paste0("|", strrep("c", runs$lengths), collapse = ""))
group_cells <- unlist(Map(function(g, k, last) {
  if (g == "") rep("", k) else sprintf("\\multicolumn{%d}{c%s}{\\textbf{%s}}", k, if (last) "" else "|", g)
}, runs$values, runs$lengths, seq_along(runs$values) == length(runs$values)), use.names = FALSE)

body_lines <- character(0)
for (pn in unique(tab$panel)) {
  rows <- as.matrix(tab[tab$panel == pn, -1, drop = FALSE])
  body_lines <- c(body_lines, "\\hline", tex_row(sprintf("\\multicolumn{%d}{c}{%s}", n_all, pn)), "\\hline",
                  unname(apply(rows, 1, tex_row)))
}

tabular <- c(
  sprintf("\\begin{tabular}{%s}", spec),
  "\\hline",
  if (any(runs$values != "")) c(tex_row(c(rep("", n_lab), group_cells)), "\\hline"),
  tex_row(c("$p$-Hacking", "$\\Pi_0$", extra_hdr, tests$tex)),
  body_lines,
  "\\bottomrule",
  "\\end{tabular}"
)
if (shrink_to_fit) {
  tabular <- c("\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%",
               head(tabular, -1), "\\end{tabular}}")
}

# Note, built from the saved settings
col_of    <- function(v) n_lab + match(v, tests$col)       # position of a test's column in the table
hack_note <- c(
  threshold = sprintf(paste("Under threshold $p$-hacking, researchers report the first $t$-statistic $T_1$ they",
                            "draw if $T_1 > %s$ or the maximum of two $t$-statistics $\\max\\{T_1, T_2\\}$ otherwise."),
                      uniq(res$cv)),
  max       = paste("Under maximization $p$-hacking, for each $H$, the researcher draws eight independent",
                    "$t$-statistics with $T \\mid H \\sim N(H, 1)$ and reports only the largest one.")
)
panels <- unique(tab$panel)
note   <- c(
  if (length(panels) == 2) {
    paste0("The upper panel corresponds to size (no $p$-hacking) and the lower panel to power",
           if (is.na(what)) "." else paste0(" (", what, ")."))
  } else if (panels == "Size") "All rows report size (no $p$-hacking)." else "All rows report power.",
  hack_note[intersect(hacked, names(hack_note))],
  sprintf(paste("Column %d is the rejection rate of the projection test, columns %d--%d are the rejection",
                "rates of the six tests from %s, and columns %d--%d are the rejection rates of the three",
                "tests of %s. All tests are run at the 5\\%% level."),
          col_of("reject_rate"), col_of("rej_lcm_EWK"), col_of("CS2B_EWK"), ekw_cite,
          col_of("Kudrin_T_inf"), col_of("Kudrin_T_H_hybrid"), kudrin_cite),
  sprintf("The sample size is $n = %s$.", uniq(res$n)),
  sprintf("For the tests in %s we use their values of the tuning parameters and their code file Tests.R.", ekw_cite),
  # kappa = 3 is the default of kudrin_tests(), which run_sims() does not change
  sprintf(paste("For the tests of %s we use the sinc kernel with $\\kappa = 3$, and the critical value of",
                "$T_{\\infty}$ comes from %s nonparametric bootstrap draws."), kudrin_cite, uniq(res$nboots)),
  # create_weights(): w_j = min(1, (1/p)^(j - J0)) for j = 0, ..., num_coeffs - 1
  sprintf("For the projection test we set $w_j = \\min\\{1, %s^{j - %s}\\}$ and $J = %s$.",
          uniq(1 / res$p), uniq(res$J0), uniq(res$num_coeffs - 1)),
  # create_basis_stable() adds the points 0 and 9999 to the grid
  sprintf(paste("We form $\\mathcal{X}$ using a grid of %s elements evenly spaced between $\\pm %s$ plus the",
                "points 0 and 9999."), uniq(res$numgrid), uniq(res$L)),
  # run_sims() passes cv, not the shift column, to run_test_weight_fast()
  sprintf("We symmetrize about zero and shift by %s as described in Remarks 12 and 13.", uniq(res$cv)),
  sprintf("We run %s bootstrap repetitions and %s simulation repetitions.", uniq(res$nboots), uniq(res$nsims)),
  "Projection used \\texttt{nnls} package in R by Mullen and van Stokkum (2024).",
  "Transparency R package available at: \\url{https://github.com/sfaridan/When_is_p_hacking_detectable}"
)

label <- paste0("tab:sims_", if (length(hacked)) paste(hacked, collapse = "_") else "size")
tex   <- c(
  sprintf("%% Made by make_sim_table.R from %s on %s", basename(results_file), format(Sys.time(), "%Y-%m-%d %H:%M")),
  "% Needs booktabs (\\bottomrule), url or hyperref (\\url) and, with shrink_to_fit, graphicx (\\resizebox)",
  "\\begin{table}[htbp]",
  "\\centering",
  sprintf("\\caption{%s}", table_caption),
  sprintf("\\label{%s}", label),
  tabular,
  "",
  "\\medskip",
  "\\begin{minipage}{\\linewidth}",
  "\\footnotesize",
  strwrap(paste(note, collapse = " "), width = 100),
  "\\end{minipage}",
  "\\end{table}"
)

if (is.null(table_file)) {
  table_file <- file.path(dirname(results_file), paste0("table_", sub("\\.RDS$", ".tex", basename(results_file))))
}
writeLines(tex, table_file)
message("LaTeX table written to ", table_file)