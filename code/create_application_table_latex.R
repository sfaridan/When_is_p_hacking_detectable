############################################################
# Standalone script: make_brodeur_results_table.R
#
# Purpose:
#   Load saved empirical application results and generate the
#   LaTeX table tab:brodeur_results, now including EKW p-values.
#
# Expected input files:
#   results/MM_tangcone_results_20_J
#   results/MM_tangcone_results_30_J
#   results/EWK_results
#
# Output:
#   results/brodeur_results.tex
############################################################

# -----------------------------
# User settings
# -----------------------------
root <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
results_dir <- file.path(root, "results")
out_file <- file.path(results_dir, "brodeur_results.tex")

# -----------------------------
# Input files
# -----------------------------
file20 <- file.path(results_dir, "MM_tangcone_results_20_J")
file30 <- file.path(results_dir, "MM_tangcone_results_30_J")
file_ewk <- file.path(results_dir, "EWK_results")

if (!file.exists(file20)) stop("Missing file: ", file20)
if (!file.exists(file30)) stop("Missing file: ", file30)
if (!file.exists(file_ewk)) stop("Missing file: ", file_ewk)

# -----------------------------
# Load saved results
# -----------------------------
res20_raw <- readRDS(file20)
res30_raw <- readRDS(file30)
ewk_raw   <- readRDS(file_ewk)

# -----------------------------
# Keep only needed fields
# -----------------------------
prep_results <- function(df) {
  out <- data.frame(
    methods   = as.character(df$methods),
    n         = as.numeric(df$n),
    articles  = as.numeric(df$articles),
    pval      = as.numeric(df$pval),
    breakdown = as.numeric(df$breakdown),
    stringsAsFactors = FALSE
  )
  out
}

res20 <- prep_results(res20_raw)
res30 <- prep_results(res30_raw)

# Paper row order
row_order <- c("RCT", "IV", "DID", "RDD")
res20 <- res20[match(row_order, res20$methods), ]
res30 <- res30[match(row_order, res30$methods), ]

# Sanity checks
stopifnot(identical(res20$methods, row_order))
stopifnot(identical(res30$methods, row_order))
stopifnot(all(res20$n == res30$n))
stopifnot(all(res20$articles == res30$articles))

# -----------------------------
# Prepare EKW results
# Expected source object:
# empirical_ewk <- cbind(ewk_rct, ewk_iv, ewk_did, ewk_rdd)
# so methods are in order: RCT, IV, DID, RDD
# and rows should be:
#   lcms_EWK, Fisher_EWK, disconts_EWK, binomial_EWK,
#   CS1_EWK_s, CS2_EWK_s, min
# -----------------------------
prep_ewk <- function(x) {
  # Convert whatever came back into a matrix/data frame
  if (is.list(x) && !is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  
  # If saved as 7x4 matrix/list-matrix, rows are tests and cols are methods
  # If saved as 4x7, transpose it.
  if (nrow(x) == 7 && ncol(x) == 4) {
    ewk_df <- x
  } else if (nrow(x) == 4 && ncol(x) == 7) {
    ewk_df <- as.data.frame(t(as.matrix(x)), stringsAsFactors = FALSE)
  } else {
    stop("Unexpected dimensions for EWK_results: ", nrow(x), " x ", ncol(x))
  }
  
  # Standardize names
  rownames(ewk_df) <- c(
    "lcms_EWK", "Fisher_EWK", "disconts_EWK", "binomial_EWK",
    "CS1_EWK_s", "CS2_EWK_s", "min"
  )
  colnames(ewk_df) <- row_order
  
  # Convert all entries to numeric
  ewk_df[] <- lapply(ewk_df, function(col) as.numeric(unlist(col)))
  
  # Return method-by-row lookup frame
  out <- data.frame(
    methods      = row_order,
    LCM          = as.numeric(ewk_df["lcms_EWK", row_order]),
    Fisher       = as.numeric(ewk_df["Fisher_EWK", row_order]),
    Disc         = as.numeric(ewk_df["disconts_EWK", row_order]),
    Binomial     = as.numeric(ewk_df["binomial_EWK", row_order]),
    CS1          = as.numeric(ewk_df["CS1_EWK_s", row_order]),
    CS2          = as.numeric(ewk_df["CS2_EWK_s", row_order]),
    Min          = as.numeric(ewk_df["min", row_order]),
    stringsAsFactors = FALSE
  )
  
  out
}

ewk <- prep_ewk(ewk_raw)
ewk <- ewk[match(row_order, ewk$methods), ]
stopifnot(identical(ewk$methods, row_order))

# -----------------------------
# Manuscript-only columns
# (not computed in empirical_application.R)
# -----------------------------

Deltas <-  readRDS(paste0(root,"/results/delta_summary.RDS"))

med_nu <- c(
  RCT = Deltas$RCT[1], #"185",
  IV  = Deltas$IV[1], #"287",
  DID = Deltas$DID[1], #"261",
  RDD = "$\\cdot$"
)

delta_nu <- c(
  RCT = round(Deltas$RCT[2],4), #".0014",
  IV  = round(Deltas$IV[2],4), #".0014",
  DID = round(Deltas$DID[2],4), #".0012",
  RDD = "$\\cdot$"
)

# -----------------------------
# Formatting helpers
# -----------------------------
fmt_no_leading_zero <- function(x, digits) {
  sub("^0", "", sprintf(paste0("%.", digits, "f"), x))
}

fmt_p <- function(p) {
  fmt_no_leading_zero(p, 2)
}


fmt_Bhat <- function(breakdown, n, pval, alpha = 0.05) {
  if (is.na(breakdown) || is.na(n) || is.na(pval) || pval >= alpha || breakdown <= 0) {
    return("$\\cdot$")
  }
  fmt_no_leading_zero(breakdown , 4)
}

make_row <- function(method_name) {
  i20 <- match(method_name, res20$methods)
  i30 <- match(method_name, res30$methods)
  ie  <- match(method_name, ewk$methods)
  
  paste0(
    method_name,
    " & ", res30$n[i30],
    " & ", res30$articles[i30],
    " & ", med_nu[method_name],
    " & ", delta_nu[method_name],
    " & ", fmt_p(res30$pval[i30]),
    " & ", fmt_Bhat(res30$breakdown[i30], res30$n[i30], res30$pval[i30]),
    " && ",
    fmt_p(res20$pval[i20]),
    " & ", fmt_Bhat(res20$breakdown[i20], res20$n[i20], res20$pval[i20]),
    " && ",
    fmt_p(ewk$LCM[ie]),
    " & ", fmt_p(ewk$Fisher[ie]),
    " & ", fmt_p(ewk$Disc[ie]),
    " & ", fmt_p(ewk$Binomial[ie]),
    " & ", fmt_p(ewk$CS1[ie]),
    " & ", fmt_p(ewk$CS2[ie]),
    " & ", fmt_p(ewk$Min[ie]),
    "\\\\"
  )
}

# -----------------------------
# Build LaTeX table
# -----------------------------
table_lines <- c(
  "\\begin{table}[h!]",
  " \\begin{center}",
  " \\caption{\\label{tab:brodeur_results} Empirical Application Results}",
  " \\resizebox{\\textwidth}{!}{%",
  " \\begin{tabular}{lccccccccccccccccc}",
  " & & & & &\\multicolumn{2}{c}{J=30} && \\multicolumn{2}{c}{J=20} && \\multicolumn{7}{c}{EKW p-values} \\\\",
  " \\cline{6-7} \\cline{9-10} \\cline{12-18}",
  " & n & Articles & Med. $\\\\nu$ & $\\\\Delta_\\\\nu$ & $p$ & $\\\\widehat{B}$ & & $p$ & $\\\\widehat{B}$ & & LCM & Fisher & Disc. & Binomial & CS1 & CS2 & Min \\\\",
  " \\hline",
  make_row("RCT"),
  make_row("IV"),
  make_row("DID"),
  make_row("RDD"),
  " \\hline",
  " \\end{tabular}%",
  " }",
  " \\end{center}",
    "\\end{table}"
)

# -----------------------------
# Save and print
# -----------------------------
writeLines(table_lines, con = out_file)
cat(paste(table_lines, collapse = "\n"))
cat("\n\nSaved LaTeX table to:\n", out_file, "\n")