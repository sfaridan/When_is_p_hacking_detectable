#### =====================================================================================================
#### cochrane_exhibits.R
#### The Cochrane exhibits for the paper, from what application_cochrane.R saved. It needs nothing in memory: run
#### application_cochrane.R, clear R, then run this. Writes, in res_dir:
####   table_cochrane<tag>.tex                  each test by de-rounding treatment, one panel per subsample
####                                            (needs \usepackage{booktabs, threeparttable})
####   Figure_cochrane_tinf_gap<tag>.pdf/.png   where T_inf's rejection comes from: the empirical CDF of p-values
####                                            minus Kudrin's deconvolution estimate, printed and de-rounded
####   Figdata_cochrane_tinf_gap<tag>.csv       the figure's curves
#### =====================================================================================================
rm(list = ls())

#### 0. Settings -------------------------------------------------------------------------------------------
root    <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"
res_dir <- file.path(root, "results", "cochrane")
tag     <- "_n100_all_trial_kappa0.5_zcounts_cont"   #the tag on application_cochrane.R's files (end of its sample report)
table_groups <- c("all", "early", "late")    #the table's panels, in order
show_sd      <- TRUE                         #standard deviations across draws, in parentheses, under the noisy treatments
show_J0_30   <- FALSE                        #a column for the projection test with J0 = 30, if it was run
kudrin_kappa <- 3                            #T_inf's frequency cutoff, as in kudrin_tests.R
fig_group    <- "all"                        #the subsample in the figure
fig_draw     <- 1                            #the uniform de-rounding draw in the figure

#### 1. Read what application_cochrane.R saved ---------------------------------------------------------------
res_file    <- file.path(res_dir, paste0("Cochrane_results", tag, ".rds"))
sample_file <- file.path(res_dir, paste0("Cochrane_sample", tag, ".rds"))
for (f in c(res_file, sample_file)) {
  if (!file.exists(f)) stop("No ", f, ": run application_cochrane.R first, or set tag to the one on its files")
}
RES <- readRDS(res_file)
SMP <- readRDS(sample_file)
st  <- RES$settings
d   <- SMP$sample
if (nrow(d) != nrow(RES$sample) || !isTRUE(all.equal(d$z, RES$sample$z))) {
  stop("The sample file and the results come from different runs (an inspect_only run since?): rerun application_cochrane.R")
}
cat(sprintf("Results saved %s: z_from = %s, kappa_cap = %s, treatments %s, %d draws, J0 %s, clustered by %s\n",
            st$saved, st$z_from, if (is.null(st$kappa_cap)) "none" else format(st$kappa_cap),
            paste(st$treatments, collapse = ", "), as.integer(st$n_draws), paste(st$J0_list, collapse = ", "), st$cluster_by))

dr  <- RES$draws                                               #one row per treatment, draw and subsample
ekw <- intersect(c("LCM", "Fisher", "Disc", "Binomial", "CS1", "CS2B"), names(dr))
dr$EKW_min <- if (length(ekw)) {
  apply(dr[, ekw, drop = FALSE], 1, function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))
} else NA_real_
cell_of <- function(tr, g, cn, f) {
  v <- dr[[cn]][dr$treatment == tr & dr$group == g]
  if (!length(v) || all(is.na(v))) NA_real_ else f(v[is.finite(v)])
}
mean_of <- function(tr, g, cn) cell_of(tr, g, cn, mean)
sd_of   <- function(tr, g, cn) cell_of(tr, g, cn, function(v) if (length(v) > 1) stats::sd(v) else NA_real_)
draws_n <- function(tr, g) sum(dr$treatment == tr & dr$group == g)
# EKW: the smallest of the six tests' averaged p-values, as in the Brodeur table (application_MM_analysis.R), with
# the standard deviation across draws of the test that attains it
ekw_min_of <- function(tr, g) {
  means <- vapply(ekw, function(cn) mean_of(tr, g, cn), 0)
  if (!length(ekw) || all(!is.finite(means))) return(c(p = NA_real_, sd = NA_real_))
  c(p = min(means, na.rm = TRUE), sd = sd_of(tr, g, ekw[which.min(means)]))
}
val_of <- function(tr, g, cn) if (cn == "EKW_min") ekw_min_of(tr, g)[["p"]] else mean_of(tr, g, cn)
dev_of <- function(tr, g, cn) if (cn == "EKW_min") ekw_min_of(tr, g)[["sd"]] else sd_of(tr, g, cn)
res_row <- function(tr, g) RES$results[RES$results$treatment == tr & RES$results$group == g, , drop = FALSE]

#### 2. The table ------------------------------------------------------------------------------------------
fp <- function(x) {                                            #a p-value: .974, <.001; bold when at most .05
  if (!is.finite(x)) return("")
  s <- if (x < 0.0005) "$<$.001" else sub("^0\\.", ".", sprintf("%.3f", x))
  if (x <= 0.05) paste0("\\textbf{", s, "}") else s
}
fsd <- function(x) if (!is.finite(x)) "" else paste0("(", sub("^0\\.", ".", sprintf("%.3f", x)), ")")
fB  <- function(x, bold = FALSE) {                             #B-hat or an allowance, x 10^3
  if (!is.finite(x)) return("")
  s <- sprintf("%.2f", 1e3 * x)
  if (bold) paste0("\\textbf{", s, "}") else s
}

TR_TEX <- c(none = "None", uniform = "Uniform", smin = "$s_{\\min}$", smin_noise = "$s_{\\min}$ + noise")
noisy  <- c("uniform", "smin_noise")
smin_t <- c("smin", "smin_noise")
trs    <- intersect(names(TR_TEX), st$treatments)
gr     <- intersect(table_groups, unique(RES$results$group))
if (!length(gr)) stop("None of table_groups is in the results: they hold ", paste(unique(RES$results$group), collapse = ", "))
p30    <- show_J0_30 && "proj_p_J0_30" %in% names(dr)
pcols  <- c("proj_p_J0_20", if (p30) "proj_p_J0_30")
kcols  <- intersect(c("T_inf", "T_H"), names(dr))
ocols  <- c(if (any(is.finite(dr$EKW_min))) "EKW_min", if ("Stouffer" %in% names(dr)) "Stouffer")
n_proj <- length(pcols) + 2                                    #p-value(s), 10^3 B-hat, allowance
ncols  <- 1 + n_proj + length(kcols) + length(ocols)

hdr1 <- paste0(" & \\multicolumn{", n_proj, "}{c}{Projection test}",
               if (length(kcols)) paste0(" & \\multicolumn{", length(kcols), "}{c}{Kudrin (2025)}") else "",
               if ("EKW_min" %in% ocols) " & EKW" else "", if ("Stouffer" %in% ocols) " & $p$-curve" else "", " \\\\")
rules <- paste0("\\cmidrule(lr){2-", 1 + n_proj, "}",
                if (length(kcols)) sprintf("\\cmidrule(lr){%d-%d}", 2 + n_proj, 1 + n_proj + length(kcols)) else "")
hdr2 <- paste0(paste(c("De-rounding",
                       if (p30) c("$p$, $J_0 = 20$", "$p$, $J_0 = 30$") else "$p$",
                       "$10^3 \\hat B$", "Allowance",
                       c(T_inf = "$T_\\infty$", T_H = "$T_H$")[kcols],
                       c(EKW_min = "smallest $p$", Stouffer = "Stouffer")[ocols]), collapse = " & "), " \\\\")

row_main <- function(tr, g) {
  rr   <- res_row(tr, g)
  smin <- tr %in% smin_t && nrow(rr) > 0
  cert <- smin && isTRUE(rr$certified_J0_20[1] == 1)
  cells <- c(TR_TEX[[tr]],
             vapply(pcols, function(cn) fp(mean_of(tr, g, cn)), ""),
             fB(mean_of(tr, g, "proj_B_J0_20"), bold = cert),
             if (smin) fB(rr$allowance_J0_20[1]) else "",
             vapply(c(kcols, ocols), function(cn) fp(val_of(tr, g, cn)), ""))
  paste0(paste(cells, collapse = " & "), " \\\\")
}
row_sd <- function(tr, g) {
  sdB   <- sd_of(tr, g, "proj_B_J0_20")
  cells <- c("",
             vapply(pcols, function(cn) fsd(sd_of(tr, g, cn)), ""),
             if (is.finite(sdB)) sprintf("(%.2f)", 1e3 * sdB) else "",
             "",
             vapply(c(kcols, ocols), function(cn) fsd(dev_of(tr, g, cn)), ""))
  paste0(paste(cells, collapse = " & "), " \\\\")
}
panel_head <- function(g, k) {
  rr <- RES$results[RES$results$group == g, , drop = FALSE][1, ]
  sprintf("\\multicolumn{%d}{l}{\\textit{Panel %s. %s: %s outcomes, %s %s}} \\\\", ncols, LETTERS[k], rr$label,
          format(rr$n, big.mark = ","), format(rr$clusters, big.mark = ","),
          if (identical(st$cluster_by, "trial")) "trials" else "clusters")
}
body <- unlist(lapply(seq_along(gr), function(k) {
  g <- gr[k]
  c(if (k > 1) "\\addlinespace", panel_head(g, k),
    unlist(lapply(trs, function(tr) {
      c(row_main(tr, g), if (show_sd && tr %in% noisy && draws_n(tr, g) > 1) row_sd(tr, g))
    })))
}))

notes <- paste0(
  sprintf("Continuous efficacy outcomes of randomized trials with at least %d participants in the Cochrane Database of Systematic Reviews (the data of van Zwet, Schwab and Senn 2021); ", as.integer(st$n_min)),
  if (identical(st$z_from, "counts")) "$z$ is the mean difference over its standard error, from the printed arm means and standard deviations" else "$z$ is the database's own $z$",
  if (!is.null(st$kappa_cap)) sprintf("; outcomes whose rounding unit in $z$-units, $\\kappa = 2a/\\mathrm{SE}_{\\min}$, exceeds %s are left out. ", format(st$kappa_cap)) else ". ",
  "De-rounding: None uses the printed values; Uniform draws each mean and standard deviation uniformly within its rounding cell; $s_{\\min}$ rescales $z$ to the smallest standard error consistent with the printed standard deviations; $s_{\\min}$ + noise also adds uniform noise at each mean's last digit. ",
  sprintf("Under Uniform and $s_{\\min}$ + noise, entries are averages over %d draws%s. ", as.integer(st$n_draws),
          if (show_sd) ", with standard deviations across draws in parentheses" else ""),
  "Projection test: $J = 40$, $J_0 = 20$, $p = 2$, centered at 1.96. Allowance: the rounding allowance $0.13\\,\\overline{\\kappa^6}$, times $10^3$ (Propositions 3 and 4); under $s_{\\min}$ and $s_{\\min}$ + noise, the projection test rejects only when $\\hat B$ exceeds it. ",
  "Every test clusters by trial. EKW: the smallest of the six tests' $p$-values of Elliott, Kudrin and W\\\"uthrich (2022), each averaged over draws as above; not itself a valid $p$-value. $p$-curve: Stouffer's test for left skew (Simonsohn, Nelson and Simmons 2014). Bold: $p$-values of at most .05, and $10^3 \\hat B$ above the allowance.")

tex <- c(sprintf("%% Written by cochrane_exhibits.R from Cochrane_results%s.rds (saved %s)", tag, st$saved),
         "% Needs \\usepackage{booktabs, threeparttable}",
         "\\begin{table}[!htbp]\\centering",
         "\\begin{threeparttable}",
         "\\caption{Tests for selective reporting in the Cochrane Database of Systematic Reviews: continuous efficacy outcomes}",
         "\\label{tab:cochrane}",
         sprintf("\\begin{tabular}{l%s}", strrep("c", ncols - 1)),
         "\\toprule", hdr1, rules, hdr2, "\\midrule", body, "\\bottomrule",
         "\\end{tabular}",
         "\\begin{tablenotes}[flushleft]\\footnotesize",
         paste0("\\item ", notes),
         "\\end{tablenotes}",
         "\\end{threeparttable}",
         "\\end{table}")
table_file <- file.path(res_dir, paste0("table_cochrane", tag, ".tex"))
writeLines(tex, table_file)
cat("\n", paste(tex, collapse = "\n"), "\n\nTable written to ", table_file, "\n", sep = "")

#### 3. The figure: the gap T_inf takes the supremum of -----------------------------------------------------
# Kudrin's deconvolution estimate of the p-value CDF with the sinc kernel (his equation 7),
# G_tilde(p) = mean_s [1 - (Si(kappa (cv(p) - |t_s|)) + Si(kappa (cv(p) + |t_s|))) / pi], cv(p) = Phi^-1(1 - p/2),
# against the empirical CDF G_hat(p) = share of |t_s| >= cv(p). The grid is in |t|, so both ends of p are dense.

# The uniform de-rounding draw r of application_cochrane.R's z_treat("uniform", r), from the saved sample and seed
z_uniform <- function(r) {
  z <- d$z
  i <- which(d$source == "continuous" & is.finite(d$z_raw))
  if (!length(i)) return(z)
  set.seed(st$seed + r)
  U <- matrix(runif(4 * length(i)), ncol = 4)
  a <- d$a_mean[i]; b <- d$b_sd[i]
  sd_draw <- function(s, u) { lo <- pmax(s - b / 2, 0); lo + (s + b / 2 - lo) * u }
  m1 <- d$m1[i] + a * (U[, 1] - 0.5); m2 <- d$m2[i] + a * (U[, 2] - 0.5)
  s1 <- sd_draw(d$s1[i], U[, 3]);     s2 <- sd_draw(d$s2[i], U[, 4])
  zc <- (m1 - m2) / sqrt(s1^2 / d$n1[i] + s2^2 / d$n2[i])
  z[i] <- ifelse(is.finite(zc), z[i] + (zc - d$z_raw[i]), z[i])
  z
}

si_x <- seq(0, 200, by = 0.001)                                #the sine integral by table up to 200,
si_f <- c(1, sin(si_x[-1]) / si_x[-1])                         #the asymptotic series beyond
si_v <- c(0, cumsum((si_f[-1] + si_f[-length(si_f)]) / 2) * 0.001)
Si <- function(x) {
  ax  <- abs(x)
  out <- approx(si_x, si_v, pmin(ax, 200))$y
  big <- ax > 200
  out[big] <- pi / 2 - cos(ax[big]) / ax[big] - sin(ax[big]) / ax[big]^2
  sign(x) * out
}
t_grid <- seq(0.002, 3.3, by = 0.002)                          #|t|; p from .998 down to .001
p_grid <- 2 * pnorm(-t_grid)
gap_curve <- function(z) {
  a  <- abs(z[is.finite(z)])
  Gt <- vapply(t_grid, function(cv) mean(1 - (Si(kudrin_kappa * (cv - a)) + Si(kudrin_kappa * (cv + a))) / pi), 0)
  Gh <- vapply(t_grid, function(cv) mean(a >= cv), 0)
  Gh - Gt
}

grp <- list(all = rep(TRUE, nrow(d)),
            early = is.finite(d$year) & d$year <= st$year_split,
            late = is.finite(d$year) & d$year > st$year_split,
            dichotomous = d$source == "dichotomous",
            continuous = d$source == "continuous")
if (!fig_group %in% names(grp)) stop("fig_group must be one of ", paste(names(grp), collapse = ", "))
s <- grp[[fig_group]]
has_unif <- "uniform" %in% st$treatments
cat("\nComputing the gap curves (a minute or so)...\n")
gap_printed <- gap_curve(d$z[s])
gap_unif    <- if (has_unif) gap_curve(z_uniform(fig_draw)[s]) else rep(NA_real_, length(t_grid))

lab_p <- function(x) if (!is.finite(x)) "n/a" else if (x < 0.0005) "< .001" else sub("^0\\.", ".", sprintf("%.3f", x))
p_none <- mean_of("none", fig_group, "T_inf")
p_draw <- {
  v <- dr$T_inf[dr$treatment == "uniform" & dr$group == fig_group & dr$draw == fig_draw]
  if (length(v)) v[1] else NA_real_
}
p_unif <- mean_of("uniform", fig_group, "T_inf")
n_unif <- draws_n("uniform", fig_group)

where <- function(g, lab) {
  if (all(is.na(g))) return(invisible(NULL))
  i <- which.max(abs(g))
  j <- which(p_grid <= 0.15); k <- j[which.max(abs(g[j]))]
  cat(sprintf("%-34s largest gap %+.4f at p = %.3f (|t| = %.2f); largest at p <= .15: %+.4f at p = %.3f\n",
              lab, g[i], p_grid[i], t_grid[i], g[k], p_grid[k]))
}
where(gap_printed, "Printed values:")
where(gap_unif, sprintf("Uniform de-rounding, draw %d:", fig_draw))

fig  <- file.path(res_dir, paste0("Figure_cochrane_tinf_gap", tag))
yl   <- range(c(gap_printed, gap_unif), na.rm = TRUE) * 1.08
leg  <- c(as.expression(bquote("Printed values: " * T[infinity] * " p = " * .(lab_p(p_none)))),
          if (has_unif) as.expression(bquote("Uniform de-rounding, draw " * .(fig_draw) * ": " * T[infinity] * " p = " *
                                               .(lab_p(p_draw)) * " (mean over " * .(n_unif) * " draws " * .(lab_p(p_unif)) * ")")))
for (ext in c("pdf", "png")) {
  if (ext == "pdf") grDevices::pdf(paste0(fig, ".pdf"), width = 5.6, height = 3.6, pointsize = 9)
  else grDevices::png(paste0(fig, ".png"), width = 5.6, height = 3.6, units = "in", res = 300, pointsize = 9)
  op <- par(mar = c(3.4, 4.6, 1.6, 0.8), mgp = c(2.3, 0.6, 0), las = 1)
  plot(p_grid, gap_printed, type = "n", log = "x", xlim = c(0.001, 1), ylim = yl, axes = FALSE,
       xlab = "p-value (log scale)", ylab = "Empirical CDF minus deconvolution estimate")
  axis(1, at = c(0.001, 0.005, 0.01, 0.05, 0.15, 0.5, 1), labels = c(".001", ".005", ".01", ".05", ".15", ".5", "1"))
  axis(2); box()
  abline(h = 0, lty = 3, col = "grey60")
  abline(v = c(0.05, 0.15), lty = 2, col = "grey40")
  mtext(c(".05", ".15"), side = 3, at = c(0.05, 0.15), line = 0.2, cex = 0.75, col = "grey30")
  lines(p_grid, gap_printed, col = "grey55", lwd = 1.2)
  if (has_unif) lines(p_grid, gap_unif, col = "black", lwd = 1.3)
  legend("bottomleft", bty = "n", legend = leg, lty = 1, lwd = c(1.2, 1.3)[seq_along(leg)],
         col = c("grey55", "black")[seq_along(leg)], cex = 0.9)
  par(op)
  grDevices::dev.off()
}
write.csv(data.frame(t = t_grid, p = p_grid, gap_printed = gap_printed, gap_uniform = gap_unif),
          file.path(res_dir, paste0("Figdata_cochrane_tinf_gap", tag, ".csv")), row.names = FALSE)
cat("Figure written to ", fig, ".pdf/.png; its curves to Figdata_cochrane_tinf_gap", tag, ".csv\n", sep = "")