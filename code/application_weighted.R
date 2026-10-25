# 1. Clear all objects from the global environment
rm(list = ls())

# 2. Force R to release the freed memory to your OS
gc(reset = TRUE)

#load in functions
root               <- "C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/sfaridani6/Documents/Research/smoothness test for p-hacking"
this_script        <- rstudioapi::getSourceEditorContext()$path

source(paste0(root,"/code/functions_all_weighted.R"))

#Load in the data
data_location <- "C:/Users/stefa/OneDrive/Documents/R/Underpowered Literatures/data"
setwd(data_location)
MM_raw <- haven::read_dta('MM data.dta')
#MM_raw <- haven::read_dta('MM data fromkudrin.dta')
MM_data <- MM_raw

#Kudrin's (2025) drop
drop          <- is.na(MM_raw$mu) | MM_raw$mu == 0
MM_dropped    <- MM_raw[!drop, ]
MM_data       <- MM_dropped

#### Prep data

#Set the seed (for consistent de-rounding)
set.seed(1) 

#de-round via mus and sd
MM_data$mu    <- de_round(MM_data$mu)
MM_data$sd    <- de_round(MM_data$sd)
MM_data$t_der <- MM_data$mu/MM_data$sd
MM_data$t[!is.na(MM_data$t_der)] <- MM_data$t_der[!is.na(MM_data$t_der)] #replace t with derounded whenever it exists

#de-round the raw t-score itself when mu, sd are not available
MM_data$t[is.na(MM_data$t_der)] <- de_round(MM_data$t[is.na(MM_data$t_der)])

if (sum(is.infinite(MM_data$t))+sum(is.na(MM_data$t))+sum(is.nan(MM_data$t)) > 0){
  stop("Tscores not clean")
}

### Run the tests

cv         <- 1.96
sigma_Y    <- 1 
shift      <- cv
L          <- 6.5 #grid support [-L,L] \cup 9999
numgrid    <- 3000 #total number of elements in grid (density)
boots      <- 1000 #bootstrap repetitions
J          <- 40
coeffs     <- J+1
p          <- 2
epsilon_U_pre <- compute_epsilons(L, nx = numgrid)$epsilon_U #pre-compute for speed

#Run the projection test with more coefficients

for (J0 in c(20)){
  rct_results <- run_test_weight_fast(MM_data[MM_data$method=="RCT",],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
  did_results <- run_test_weight_fast(MM_data[MM_data$method=="DID",],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
  rdd_results <- run_test_weight_fast(MM_data[MM_data$method=="RDD",],numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
  iv_results  <- run_test_weight_fast(MM_data[MM_data$method=="IV",] ,numcoeffs=coeffs,epsilon_U_pre=epsilon_U_pre,sigma_Y=1,shift=shift,L=L,numgrid=numgrid,boots=boots,J0=J0,p=p)
  
  #Save results
  empirical_results <- as.data.frame(do.call(rbind, list(rct_results,iv_results, did_results,rdd_results )))
  empirical_results$methods <- c("RCT","IV","DID","RDD")
  empirical_results<- empirical_results[,c("methods","n","articles","pval","breakdown","numcoeffs","J0")]
  setwd(paste0(root,"/results"))
  name <- paste0("MM_weights_results_",J0,"_J0")
  saveRDS(empirical_results, file = name)
  print(readRDS(name))
}
print(c(rct_results$runtime,rct_results$init_time))

#Run EWK Tests

source(paste0(root,"/code/Tests.R"))

run_ewk <- function(ps,ids=1:length(ps),pmax=0.15){
  pmin <- min(ps[ps>0])
  out <- list()
  
  out$lcms_EWK     <- LCM(ps,pmin,pmax)
  out$Fisher_EWK   <- Fisher(ps,pmin,pmax)
  out$disconts_EWK <- Discontinuity_test(ps,.05)
  out$binomial_EWK <- Binomial(ps,.04, .05, "c")
  out$CS1_EWK_s    <- CoxShi(ps, ids, pmin, pmax, 15, 1, 0)
  out$CS2_EWK_s    <- CoxShi(ps, ids, pmin, pmax, 15, 2, 1)
  
  out$min <- min(out$lcms_EWK,out$Fisher_EWK,out$disconts_EWK, out$binomial_EWK,  out$CS1_EWK_s, out$CS2_EWK_s)
  
  return(out)
}
rct     <- MM_data$method == "RCT"
ewk_rct <- run_ewk(2 * pnorm(-abs(MM_data$t[rct])), ids = as.integer(factor(MM_data$title[rct])))
iv      <- MM_data$method == "IV"
ewk_iv  <- run_ewk(2 * pnorm(-abs(MM_data$t[iv])),  ids = as.integer(factor(MM_data$title[iv])))
did     <- MM_data$method == "DID"
ewk_did <- run_ewk(2 * pnorm(-abs(MM_data$t[did])), ids = as.integer(factor(MM_data$title[did])))
rdd     <- MM_data$method == "RDD"
ewk_rdd <- run_ewk(2 * pnorm(-abs(MM_data$t[rdd])), ids = as.integer(factor(MM_data$title[rdd])))
empirical_ewk <-cbind(ewk_rct,ewk_iv,ewk_did,ewk_rdd)
setwd(paste0(root,"/results"))
saveRDS(empirical_ewk, file = paste0("EWK_results"))
print(empirical_ewk)

#Run Kudrin (2025) Tests: T_inf, T_H, T_H_hybrid

source(paste0(root,"/code/kudrin_tests.R"))
#kudrin_tests_selfcheck() #optional one-time check of kudrin_tests.R

run_kudrin <- function(data){
  kt <- kudrin_tests(data$t, cluster = data$title, boots = boots) #cluster by article, as in the projection test
  data.frame(n = length(data$t), articles = length(unique(data$title)), as.list(kt$pval),
             min = min(kt$pval), restricted = kt$restricted) #restricted = FALSE means T_H_hybrid = T_H
}
run_kudrin_nocluster <- function(data){
  kt <- kudrin_tests(data$t, cluster = NULL, boots = boots) #cluster by article, as in the projection test
  data.frame(n = length(data$t), articles = length(unique(data$title)), as.list(kt$pval),
             min = min(kt$pval), restricted = kt$restricted) #restricted = FALSE means T_H_hybrid = T_H
}


set.seed(1) #T_inf is bootstrapped; reset the seed so these results don't depend on the code above
kudrin_rct <- run_kudrin(MM_data[MM_data$method=="RCT",])
kudrin_iv  <- run_kudrin(MM_data[MM_data$method=="IV",])
kudrin_did <- run_kudrin(MM_data[MM_data$method=="DID",])
kudrin_rdd <- run_kudrin(MM_data[MM_data$method=="RDD",])

empirical_kudrin <- cbind(methods = c("RCT","IV","DID","RDD"), rbind(kudrin_rct,kudrin_iv,kudrin_did,kudrin_rdd))
setwd(paste0(root,"/results"))
saveRDS(empirical_kudrin, file = "Kudrin_results")
print(readRDS("Kudrin_results"))

#### Table: projection test, EKW and Kudrin (2025) tests (same layout as Table 3 in the paper)
# Reads the results files saved above (so after one full run you can tweak and re-run this block alone),
# adds a column group for Kudrin's three tests, prints the table and writes results/table_application.tex.
# LaTeX needs booktabs, graphicx (shrink-to-fit wrapper) and hyperref or url (for \url in the notes).

make_application_table <- function(J0_values = c(20),            # one (p, B-hat) group per J0, in this order; each must be in the J0 loop above
                                   res_dir   = paste0(root, "/results"),
                                   out_file  = "table_application.tex") {
  meth <- c("RCT", "IV", "DID", "RDD")
  
  # Effective sample sizes from the hand-collected sample of 100 articles per method (typed in from the paper)
  med_nu   <- c(RCT = 185,    IV = 262,    DID = 293.5,  RDD = NA)
  Delta_nu <- c(RCT = 0.0014, IV = 0.0013, DID = 0.0014, RDD = NA)
  
  # Citations used in the header and notes; swap in "\\citet{key}" to get hyperlinked citations
  cite <- c(brodeur = "Brodeur et al. (2020)", ekw  = "Elliott et al. (2022)", kudrin = "Kudrin (2025)",
            coxshi  = "Cox and Shi (2023)",    osqp = "Stellato et al. (2024)")
  
  num  <- function(x) vapply(x, function(v) as.numeric(v)[1], numeric(1), USE.NAMES = FALSE) #results files hold list-columns
  pick <- function(d, what) {                                    #rows of d in the order of meth
    i <- match(meth, d$methods)
    if (anyNA(i)) stop(what, " has no results for ", paste(meth[is.na(i)], collapse = ", "))
    d[i, ]
  }
  
  #--- Read the saved results
  proj <- lapply(J0_values, function(j0) {
    f <- file.path(res_dir, paste0("MM_weights_results_", j0, "_J0"))
    if (!file.exists(f)) stop("No projection results for J0 = ", j0, ": add it to the J0 loop above")
    r <- pick(readRDS(f), basename(f))
    if (any(num(r$numcoeffs) != coeffs)) warning(basename(f), " was run with numcoeffs != ", coeffs, " (stale file?)")
    data.frame(n = num(r$n), articles = num(r$articles), p = num(r$pval), B = num(r$breakdown))
  })
  ewk <- readRDS(file.path(res_dir, "EWK_results"))              #list-matrix: tests x methods
  ewk <- matrix(num(ewk), nrow(ewk), dimnames = dimnames(ewk))[, paste0("ewk_", tolower(meth))]
  kud <- pick(readRDS(file.path(res_dir, "Kudrin_results")), "Kudrin_results")
  
  n_all <- cbind(sapply(proj, function(d) d$n),        num(kud$n))   #same samples in every file?
  a_all <- cbind(sapply(proj, function(d) d$articles), num(kud$articles))
  if (any(n_all != n_all[, 1]) || any(a_all != a_all[, 1]))
    warning("n or article counts differ across results files: were they all run on the same MM_data?")
  
  #--- Numbers
  tab <- data.frame(n = proj[[1]]$n, articles = proj[[1]]$articles,
                    med_nu = med_nu[meth], Delta_nu = Delta_nu[meth], row.names = meth)
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
  cells <- cbind(meth, sprintf("%.0f", tab$n), sprintf("%.0f", tab$articles),
                 ifelse(is.na(tab$med_nu), "$\\cdot$", as.character(tab$med_nu)), fmt(tab$Delta_nu, 4, lead0 = TRUE))
  for (j0 in J0_values) cells <- cbind(cells, fmt(tab[[paste0("p_J0_", j0)]], 2), fmt(tab[[paste0("B_J0_", j0)]], 4))
  for (k in c(names(ekw_rows), names(kud_cols))) cells <- cbind(cells, fmt(tab[[k]], 2))
  unr <- !tab$hybrid_restricted                                  #dagger: hybrid restriction did not bind
  if (any(unr)) cells[unr, ncol(cells) - 1] <- paste0(cells[unr, ncol(cells) - 1], "$^{\\dagger}$")
  
  #--- Notes (one sentence per line in the .tex; edit freely)
  rdd_insig <- isTRUE(all(unlist(tab["RDD", paste0("p_J0_", J0_values)]) > 0.05))
  if (!rdd_insig) warning("RDD projection-test p-value is now <= .05, so the notes' reason for omitting its Delta_nu no longer holds")
  notes <- c(
    paste0("Results for the ", cite[["brodeur"]], " data. Following ", cite[["kudrin"]], ", we drop tests whose reported coefficient is missing or zero."),
    paste0("Columns headed $J_0$ report $p$-values and breakdown statistics $\\widehat{B}$ of the projection test with $J=", J, "$",
           if (length(J0_values) > 1) "; several values of $J_0$ are included to show low sensitivity to the smoothing parameter", "."),
    "Med.\\ $\\nu$ is the median effective sample size as measured by the median number of clusters in a random sample of 100 articles for each causal inference method.",
    "In ambiguous cases, the smallest cluster count was always used.",
    paste0("$\\Delta_\\nu=\\frac{1}{100}\\sum_{k=1}^{100}\\sup_h\\|\\varphi(t-h)-t_{\\nu_k,h}(t)\\|_T$ is the average difference between the PDFs of the standard normal ",
           "and Student-$t$ over the sample of 100 effective sample sizes collected (computed numerically)."),
    if (rdd_insig) "Since the $p$-value for RDDs is insignificant, we did not compute $\\Delta_\\nu$ for them.",
    sprintf("We form $\\mathcal{X}$ using a grid of %d elements evenly spaced between $\\pm %s$ plus the scalar 9999.", numgrid, L),
    sprintf("Discretization incurs error $\\mathcal{E}(L,\\delta)=%.5f$, which is already accounted for in $\\widehat{B}$ and the $p$-value.", epsilon_U_pre),
    "The discretization was chosen to be made of round numbers and result in error relatively small compared to $\\Delta_\\nu$.",
    paste0("EKW $p$-values are for the tests of ", cite[["ekw"]], ": the least concave majorant test (LCM) and Fisher's test (Fisher) on $p$-values in $(0,0.15]$, ",
           "the density discontinuity test at 0.05 (Disc.), the binomial test on $[0.04,0.05]$ (Binomial), and the ", cite[["coxshi"]],
           " histogram-based tests on $(0,0.15]$ with 15 bins for monotonicity (CS1) and for 2-monotonicity with bounds (CS2B)."),
    paste0(cite[["kudrin"]], " $p$-values are for his KS-type test ($T_\\infty$), histogram-based test ($T_H$) and hybrid histogram-based test ($T_H^{\\mathrm{hybrid}}$), ",
           "using the sinc kernel with $\\kappa=3$ and, for $T_H$ and $T_H^{\\mathrm{hybrid}}$, the points $0.01,0.02,\\ldots,0.15$."),
    "All three cluster by article: $T_H$ and $T_H^{\\mathrm{hybrid}}$ use cluster-robust covariance matrices and $T_\\infty$ uses an article-level bootstrap.",
    if (any(unr)) "$^\\dagger$The deconvolution estimate already satisfies the hybrid test's shape restrictions, so $T_H^{\\mathrm{hybrid}}=T_H$.",
    "Min is the smallest $p$-value in each group of tests.",
    sprintf("The projection test and $T_\\infty$ use %d bootstrap repetitions.", boots),
    paste0("Computation used osqp package in R by ", cite[["osqp"]], ". Data: ", cite[["brodeur"]], "."),
    "Transparency R package available at: \\url{https://github.com/sfaridan/When_is_p_hacking_detectable}")
  
  #--- LaTeX
  width <- c(rep(2, length(J0_values)), 7, 4)                    #column groups after the 5 leading columns
  last  <- 5 + cumsum(width)
  grp   <- c(sprintf("$J_0=%s$", J0_values), "EKW $p$-values", paste(cite[["kudrin"]], "$p$-values"))
  head2 <- c("", "$n$", "Articles", "Med.\\ $\\nu$", "$\\Delta_\\nu$", rep(c("$p$", "$\\widehat{B}$"), length(J0_values)),
             "LCM", "Fisher", "Disc.", "Binomial", "CS1", "CS2B", "Min",
             "$T_\\infty$", "$T_H$", "$T_H^{\\mathrm{hybrid}}$", "Min")
  stopifnot(length(head2) == ncol(cells), last[length(last)] == ncol(cells))
  tex_row <- function(x) paste0(paste(x, collapse = " & "), " \\\\")
  
  tex <- c("\\begin{table}[htbp]",
           "\\centering",
           "\\caption{Empirical Application Results}",
           "\\label{tab:application}",
           "\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{% shrinks the table only if it is wider than the text",
           paste0("\\begin{tabular}{l", strrep("c", ncol(cells) - 1), "}"),
           tex_row(c(rep("", 5), paste0("\\multicolumn{", width, "}{c}{", grp, "}"))),
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
  saveRDS(tab, file.path(res_dir, "Table_application_numbers"))
  
  #--- Console version
  txt <- gsub("$^{\\dagger}$", "+", gsub("$\\cdot$", "-", cells[, -1], fixed = TRUE), fixed = TRUE)
  dimnames(txt) <- list(meth, c("n", "Articles", "Med.nu", "Delta_nu", paste0(c("p_J0=", "B_J0="), rep(J0_values, each = 2)),
                                "LCM", "Fisher", "Disc.", "Binomial", "CS1", "CS2B", "EKW_min", "T_inf", "T_H", "T_H_hyb", "Kud_min"))
  print(noquote(txt), right = TRUE)
  if (any(unr)) cat("+ hybrid restriction did not bind, so T_H_hybrid = T_H\n")
  cat("LaTeX table written to", file.path(res_dir, out_file), "\n")
  invisible(list(numbers = tab, latex = tex))
}

table_application <- make_application_table(J0_values = c(20))   #e.g. c(30, 20) to mirror the paper's two columns
