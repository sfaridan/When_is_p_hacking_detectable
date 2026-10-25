################################################################################
#  kudrin_cleaning.R
#
#  Data preparation for the Brodeur, Cook & Heyes (2020) "Methods Matter" data,
#  following Kudrin (2025) and the EKW replication code
#  (Step1_Methods_Matter_Data_Work.R, github.com/nvkudrin/pHackingPower).
#  Written for Kudrin's copy of the data (columns t, mu_orig, sd_orig, pv, report).
#
#  prep_kudrin(file)   Kudrin's sample with raw t-statistics: the input to his Table 2
#    * drops tests whose coefficient is missing or exactly zero, which gives his
#      5,780 / 5,158 / 7,101 / 3,117 tests (DID / IV / RCT / RDD) from 680 papers
#    * uses his raw p-value, p = 2 * (1 - pnorm(|t|)). That is exactly 0 for
#      |t| > 8.29, so those rows get |t| = 1e6: same p, and xi = 1 in kudrin_tests()
#    * adds paper ids: `paper` (title, matches his paper counts) and `paper_ja`
#      (journal + article, the cluster id in the EKW code)
#    * keeps the stored t in `t_stored` (use abs() of it for tests that are not his)
#
#  deround_kudrin(df)  one draw of the EKW de-rounding on that sample (set the seed first)
#    * de-rounds whatever the paper reported, by `report`: t-stat (t), p-value (p),
#      coefficient and SE (s), confidence interval (ci; SE kept fixed)
#    * uniform noise of one unit in the last decimal, with decimals read off the
#      stored values, as in his code
#    * fallback = "kudrin" (default): an input the report type needs but is missing
#      is set to 999 and a de-rounded p above 1 is set to 1, as in his code;
#      fallback = "t": such rows de-round the reported t instead
#
#  Both return the data with t = |t|. His Table 2 uses no de-rounding, and de-rounded
#  draws can only match his in distribution, since his draws depend on everything
#  else his script does with the random number stream.
################################################################################

prep_kudrin <- function(file) {
  num  <- function(x) suppressWarnings(as.numeric(x))
  df   <- haven::read_dta(file)
  need <- c("t", "mu_orig", "sd_orig", "pv", "report", "method", "title", "journal", "article")
  if (length(setdiff(need, names(df))) > 0)
    stop("missing columns: ", paste(setdiff(need, names(df)), collapse = ", "))
  
  mu <- num(df$mu_orig)
  df <- df[!(is.na(mu) | mu == 0), ]                          # "observations with missing information"
  
  df$method   <- as.character(df$method)
  df$title    <- as.character(df$title)
  df$report   <- as.character(df$report)
  df$t_stored <- num(df$t)
  df$p_raw    <- 2 * (1 - pnorm(abs(df$t_stored)))            # his raw p-value
  df$t        <- ifelse(df$p_raw == 0, 1e6, abs(df$t_stored))  # same p in t form
  df$paper    <- as.integer(factor(df$title))
  df$paper_ja <- as.integer(factor(paste0(df$journal, df$article)))
  
  m   <- c("DID", "IV", "RCT", "RDD")
  chk <- rbind(tests  = as.vector(table(df$method)[m]),
               papers = as.vector(tapply(df$title, df$method, function(x) length(unique(x)))[m]))
  colnames(chk) <- m
  print(chk)                                                  # want 5780 5158 7101 3117 / 241 281 141 85
  if (!all(chk["tests", ] == c(5780, 5158, 7101, 3117)) || !all(chk["papers", ] == c(241, 281, 141, 85)))
    warning("sample does not match Kudrin's Table 2 counts")
  if (anyNA(df$t)) warning(sum(is.na(df$t)), " rows have no usable t")
  df
}


deround_kudrin <- function(df, fallback = c("kudrin", "t")) {
  fallback <- match.arg(fallback)
  num <- function(x) suppressWarnings(as.numeric(x))
  dec <- function(x) ifelse(grepl("\\.", x), nchar(sub(".*\\.", "", x)), 0)   # his decimal_places()
  N   <- nrow(df)
  rt  <- df$report
  t0  <- df$t_stored; pv0 <- num(df$pv); mu0 <- num(df$mu_orig); se0 <- num(df$sd_orig)
  
  # Noise in his order (t, p-value, coefficient, SE); missing inputs become 999 as in his loop
  t_dr  <- abs(t0 + (runif(N) - 0.5) * 10^(-dec(t0)))
  pv_dr <- abs(ifelse(is.na(pv0), 999, pv0) + (runif(N) - 0.5) * 10^(-dec(pv0)))
  mu_dr <-     ifelse(is.na(mu0), 999, mu0) + (runif(N) - 0.5) * 10^(-dec(mu0))
  se_dr <- abs(ifelse(is.na(se0), 999, se0) + (runif(N) - 0.5) * 10^(-dec(se0)))
  se_ci <-     ifelse(is.na(se0), 999, se0)
  
  t_new <- ifelse(rt %in% "t",  t_dr,
                  ifelse(rt %in% "p",  qnorm(pmin(pv_dr, 1) / 2, lower.tail = FALSE),   # p > 1 -> 1, i.e. t = 0
                         ifelse(rt %in% "s",  abs(mu_dr / se_dr),
                                ifelse(rt %in% "ci", abs(mu_dr / se_ci), t_dr))))
  
  needs_999 <- (rt %in% "p" & is.na(pv0)) | (rt %in% c("s", "ci") & is.na(se0))
  other     <- !(rt %in% c("t", "p", "s", "ci"))
  if (fallback == "t") t_new[needs_999] <- t_dr[needs_999]
  if (any(needs_999)) message(sum(needs_999), " rows lack the input their report type needs (",
                              if (fallback == "kudrin") "set to 999 as in his code" else "de-rounded from t", ")")
  if (any(other)) message(sum(other), " rows have no usable report type; de-rounded from t")
  
  t_new[is.infinite(t_new)] <- 1e6                            # p = 0, e.g. a CI row with SE 0
  if (anyNA(t_new)) warning(sum(is.na(t_new)), " de-rounded t are NA")
  df$t <- t_new
  df$p <- 2 * pnorm(-t_new)
  df
}