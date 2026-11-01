#### =====================================================================================================
#### ekw_safe.R: the six EKW (2022) tests of Tests.R, guarded so that a test that fails never enters an average
#### as if it were a p-value. Source it after Tests.R (application_cochrane.R and application_MM_analysis.R do).
#### Tests.R's CoxShi() returns the code 999 when its optimizer does not converge and 888 when its covariance
#### matrix is singular, and any of the six tests can stop with an error. ekw_guard() turns each of these into
#### failed_p (NA by default) and records why in the attribute "status":
####   0    ran
####   1    error (the run now continues instead of stopping)
####   2    not a p-value: missing, outside [0, 1], or not a single number
####   888  Cox-Shi: singular covariance matrix
####   999  Cox-Shi: optimizer did not converge
#### Sourcing this file runs no test. ekw_safe_selfcheck() is a one-time check to run by hand.
#### =====================================================================================================

EKW_STATUS <- c("0" = "ran", "1" = "error", "2" = "not a p-value",
                "888" = "Cox-Shi: singular covariance matrix", "999" = "Cox-Shi: optimizer did not converge")

# One test. expr is the test's call; it is evaluated here, inside tryCatch(), so an error is caught. Returns the
# p-value, or failed_p if the test failed, with its status code (above) in the attribute "status".
ekw_guard <- function(expr, failed_p = NA_real_, label = "EKW test") {
  msg <- NULL
  v   <- tryCatch(as.numeric(expr), error = function(e) { msg <<- conditionMessage(e); NULL })
  code <- if (is.null(v)) 1L else if (length(v) != 1L || is.na(v)) 2L else if (v == 999) 999L else
    if (v == 888) 888L else if (v < 0 || v > 1) 2L else 0L
  if (code != 0L) {
    message(sprintf("%s failed (%s%s); its p-value is set to %s", label, EKW_STATUS[[as.character(code)]],
                    if (is.null(msg)) "" else paste0(": ", msg), format(failed_p)))
    v <- failed_p
  }
  structure(as.numeric(v), status = code)
}

# The six tests exactly as the two applications call them, on two-sided p-values ps (ids: the clusters of the
# Cox-Shi tests). Returns list(p = the six p-values, status = their codes), both named
# LCM, Fisher, Disc, Binomial, CS1, CS2B.
ekw_six <- function(ps, ids = seq_along(ps), pmax = 0.15, failed_p = NA_real_) {
  pmin <- min(ps[ps > 0])
  g    <- function(label, expr) ekw_guard(expr, failed_p, paste("EKW", label))
  out  <- list(LCM      = g("LCM",      LCM(ps, pmin, pmax)),
               Fisher   = g("Fisher",   Fisher(ps, pmin, pmax)),
               Disc     = g("Disc",     Discontinuity_test(ps, .05)),
               Binomial = g("Binomial", Binomial(ps, .04, .05, "c")),
               CS1      = g("CS1",      CoxShi(ps, ids, pmin, pmax, 15, 1, 0)),
               CS2B     = g("CS2B",     CoxShi(ps, ids, pmin, pmax, 15, 2, 1)))
  list(p      = vapply(out, as.numeric, numeric(1)),
       status = vapply(out, function(x) attr(x, "status"), integer(1)))
}

# Summaries over draws that leave out the draws in which a test failed (NA); NA if it failed in every draw
mean_finite <- function(v) { v <- v[is.finite(v)]; if (length(v)) mean(v) else NA_real_ }
min_finite  <- function(v) { v <- v[is.finite(v)]; if (length(v)) min(v) else NA_real_ }

# One-time check by hand, after sourcing Tests.R and this file: ekw_safe_selfcheck(). It calls set.seed().
ekw_safe_selfcheck <- function() {
  is_p  <- function(x, p, s) isTRUE(all.equal(as.numeric(x), p)) && identical(attr(x, "status"), s)
  is_na <- function(x, s) is.na(as.numeric(x)) && identical(attr(x, "status"), s)
  run1  <- function(expr, failed_p = NA_real_) suppressMessages(ekw_guard(expr, failed_p))
  ok <- c(ran            = is_p(run1(0.25), 0.25, 0L),
          no_convergence = is_na(run1(999), 999L),
          singular_omega = is_na(run1(888), 888L),
          error          = is_na(run1(stop("boom")), 1L),
          missing        = is_na(run1(NA_real_), 2L),
          above_one      = is_na(run1(1.5), 2L),
          as_nonreject   = is_p(run1(999, failed_p = 1), 1, 999L))
  if (all(vapply(c("LCM", "Fisher", "Discontinuity_test", "Binomial", "CoxShi"), exists, TRUE))) {
    set.seed(1)
    ps <- 2 * pnorm(-abs(rnorm(3000) + rexp(3000)))
    e  <- suppressMessages(ekw_six(ps, ids = rep(seq_len(600), each = 5)))
    ok <- c(ok, six_tests = identical(names(e$p), c("LCM", "Fisher", "Disc", "Binomial", "CS1", "CS2B")) &&
              identical(names(e$status), names(e$p)) && all(is.na(e$p) == (e$status != 0L)))
    message("EKW p-values on a simulated sample: ", paste(sprintf("%s %.3f", names(e$p), e$p), collapse = ", "))
  } else {
    message("Tests.R is not loaded, so only ekw_guard() was checked")
  }
  if (!all(ok)) stop("ekw_safe_selfcheck failed: ", paste(names(ok)[!ok], collapse = ", "))
  message("ekw_safe_selfcheck: all ", length(ok), " checks passed")
  invisible(ok)
}