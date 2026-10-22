################################################################################
#  smooth_pub_dgp.R
#
#  Smooth publication-bias DGP (the counterexample to Kudrin's tests):
#    H ~ N(0, tau2),  T* = H + Z,  Z ~ t(nu)  (nu = 99999 or Inf: standard normal),
#    each study is published with probability  w(t) = 1 - theta * exp(-t^2 / (2 * sw2)).
#  Published density (normal Z):
#    f(t) = [phi_s2(t) - rho * phi_v(t)] / (1 - rho),   s2 = 1 + tau2,
#    v = 1 / (1/s2 + 1/sw2),   rho = theta * sqrt(v / s2),   P(R) = 1 - rho.
#  Defaults tau2 = 2, sw2 = 6, theta = 1: P(R) = 0.18 and f(0) = 0, so H0 fails, yet all
#  three Kudrin (2025) tests have asymptotic power zero.
#  For theta <= smooth_pub_threshold(tau2, sw2) (0.866 at the defaults) the published
#  t-curve is an honest Gaussian mixture: publication bias is present but undetectable,
#  so every valid test should reject at most 5% of the time (a size check).
#
#  Usage: source after functions_all_weighted.R, which needs the two edits marked
#  "smooth_pub" in run_sims(). Rows use pi0_shape = "smooth_pub" and the optional columns
#  pub_tau2, pub_sw2, pub_theta (defaults 2, 6, 1). The DGP replaces the hacking and step
#  truncation of run_sims(), so its rows should have prob_hack = 0 and theta = 1;
#  smooth_pub_rows() and add_smooth_pub_rows() set this up.
################################################################################


# ---------------------------------------------------------------------------- #
#  Closed forms (normal noise)
# ---------------------------------------------------------------------------- #

smooth_pub_consts <- function(tau2 = 2, sw2 = 6, theta = 1) {
  stopifnot(tau2 > 0, sw2 > 0, theta >= 0, theta <= 1)
  s2  <- 1 + tau2                      # variance of the latent T*
  v   <- 1 / (1 / s2 + 1 / sw2)        # variance of the subtracted normal component
  rho <- theta * sqrt(v / s2)
  list(s2 = s2, v = v, rho = rho, PR = 1 - rho)
}

# P(published | T* = t)
smooth_pub_w <- function(t, sw2 = 6, theta = 1) 1 - theta * exp(-t^2 / (2 * sw2))

# Overall publication probability P(R)
smooth_pub_prob <- function(tau2 = 2, sw2 = 6, theta = 1) smooth_pub_consts(tau2, sw2, theta)$PR

# Density and CDF of a published t-statistic
smooth_pub_density <- function(t, tau2 = 2, sw2 = 6, theta = 1) {
  k <- smooth_pub_consts(tau2, sw2, theta)
  (dnorm(t, sd = sqrt(k$s2)) - k$rho * dnorm(t, sd = sqrt(k$v))) / k$PR
}
smooth_pub_cdf <- function(t, tau2 = 2, sw2 = 6, theta = 1) {
  k <- smooth_pub_consts(tau2, sw2, theta)
  (pnorm(t, sd = sqrt(k$s2)) - k$rho * pnorm(t, sd = sqrt(k$v))) / k$PR
}

# The published t-curve is honest (undetectable) iff theta <= this threshold
smooth_pub_threshold <- function(tau2 = 2, sw2 = 6) {
  k <- smooth_pub_consts(tau2, sw2, 1)
  if (k$v < 1) 0 else sqrt(k$s2 * (k$v - 1) / (tau2 * k$v))
}


# ---------------------------------------------------------------------------- #
#  Sampler and run_sims() plumbing
# ---------------------------------------------------------------------------- #

# Exactly n published t-statistics: latent studies are drawn in batches until n are published
r_smooth_pub <- function(n, tau2 = 2, sw2 = 6, theta = 1, nu = Inf) {
  pr  <- smooth_pub_prob(tau2, sw2, theta)
  out <- numeric(0)
  while (length(out) < n) {
    m    <- ceiling(1.1 * (n - length(out)) / pr) + 100
    tl   <- rnorm(m, sd = sqrt(tau2)) + rt(m, df = nu)       # latent t-statistics T* = H + Z
    keep <- runif(m) < smooth_pub_w(tl, sw2, theta)           # publication
    out  <- c(out, tl[keep])
  }
  out[seq_len(n)]
}

# DGP parameters of row `parm` of a run_sims() grid; missing columns or NAs give the defaults
smooth_pub_par <- function(parms, parm) {
  val <- function(nm, default) {
    v <- if (nm %in% names(parms)) parms[[nm]][parm] else NA
    v <- suppressWarnings(as.numeric(as.character(v)))
    if (length(v) != 1 || is.na(v)) default else v
  }
  list(tau2 = val("pub_tau2", 2), sw2 = val("pub_sw2", 6), theta = val("pub_theta", 1))
}

# New rows for run_sims(): copies one row of an existing grid (template) and sets the DGP.
# Vector arguments are crossed. Default pub_theta = c(0.8, 1): size check and counterexample.
smooth_pub_rows <- function(template, n = template$n[1], pub_theta = c(0.8, 1), pub_tau2 = 2,
                            pub_sw2 = 6, kudrin_tests = TRUE, omit_EWK = FALSE, omit_proj = FALSE) {
  g    <- expand.grid(n = n, pub_theta = pub_theta, pub_tau2 = pub_tau2, pub_sw2 = pub_sw2,
                      KEEP.OUT.ATTRS = FALSE)
  rows <- template[rep(1, nrow(g)), , drop = FALSE]
  rows[] <- lapply(rows, function(x) if (is.factor(x)) as.character(x) else x)
  for (nm in names(g)) rows[[nm]] <- g[[nm]]
  rows$pi0_shape      <- "smooth_pub"
  rows$prob_hack      <- 0
  rows$theta          <- 1
  rows$kudrin_tests   <- kudrin_tests
  rows$omit_EWK       <- omit_EWK
  rows$omit_proj      <- omit_proj
  rows$pub_PR         <- mapply(smooth_pub_prob, g$pub_tau2, g$pub_sw2, g$pub_theta)
  rows$pub_detectable <- g$pub_theta > mapply(smooth_pub_threshold, g$pub_tau2, g$pub_sw2)
  rownames(rows) <- NULL
  rows
}

# Appends smooth_pub rows to an existing grid; columns missing on either side are filled with NA
add_smooth_pub_rows <- function(parms, ...) {
  new <- smooth_pub_rows(parms[1, , drop = FALSE], ...)
  parms[] <- lapply(parms, function(x) if (is.factor(x)) as.character(x) else x)
  for (nm in setdiff(names(new), names(parms))) parms[[nm]] <- NA
  for (nm in setdiff(names(parms), names(new))) new[[nm]] <- NA
  out <- rbind(parms, new[, names(parms), drop = FALSE])
  rownames(out) <- NULL
  out
}


# ---------------------------------------------------------------------------- #
#  Pre-flight check (a few seconds)
# ---------------------------------------------------------------------------- #

smooth_pub_selfcheck <- function(n = 2e5, tau2 = 2, sw2 = 6, theta = 1, seed = 1) {
  set.seed(seed)
  x    <- r_smooth_pub(n, tau2, sw2, theta, nu = Inf)
  ks   <- ks.test(x, smooth_pub_cdf, tau2 = tau2, sw2 = sw2, theta = theta)
  tl   <- rnorm(n, sd = sqrt(tau2)) + rnorm(n)
  acc  <- mean(runif(n) < smooth_pub_w(tl, sw2, theta))
  mass <- integrate(smooth_pub_density, -Inf, Inf, tau2 = tau2, sw2 = sw2, theta = theta)$value
  thr  <- smooth_pub_threshold(tau2, sw2)
  PR   <- smooth_pub_prob(tau2, sw2, theta)
  message(sprintf(paste("smooth_pub_selfcheck: P(R) = %.4f (simulated %.4f); density integrates to %.6f;",
                        "f(0) = %.1e; KS p-value %.3f; threshold %.3f, so theta = %.2f is %s."),
                  PR, acc, mass, smooth_pub_density(0, tau2, sw2, theta), ks$p.value, thr, theta,
                  if (theta > thr) "detectable" else "undetectable (honest t-curve)"))
  stopifnot(abs(acc - PR) < 5 * sqrt(PR * (1 - PR) / n), abs(mass - 1) < 1e-6, ks$p.value > 1e-4)
  invisible(list(ks = ks, accept = acc, PR = PR, threshold = thr))
}