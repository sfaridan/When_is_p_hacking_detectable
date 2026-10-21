################################################################################
#  kudrin_dgp.R
#
#  Covariate-selection p-hacking DGPs of Kudrin (2025, "Testing for and
#  Evaluating the Extent of Selective Reporting", Sec. 5.1), which adopts the
#  design of Elliott, Kudrin & Wuthrich (2025, "The Power of Tests for Detecting
#  p-Hacking", Secs. 5.1.1 and 5.3.1).
#
#  Source this file AFTER functions_all_weighted.R. It plugs into run_sims()
#  through two short hooks (marked "Kudrin" in the patched run_sims()):
#     pi0_shape = "kudrin_gamma"  -> h ~ Gamma(shape 1.4134, rate 1.2980)
#     hack_type = "kudrin_min"    -> minimum approach
#     hack_type = "kudrin_thresh" -> threshold approach (general-to-specific)
#  prob_hack plays the role of Kudrin's tau (fraction of p-hacking studies).
#  Point masses h = 1, 2 use the existing pi0_shape = "point" with h_center.
#
#  DGP for study s (one reported t-statistic per study):
#     X_i, u_i, e_ki iid N(0,1), i = 1..N (N = 200);  gamma_k ~ U[-0.8, 0.8]
#     W_ki = gamma_k X_i + sqrt(1 - gamma_k^2) e_ki,   k = 1..K
#     Y_i  = X_i h_s / sqrt(N) + u_i
#  Researchers regress Y on X and each of the 2^K subsets of W (OLS, no
#  constant, homoskedastic SEs, Student-t p-values as printed by lm()) and
#  test beta = 0 against a two-sided alternative.
#     honest    : the first specification, with all K controls. This is g^np in
#                 EKW Sec. 5.3.1 ("the first p-value that the researchers
#                 obtain") and the single G0 shared by both searches in
#                 Kudrin's Fig. 7.
#     minimum   : the smallest p-value over all 2^K specifications.
#     threshold : all K controls; if p > 0.05, the best of all (K-1)-subsets;
#                 ... stop at the first subset size whose best p is <= 0.05;
#                 if nothing is significant, the smallest p-value overall.
#  The meta-analyst observes the t-statistic of the reported specification.
#
#  Optional parms columns [defaults]:
#     kudrin_K [3], kudrin_N [200], kudrin_gamma_max [0.8], kudrin_alpha [0.05],
#     kudrin_dgp ["regression" = Kudrin's finite-sample regressions;
#                 "asymptotic" = the N -> infinity limit experiment, in which
#                 T | h, gamma is exactly normal (no Student-t distortion)].
################################################################################

KUDRIN_GAMMA <- c(shape = 1.4134, rate = 1.2980)   # EKW (2025) MLE on Brodeur et al. (2020) RCTs


# 2^K x K matrix of 0/1 subset indicators. Row c holds the controls whose bit is set
# in c - 1: row 1 = no controls, row 2^K = all K controls.
.kd_masks <- function(K) {
  matrix(sapply(seq_len(K), function(k) ((0:(2^K - 1)) %/% 2^(k - 1)) %% 2),
         nrow = 2^K, ncol = K)
}


# Raw data for n studies: list of n x N matrices (row s = study s) X, W_1, ..., W_K, Y.
.kd_simulate_data <- function(h, K, N, gamma_max) {
  n   <- length(h)
  gam <- matrix(runif(n * K, -gamma_max, gamma_max), n, K)
  X   <- matrix(rnorm(n * N), n, N)
  W   <- lapply(seq_len(K), function(k)
    X * gam[, k] + matrix(rnorm(n * N), n, N) * sqrt(1 - gam[, k]^2))
  Y   <- X * (h / sqrt(N)) + matrix(rnorm(n * N), n, N)      # beta_s = h_s / sqrt(N)
  c(list(X), W, list(Y))
}


# Residual cross-products of (X, Y) after partialling out every subset of W.
# Depth-first walk over the subset lattice: each subset is reached from its parent
# by one Gaussian-elimination (sweep) step, vectorised across studies.
# A: n x q x q array of cross-products; labs: 0 = X, 1..K = W_k, K + 1 = Y.
.kd_sweep_tree <- function(A, labs, S, K, env) {
  ix  <- match(0L, labs)
  iy  <- match(K + 1L, labs)
  col <- 1L + sum(2L^(S - 1L))
  env$Sxx[, col] <- A[, ix, ix]
  env$Sxy[, col] <- A[, ix, iy]
  env$Syy[, col] <- A[, iy, iy]
  for (k in labs[labs >= 1L & labs <= K]) {        # only controls above max(S) remain in labs
    j    <- match(k, labs)
    keep <- which(labs == 0L | labs == K + 1L | (labs > k & labs <= K))
    q    <- length(keep)
    B    <- array(NA_real_, c(dim(A)[1], q, q))
    for (a in seq_len(q)) for (b in a:q) {
      B[, a, b] <- B[, b, a] <-
        A[, keep[a], keep[b]] - A[, keep[a], j] * A[, keep[b], j] / A[, j, j]
    }
    .kd_sweep_tree(B, labs[keep], c(S, k), K, env)
  }
  invisible(NULL)
}


# OLS t-statistics (and p-values) on X for all 2^K specifications.
.kd_tstats_from_data <- function(V, K, N) {
  n <- nrow(V[[1]])
  p <- K + 2L
  C <- array(NA_real_, c(n, p, p))
  for (a in seq_len(p)) for (b in a:p) C[, a, b] <- C[, b, a] <- rowSums(V[[a]] * V[[b]])
  env <- new.env()
  env$Sxx <- env$Sxy <- env$Syy <- matrix(NA_real_, n, 2^K)
  .kd_sweep_tree(C, c(0L, seq_len(K), K + 1L), integer(0), K, env)
  size <- rowSums(.kd_masks(K))
  df   <- N - size - 1                                       # residual df of each specification
  rss  <- env$Syy - env$Sxy^2 / env$Sxx
  tst  <- env$Sxy / sqrt(env$Sxx * sweep(rss, 2, df, "/"))
  pv   <- matrix(2 * pt(-abs(tst), df = rep(df, each = n)), n, 2^K)
  list(t = tst, p = pv, size = size, df = df)
}


# Limit experiment of the same design (N -> infinity):
#   T_S = (h + G_0 - sum_{k in S} d_k G_k) / sqrt(1 + sum_{k in S} d_k^2),
#   d_k = gamma_k / sqrt(1 - gamma_k^2),  G_0, ..., G_K iid N(0, 1).
.kd_tstats_asymptotic <- function(h, K, gamma_max) {
  n   <- length(h)
  M   <- t(.kd_masks(K))                                     # K x 2^K
  gam <- matrix(runif(n * K, -gamma_max, gamma_max), n, K)
  d   <- gam / sqrt(1 - gam^2)
  G   <- matrix(rnorm(n * K), n, K)
  tst <- ((h + rnorm(n)) - (d * G) %*% M) / sqrt(1 + d^2 %*% M)
  list(t = tst, p = 2 * pnorm(-abs(tst)), size = colSums(M), df = rep(Inf, 2^K))
}


# t-statistics and p-values of all 2^K specifications for studies with local effects h.
kudrin_spec_tstats <- function(h, K = 3, N = 200, gamma_max = 0.8,
                               dgp = c("regression", "asymptotic")) {
  dgp <- match.arg(dgp)
  K   <- as.integer(K)
  if (dgp == "asymptotic") return(.kd_tstats_asymptotic(h, K, gamma_max))
  .kd_tstats_from_data(.kd_simulate_data(h, K, N, gamma_max), K, N)
}


# Column (specification) each study reports under a search rule.
kudrin_choose <- function(st, search = c("minimum", "threshold"), alpha = 0.05) {
  search <- match.arg(search)
  n <- nrow(st$p)
  if (search == "minimum") return(max.col(-st$p, ties.method = "first"))
  j <- rep(NA_integer_, n)
  for (s in max(st$size):0) {                                # general-to-specific
    cols   <- which(st$size == s)
    jbest  <- cols[max.col(-st$p[, cols, drop = FALSE], ties.method = "first")]
    hit    <- is.na(j) & st$p[cbind(seq_len(n), jbest)] <= alpha
    j[hit] <- jbest[hit]
  }
  miss <- is.na(j)                                           # nothing significant: global minimum
  if (any(miss)) j[miss] <- max.col(-st$p[miss, , drop = FALSE], ties.method = "first")
  j
}


# Honest and p-hacked t-statistics (n x 3 matrix) for local effects h, in chunks.
kudrin_draw <- function(h, K = 3, N = 200, gamma_max = 0.8, alpha = 0.05,
                        dgp = "regression", chunk = 5000) {
  out <- matrix(NA_real_, length(h), 3,
                dimnames = list(NULL, c("honest", "threshold", "minimum")))
  for (idx in split(seq_along(h), ceiling(seq_along(h) / chunk))) {
    st   <- kudrin_spec_tstats(h[idx], K, N, gamma_max, dgp)
    rows <- seq_along(idx)
    out[idx, "honest"]    <- st$t[, ncol(st$t)]              # all K controls
    out[idx, "threshold"] <- st$t[cbind(rows, kudrin_choose(st, "threshold", alpha))]
    out[idx, "minimum"]   <- st$t[cbind(rows, kudrin_choose(st, "minimum", alpha))]
  }
  out
}


# Hook called from run_sims(): the reported t-statistics of one meta-sample.
#   hs        : local effects h_s of the n studies (built from pi0_shape / h_center)
#   hacker    : logical, TRUE for the studies that p-hack (rands <= prob_hack)
#   parms_row : the current row of parms (supplies the optional kudrin_* columns)
kudrin_tstats <- function(hs, hacker, hack_type, parms_row = NULL) {
  arg <- function(name, default) {
    v <- if (!is.null(parms_row) && name %in% names(parms_row)) parms_row[[name]] else NULL
    if (is.null(v) || length(v) == 0 || is.na(v[1])) return(default)
    if (is.factor(v)) as.character(v[1]) else v[1]
  }
  search <- switch(as.character(hack_type),
                   kudrin_min    = "minimum",
                   kudrin_thresh = "threshold",
                   stop("Unknown Kudrin hack_type: ", hack_type))
  d <- kudrin_draw(hs, K = arg("kudrin_K", 3), N = arg("kudrin_N", 200),
                   gamma_max = arg("kudrin_gamma_max", 0.8),
                   alpha = arg("kudrin_alpha", 0.05), dgp = arg("kudrin_dgp", "regression"))
  ifelse(hacker, d[, search], d[, "honest"])
}


# ---------------------------------------------------------------------------- #
#  Checks (run once; neither is needed by run_sims)
# ---------------------------------------------------------------------------- #

# Pre-flight check. (1) The vectorised sweep t-statistics equal lm() on the same simulated
# data for every specification. (2) The run_sims() hook runs for both searches and both
# DGPs, with factor-valued parms columns as produced by expand.grid().
kudrin_selfcheck <- function(K = 3, N = 200, nstudies = 4, seed = 1) {
  set.seed(seed)
  K  <- as.integer(K)
  h  <- rgamma(nstudies, shape = KUDRIN_GAMMA[["shape"]], rate = KUDRIN_GAMMA[["rate"]])
  V  <- .kd_simulate_data(h, K, N, 0.8)
  st <- .kd_tstats_from_data(V, K, N)
  masks <- .kd_masks(K)
  err   <- 0
  for (s in seq_len(nstudies)) for (cc in seq_len(2^K)) {
    Z <- matrix(V[[1]][s, ], ncol = 1)                       # column 1 = X
    for (k in which(masks[cc, ] == 1)) Z <- cbind(Z, V[[k + 1]][s, ])
    y   <- V[[K + 2]][s, ]
    fit <- lm(y ~ 0 + Z)
    stopifnot(fit$df.residual == st$df[cc])
    # by position: no constant, so row 1 is X (lm's names for matrix terms depend on ncol)
    err <- max(err, abs(coef(summary(fit))[1, "t value"] - st$t[s, cc]))
  }
  if (err > 1e-8) stop(sprintf("Sweep t-statistics differ from lm() by %.2e.", err))
  for (ht in c("kudrin_min", "kudrin_thresh")) for (dg in c("regression", "asymptotic")) {
    row <- data.frame(kudrin_K = K, kudrin_dgp = dg, stringsAsFactors = TRUE)
    tt  <- kudrin_tstats(rep(2, 50), hacker = rep(c(TRUE, FALSE), 25),
                         hack_type = factor(ht), parms_row = row)
    stopifnot(length(tt) == 50, all(is.finite(tt)))
  }
  message(sprintf(paste("kudrin_selfcheck(K = %d): max |t(sweep) - t(lm)| over %d regressions",
                        "= %.1e; hook OK for both searches and both DGPs."),
                  K, nstudies * 2^K, err))
  invisible(err)
}


# Reproduces Kudrin (2025) Fig. 7: CDFs of honest (G0) and p-hacked (G1) p-values,
# p = 2 Phi(-|t|). Takes about a minute with the default nstudies.
kudrin_fig7_check <- function(nstudies = 1e5, K_vals = c(3, 5, 7), seed = 1) {
  set.seed(seed)
  hd <- list("h = 1"     = function(m) rep(1, m),
             "h = 2"     = function(m) rep(2, m),
             "h ~ Gamma" = function(m) rgamma(m, shape = KUDRIN_GAMMA[["shape"]],
                                              rate = KUDRIN_GAMMA[["rate"]]))
  pg  <- seq(0, 1, by = 0.005)
  col <- c("steelblue", "red", "orange")
  op  <- par(mfrow = c(length(K_vals), length(hd)), mar = c(4, 4, 2, 1))
  on.exit(par(op))
  for (K in K_vals) for (nm in names(hd)) {
    pv   <- kudrin_draw(hd[[nm]](nstudies), K = K)
    pv[] <- 2 * pnorm(-abs(pv))
    cdf  <- sapply(colnames(pv), function(cc) ecdf(pv[, cc])(pg))
    matplot(pg, cdf, type = "l", lty = c(1, 3, 2), lwd = 2, col = col,
            xlab = "p-value", ylab = "CDF", main = paste0(nm, ", K = ", K))
    legend("bottomright", c("G0", "G1 (thresholding)", "G1 (minimum)"),
           lty = c(1, 3, 2), lwd = 2, col = col, bty = "n", cex = 0.8)
  }
}