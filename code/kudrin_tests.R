################################################################################
#  kudrin_tests.R
#
#  The three tests of Kudrin (2025, "Testing for and Evaluating the Extent of
#  Selective Reporting", Sec. 3.3), sinc kernel, from reported t-statistics:
#    T_inf      KS-type test: sup_p |G_hat(p) - G_tilde(p)| against bootstrap
#               quantiles shifted by the worst-case deconvolution bias (Sec. 3.3.1)
#    T_H        Cox-Shi (2023) CC test of  B_lo <= E[G_hat_J - G_tilde_J] <= B_hi
#               on J = 15 points x_j = 0.01, ..., 0.15 (eq. 11)
#    T_H_hybrid T_H after projecting G_tilde_J onto shape restrictions implied
#               by the null (eqs. 12-13)
#  Source after functions_all_weighted.R. Needs only the nnls package, which
#  run_test_weight_fast() already uses.
#
#  Implementation choices (validated against Kudrin's Fig. 2 in a Python mirror):
#  * G_tilde(p) = mean_s xi_s(p),  xi_s(p) = 1 - [Si(k(c_p - |t_s|)) + Si(k(c_p + |t_s|))] / pi,
#    c_p = Phi^{-1}(1 - p/2), kappa = 3 (eq. 7).
#  * Bias B(p, h) = G_h(p) - E[G_tilde(p) | Pi = delta_h] is computed directly. It is the
#    negative of the closed form printed in Kudrin's App. A.2, whose sign is flipped relative
#    to his definition B = E[G0 - G_tilde] (checked by numerical integration).
#  * T_inf: sup over c = |t| in [0, 8], step 0.01 (p from 1 down to ~1e-15); nonparametric
#    (cluster) bootstrap with `boots` draws, in place of the warp-speed bootstrap Kudrin
#    uses for his Monte Carlo.
#  * Hybrid: the metric Omega_tilde^{-1} of eq. (12) is numerically singular (G_tilde is
#    band-limited, so its 15 values are nearly collinear). The projection uses Omega_tilde
#    with eigenvalues floored at eig_floor * max, i.e. it moves G_tilde only along its
#    well-determined directions; this reproduces Kudrin's published hybrid power in Fig. 2,
#    whereas identity or Omega_hat metrics do not. Restrictions: bin probabilities and
#    their first and second differences (EKW 2022 CS2B type), each bounded by its min and
#    max over point masses h in [0, 40].
#  * Quadratic programs are solved exactly as least-distance problems through NNLS
#    (Lawson & Hanson 1974, Ch. 23); positive NNLS multipliers give the Cox-Shi active set.
################################################################################


# ---------------------------------------------------------------------------- #
#  Numerical building blocks
# ---------------------------------------------------------------------------- #

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
.kt_gauss_legendre <- function(m) {
  j <- seq_len(m - 1)
  J <- matrix(0, m, m)
  J[cbind(j, j + 1)] <- J[cbind(j + 1, j)] <- j / sqrt(4 * j^2 - 1)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}

# Sine integral: cubic spline through a table on [0, 60] (cumulative 8-point
# Gauss-Legendre, step 0.01) and the asymptotic expansion beyond; |error| < 1e-11.
.kt_SI_MAX <- 60
.kt_si_table <- local({
  h   <- 0.01
  x   <- seq(0, .kt_SI_MAX, by = h)
  gl  <- .kt_gauss_legendre(8)
  y   <- outer(x[-1] - h / 2, gl$x * h / 2, "+")
  inc <- as.vector((sin(y) / y) %*% (gl$w * h / 2))
  splinefun(x, c(0, cumsum(inc)), method = "fmm")
})
.kt_Si <- function(x) {
  ax  <- abs(x)
  out <- .kt_si_table(pmin(ax, .kt_SI_MAX))
  big <- which(ax > .kt_SI_MAX)
  if (length(big)) {
    z  <- ax[big]
    z2 <- 1 / z^2
    f  <- (1 - z2 * (2 - z2 * (24 - z2 * (720 - z2 * 40320)))) / z
    g  <- (1 - z2 * (6 - z2 * (120 - z2 * (5040 - z2 * 362880)))) * z2
    out[big] <- pi / 2 - f * cos(z) - g * sin(z)
  }
  sign(x) * out
}

# xi_s(p) for c = cv(p) in cvals (columns) and |t_s| in a (rows).
.kt_xi <- function(cvals, a, kappa) {
  A  <- matrix(a, length(a), length(cvals))
  Cm <- matrix(cvals, length(a), length(cvals), byrow = TRUE)
  1 - (.kt_Si(kappa * (Cm - A)) + .kt_Si(kappa * (Cm + A))) / pi
}

# I(a) = int_0^kappa sin(w a) exp(-w^2 / 2) / w dw (odd in a), spline table on |a| <= 50.
.kt_make_I <- function(kappa, amax = 50, step = 0.005) {
  gl <- .kt_gauss_legendre(256)
  w  <- (gl$x + 1) * kappa / 2
  wt <- gl$w * kappa / 2 * exp(-w^2 / 2) / w
  a  <- seq(0, amax, by = step)
  sf <- splinefun(a, as.vector(sin(outer(a, w)) %*% wt), method = "fmm")
  function(x) {
    if (any(abs(x) > amax)) stop("I(a) requested outside its table")
    sign(x) * sf(abs(x))
  }
}

# Point masses: G_h(p) and m(p, h) = E[G_tilde(p) | Pi = delta_h], for c = cv(p) in cvals (rows)
# and h in hgrid (columns). The bias is B = G - m.
.kt_point_mass <- function(cvals, hgrid, Ifun) {
  C <- matrix(cvals, length(cvals), length(hgrid))
  H <- matrix(hgrid, length(cvals), length(hgrid), byrow = TRUE)
  list(G = pnorm(H - C) + pnorm(-H - C),
       m = 1 - (Ifun(C - H) + Ifun(C + H)) / pi)
}

# Least-distance programming: min ||x|| s.t. G x >= h (Lawson & Hanson 1974, Alg. 23.27).
# Returns x and the constraints with positive NNLS multipliers (the active set).
.kt_ldp <- function(G, h) {
  n <- ncol(G)
  s <- max(abs(h))
  if (s == 0) return(list(x = rep(0, n), active = rep(FALSE, nrow(G))))
  E <- rbind(t(G), h / s)
  f <- c(rep(0, n), 1)
  y <- nnls::nnls(E, f)$x
  r <- as.vector(E %*% y) - f
  if (abs(r[n + 1]) < 1e-14) stop("Kudrin tests: incompatible constraints in a quadratic program")
  list(x = -r[seq_len(n)] / r[n + 1] * s, active = y > 0)
}


# ---------------------------------------------------------------------------- #
#  Grids, bias bounds and shape restrictions (computed once per setting, cached)
# ---------------------------------------------------------------------------- #

.kt_cache <- new.env()

.kt_setup <- function(kappa = 3, xJ = seq(0.01, 0.15, by = 0.01), ks_cmax = 8,
                      ks_step = 0.01, orders = 1:3) {
  key <- paste(kappa, paste(xJ, collapse = ","), ks_cmax, ks_step, paste(orders, collapse = ","))
  if (!is.null(.kt_cache[[key]])) return(.kt_cache[[key]])
  Ifun <- .kt_make_I(kappa)
  cg   <- seq(0, ks_cmax, by = ks_step)                      # T_inf grid in c = |t|
  hks  <- seq(0, 40, by = 0.01)
  bbar <- 0                                                  # sup-norm worst-case bias
  for (idx in split(seq_along(cg), ceiling(seq_along(cg) / 100))) {
    pm   <- .kt_point_mass(cg[idx], hks, Ifun)
    bbar <- max(bbar, abs(pm$G - pm$m))
  }
  cJ  <- qnorm(xJ / 2, lower.tail = FALSE)
  pm  <- .kt_point_mass(cJ, seq(0, 40, by = 0.002), Ifun)
  B   <- pm$G - pm$m
  J   <- length(xJ)
  E   <- rbind(0, diag(J))                                   # (G(0) = 0, G(x_1), ..., G(x_J))
  LM  <- do.call(rbind, lapply(orders, function(k) diff(E, differences = k)))
  LMm <- LM %*% pm$m
  st  <- list(kappa = kappa, xJ = xJ, cJ = cJ, cg = cg, bbar = bbar,
              blo = apply(B, 1, min), bhi = apply(B, 1, max),
              LM = LM, lo = apply(LMm, 1, min), hi = apply(LMm, 1, max))
  assign(key, st, envir = .kt_cache)
  st
}


# ---------------------------------------------------------------------------- #
#  The tests
# ---------------------------------------------------------------------------- #

# tstats   : reported t-statistics (signs are irrelevant)
# cluster  : optional article identifiers (cluster-robust covariances and cluster bootstrap)
# boots    : bootstrap draws for T_inf
# Returns p-values, rejections at level alpha, statistics, Cox-Shi degrees of freedom,
# the sup-norm bias bound and whether the hybrid restriction was active.
kudrin_tests <- function(tstats, cluster = NULL, kappa = 3, boots = 199, alpha = 0.05,
                         xJ = seq(0.01, 0.15, by = 0.01), ks_cmax = 8, ks_step = 0.01,
                         orders = 1:3, eig_floor = 1e-8) {
  st <- .kt_setup(kappa, xJ, ks_cmax, ks_step, orders)
  a  <- abs(as.numeric(tstats))
  S  <- length(a)
  if (is.null(cluster)) cluster <- seq_len(S)
  cluster <- as.integer(factor(cluster))
  ng <- as.vector(rowsum(rep(1, S), cluster))
  
  # T_inf
  P <- length(st$cg)
  Z <- matrix(0, S, P)                                       # 1{|t_s| >= c} - xi_s(c)
  for (idx in split(seq_len(P), ceiling(seq_len(P) / 100))) {
    Z[, idx] <- outer(a, st$cg[idx], ">=") - .kt_xi(st$cg[idx], a, kappa)
  }
  Zg    <- rowsum(Z, cluster)
  D     <- colSums(Zg) / S
  T_inf <- max(abs(D))
  cnt   <- rmultinom(boots, size = length(ng), prob = rep(1, length(ng)))   # clusters x boots
  Ds    <- crossprod(cnt, Zg) / as.vector(crossprod(cnt, ng))
  Ts    <- apply(abs(sweep(Ds, 2, D)), 1, max)
  p_inf <- mean(T_inf < Ts + st$bbar)
  
  # T_H
  XiJ <- .kt_xi(st$cJ, a, kappa)                             # S x J
  IJ  <- outer(2 * pnorm(-a), st$xJ, "<=") * 1               # 1{p_s <= x_j}
  zJ  <- IJ - XiJ
  Bh  <- colMeans(zJ)
  Cl  <- t(chol(crossprod(rowsum(sweep(zJ, 2, Bh), cluster)) / S))   # Omega_hat = Cl Cl'
  cc_test <- function(b) {                                   # Cox-Shi CC test on the bias box
    s    <- .kt_ldp(rbind(Cl, -Cl), c(st$blo - b, b - st$bhi))
    k    <- sum(s$active)
    stat <- S * sum(s$x^2)
    c(stat = stat, df = k, p = if (k == 0) 1 else pchisq(stat, k, lower.tail = FALSE))
  }
  H <- cc_test(Bh)
  
  # T_H_hybrid
  Gt  <- colMeans(XiJ)
  ev  <- eigen(crossprod(rowsum(sweep(XiJ, 2, Gt), cluster)) / S, symmetric = TRUE)
  Bm  <- sweep(ev$vectors, 2, sqrt(pmax(ev$values, eig_floor * max(ev$values))), "*")
  r0  <- as.vector(st$LM %*% Gt)
  GR  <- Gt
  restricted <- any(r0 < st$lo | r0 > st$hi)
  if (restricted) {
    LB <- st$LM %*% Bm
    u  <- .kt_ldp(rbind(LB, -LB), c(st$lo - r0, r0 - st$hi))$x
    GR <- Gt + as.vector(Bm %*% u)
  }
  Hy <- cc_test(colMeans(IJ) - GR)
  
  pval <- c(T_inf = p_inf, T_H = H[["p"]], T_H_hybrid = Hy[["p"]])
  list(pval = pval, reject = pval <= alpha,
       stat = c(T_inf = T_inf, T_H = H[["stat"]], T_H_hybrid = Hy[["stat"]]),
       df = c(T_H = H[["df"]], T_H_hybrid = Hy[["df"]]),
       bias_bound_inf = st$bbar, restricted = restricted)
}


# ---------------------------------------------------------------------------- #
#  Pre-flight check (run once)
# ---------------------------------------------------------------------------- #

kudrin_tests_selfcheck <- function(kappa = 3, seed = 1) {
  set.seed(seed)
  # 1. Sine integral vs integrate() on unit pieces
  si_ref <- function(x) {
    br <- unique(c(seq(0, abs(x), by = 1), abs(x)))
    sign(x) * sum(mapply(function(lo, hi)
      integrate(function(y) ifelse(y == 0, 1, sin(y) / y), lo, hi, rel.tol = 1e-13)$value,
      br[-length(br)], br[-1]))
  }
  xs     <- c(-75, -7.3, -0.2, 0.05, 1.7, 13.1, 59.99, 60.01, 140)
  err_si <- max(abs(.kt_Si(xs) - sapply(xs, si_ref)))
  # 2. I(a) table vs integrate(); expectation of the estimator under point masses
  Ifun  <- .kt_make_I(kappa)
  as_   <- c(0.3, 2, 9.7, 45)
  err_I <- max(abs(Ifun(as_) - sapply(as_, function(a0)
    integrate(function(w) sin(w * a0) * exp(-w^2 / 2) / w, 0, kappa,
              subdivisions = 1000L, rel.tol = 1e-12)$value)))
  ch    <- rbind(c(1.96, 0), c(1.44, 1), c(2.58, 2), c(0.5, 3))
  err_m <- max(apply(ch, 1, function(v) {
    Exi <- integrate(function(tt) as.vector(.kt_xi(v[1], abs(tt), kappa)) * dnorm(tt - v[2]),
                     v[2] - 12, v[2] + 12, subdivisions = 1000L, rel.tol = 1e-12)$value
    abs(Exi - (1 - (Ifun(v[1] - v[2]) + Ifun(v[1] + v[2])) / pi))
  }))
  # 3. Least-distance QP on a random box problem: feasibility and KKT conditions
  J  <- 15
  A0 <- matrix(rnorm(J^2), J); Om <- crossprod(A0) / J + diag(0.01, J)
  b  <- rnorm(J, 0, 3e-3); lo <- -runif(J, 1e-4, 7e-4); hi <- runif(J, 1e-4, 7e-4)
  Cl <- t(chol(Om))
  v  <- .kt_ldp(rbind(Cl, -Cl), c(lo - b, b - hi))$x
  x  <- b + as.vector(Cl %*% v)
  g  <- as.vector(solve(Om, x - b))                          # gradient of (x-b)' Om^{-1} (x-b) / 2
  tl <- 1e-9 * max(abs(g), 1e-12); tb <- 1e-12
  kkt_ok <- all(x >= lo - tb & x <= hi + tb) &&
    all(ifelse(abs(x - lo) <= tb, g >= -tl, ifelse(abs(x - hi) <= tb, g <= tl, abs(g) <= tl)))
  # 4. Full runs: honest sample, threshold-hacked sample, singleton clusters
  n  <- 3000
  t0 <- rnorm(n) + rgamma(n, shape = 1.4134, rate = 1.2980)
  t1 <- rnorm(n) + 1; t2 <- rnorm(n) + 1; th <- ifelse(t1 > 1.96, t1, pmax(t1, t2))
  k0 <- kudrin_tests(t0, boots = 50, kappa = kappa)
  kc <- kudrin_tests(t0, cluster = sample(paste0("a", seq_len(n))), boots = 50, kappa = kappa)  # singletons, relabelled
  kh <- kudrin_tests(th, boots = 50, kappa = kappa)
  stopifnot(err_si < 1e-8, err_I < 1e-8, err_m < 1e-7, kkt_ok,
            all(is.finite(c(k0$pval, kh$pval))), all(c(k0$pval, kh$pval) >= 0 & c(k0$pval, kh$pval) <= 1),
            isTRUE(all.equal(k0$stat, kc$stat)), isTRUE(all.equal(k0$pval[-1], kc$pval[-1])), names(k0$pval) == c("T_inf", "T_H", "T_H_hybrid"))
  if (!kh$restricted) warning("Hybrid restriction did not bind in the threshold-hacked test sample.")
  message(sprintf(paste("kudrin_tests_selfcheck: Si err %.1e, I err %.1e, E[G_tilde] err %.1e, QP KKT ok;",
                        "sup bias bound %.2e; p-values honest (%s), threshold-hacked (%s)."),
                  err_si, err_I, err_m, k0$bias_bound_inf,
                  paste(sprintf("%.3f", k0$pval), collapse = ", "),
                  paste(sprintf("%.3f", kh$pval), collapse = ", ")))
  invisible(list(honest = k0, hacked = kh))
}