#!/usr/bin/env Rscript
# compute_C.R  --  run with:  Rscript compute_C.R
#
# Computes the constant C in the rounding bound  d_w(theta_0) <= C * E[kappa^6]
# for the projection test, with and without symmetrization.
#
#   C         = (1/30240) * sum_{i=0}^{6} choose(6, i) * Gamma_i * ||phi^{(6-i)}||_2
#   Gamma_i   = sup_{||beta|| = 1} || sum_j beta_j w_j xi_j^{(i)} ||_2
#             = sqrt(largest eigenvalue of G_i),
#   G_i[j, k] = w_j w_k * integral xi_j^{(i)}(t) xi_k^{(i)}(t) dt,   j, k = 0, ..., J
#   w_j       = min(1, 2^(J0 - j))
#
# Coefficient functions, with psi_j = He_j / sqrt(j!) (normalized probabilists' Hermite):
#   unsymmetrized:              xi_j(t) = psi_j(t) phi(t)
#   symmetrized, shifted by c:  xi_j(t) = [psi_j(t - c) phi(t - c) + psi_j(-t - c) phi(-t - c)] / 2
# Shifting alone does not change C, because L2 norms are shift-invariant, so
# "unsymmetrized" covers both the shifted and the unshifted unsymmetrized test.
#
# Base R only. Expected output with the default settings (J = 60):
#   J0 = 20: symmetrized C = 0.1253, unsymmetrized C = 0.1688
#   J0 = 30: symmetrized C = 0.3060, unsymmetrized C = 0.4298
# Other settings: call compute_C(J0, symmetrized, shift, J, t_max, dt) directly.

# g_j(x) = psi_j(x) * phi(x) for j = 0, ..., n - 1 (stored in row j + 1), via the recurrence
#   g_{j+1}(x) = (x g_j(x) - sqrt(j) g_{j-1}(x)) / sqrt(j + 1)
hermite_fns <- function(x, n) {
  G <- matrix(0, nrow = n, ncol = length(x))
  G[1, ] <- dnorm(x)
  if (n > 1) G[2, ] <- x * G[1, ]
  if (n > 2) {
    for (j in 1:(n - 2)) {
      G[j + 2, ] <- (x * G[j + 1, ] - sqrt(j) * G[j, ]) / sqrt(j + 1)
    }
  }
  G
}

# Rows j = 0, ..., J of the i-th derivatives, using
#   d^i/dx^i g_j(x) = (-1)^i * sqrt((j + 1)(j + 2)...(j + i)) * g_{j+i}(x)
deriv_rows <- function(G, i, J) {
  j <- 0:J
  fac <- if (i == 0) rep(1, J + 1) else sapply(j, function(jj) sqrt(prod((jj + 1):(jj + i))))
  (-1)^i * fac * G[j + i + 1, , drop = FALSE]
}

# Gamma_i from the functions xi_j^{(i)}(t) evaluated on the grid (one row per j)
gamma_from_rows <- function(A, w, dt) {
  M <- w * A                          # scale row j by w_j
  G <- tcrossprod(M) * dt             # G[j, k] = w_j w_k * integral, as a Riemann sum
  sqrt(max(eigen(G, symmetric = TRUE, only.values = TRUE)$values))
}

# ||phi^{(m)}||_2 in closed form: sqrt((2m - 1)!! / (2^(m + 1) sqrt(pi))), with (-1)!! = 1
dfact <- function(n) if (n <= 0) 1 else prod(seq(n, 1, by = -2))
phi_norm <- function(m) sqrt(dfact(2 * m - 1) / (2^(m + 1) * sqrt(pi)))

compute_C <- function(J0, symmetrized = TRUE, shift = 1.96, J = 60,
                      t_max = 40, dt = 0.0025) {
  tt <- seq(-t_max, t_max, by = dt)
  n  <- J + 7                         # derivatives up to order 6 need g_j up to j = J + 6
  w  <- pmin(1, 2^(J0 - (0:J)))
  if (symmetrized) {
    G_plus  <- hermite_fns(tt - shift, n)    # g_j(t - c)
    G_minus <- hermite_fns(-tt - shift, n)   # g_j(-t - c)
  } else {
    G_0 <- hermite_fns(tt, n)
  }
  Gamma <- sapply(0:6, function(i) {
    if (symmetrized) {
      # d^i/dt^i of g_j(-t - c) is (-1)^i g_j^{(i)}(-t - c)
      A <- 0.5 * (deriv_rows(G_plus, i, J) + (-1)^i * deriv_rows(G_minus, i, J))
    } else {
      A <- deriv_rows(G_0, i, J)
    }
    gamma_from_rows(A, w, dt)
  })
  names(Gamma) <- paste0("Gamma_", 0:6)
  C <- sum(choose(6, 0:6) * Gamma * sapply(6 - (0:6), phi_norm)) / 30240
  list(C = C, Gamma = Gamma)
}

# Sanity check: the recurrence and the grid reproduce the closed-form ||phi^{(m)}||_2,
# using phi^{(m)} = (-1)^m sqrt(m!) g_m
tt  <- seq(-40, 40, by = 0.0025)
G_c <- hermite_fns(tt, 7)
num <- sapply(0:6, function(m) sqrt(factorial(m) * sum(G_c[m + 1, ]^2) * 0.0025))
stopifnot(max(abs(num - sapply(0:6, phi_norm))) < 1e-8)

res <- NULL
for (J0 in c(20, 30)) {
  for (sym in c(TRUE, FALSE)) {
    out <- compute_C(J0, symmetrized = sym)
    res <- rbind(res, data.frame(
      J0   = J0,
      case = if (sym) "symmetrized + shifted" else "unsymmetrized",
      C    = round(out$C, 4),
      t(signif(out$Gamma, 4))
    ))
  }
}
print(res, row.names = FALSE)