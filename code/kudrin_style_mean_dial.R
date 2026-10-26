################################################################################
#  kudrin_style_figures.R
#
#  New versions of Kudrin (2025) Figs. 2 and 3 with the distributions of true effects (Pi0)
#  and the p-hacking methods of your choice: one figure per sample size S (his Fig. 2 has
#  S = 1,000 and his Fig. 3 S = 5,000), one row per p-hacking method, one panel per Pi0, and
#  tau, the share of results that are p-hacked, on the horizontal axis. Same method, tests and
#  engine (kudrin_figures_engine.R) as reproduce_kudrin_fig2_fig3.R: for each (method family,
#  panel), a pool of studies is simulated once, each with an honest and a p-hacked t-statistic;
#  each meta-sample draws S studies from the pool with replacement and reports, for each, the
#  p-hacked t-statistic with probability tau and the honest one otherwise.
#
#  p-hacking methods (CFG$methods), and the honest t-statistic each is paired with:
#    "kudrin_threshold"  Kudrin's covariate search, thresholding approach   } your kudrin_draw();
#    "kudrin_minimum"    Kudrin's covariate search, minimum approach        } honest = all K
#                        controls. Settings CFG$K, N, gamma_max, dgp, search_alpha.
#    "threshold"         your run_sims() "threshold": t1 = h + e1; if t1 <= cv, report
#                        max(t1, t2)                                        } honest = t1, the
#    "max"               your run_sims() "max": the largest of max_draws     } same draw for both;
#                        t-statistics h + e_k (8 in run_sims())             } e_k ~ Student-t(nu)
#    "smooth_pub"        your smooth publication bias: r_smooth_pub() from smooth_pub_dgp.R,
#                        called as in run_sims(); honest = h + e1. r_smooth_pub() draws its own
#                        effects h ~ N(0, tau2), so every panel must then have pi0 = N(0, tau2),
#                        and tau2 is read from it (the script checks). tau is the share of results
#                        that went through the publication filter.
#  For "threshold", "max" and "smooth_pub" these are the DGPs of run_sims(); what changes is that
#  meta-samples are drawn from a pool, as in Kudrin's design. At tau = 1 a meta-sample is a
#  run_sims() meta-sample with prob_hack = 1.
#
#  Pi0 (CFG$panels): each panel is list(pi0 = ..., id = ..., label = ...); id and label are
#  optional. pi0 is one of
#    - a pi0_shape name from run_sims(): "null", "onehalf", "two", "standard_normal",
#      "normal_11", "hard" (N(2, 0.7^2)), "unif_half", "poisson(2)", "chi2(2)", "fitted"
#      (z-curve fit), "kudrin_gamma" (Gamma(1.4134, RATE 1.2980), as in run_sims()), or
#      "kudrin_pi_hat" (Gamma(1.4134, SCALE 1.2980), which reproduces Kudrin's figures);
#    - a list with a shape and its parameters, in run_sims() conventions:
#        list(shape = "point",   h_center = 1)
#        list(shape = "normal",  h_center = 1, sigma_h = 0.5)
#        list(shape = "uniform", h_center = 1, sigma_h = 2)   # U[h_center - sigma_h/2, h_center + sigma_h/2]
#        list(shape = "chi2",    h_center = 2)                # degrees of freedom = h_center
#        list(shape = "poisson", h_center = 2)
#        list(shape = "double_normal", h_center = 1, sigma_h = 0.5)
#        list(shape = "gamma",   gamma_shape = 2, gamma_scale = 1)   # or gamma_rate = ...
#    - any function(n) returning n draws of h.
#
#  Uses your code only:
#    reproduce_kudrin_simulations_functions.R  kudrin_draw()            Kudrin's DGP
#    kudrin_tests.R                            kudrin_tests()           T_inf, T_H, T_H^hybrid
#    functions_all_weighted.R                  run_test_weight_fast()   projection test
#    smooth_pub_dgp.R                          r_smooth_pub()           only for "smooth_pub"
#    kudrin_figures_engine.R                   meta-samples, test calls, cells, figures
#    Tests.R (EKW)                             only if CFG$run_ekw = TRUE
#  All in CODE_DIR.
#
#  Outputs in OUT_DIR (NAME = CFG$name):
#    NAME_S1000.pdf/.png, NAME_S5000.pdf/.png  the figures, one per S, each written as soon as
#                                              its cells finish
#    results_NAME.csv/.rds                     rejection rates by S, method, panel, tau and test
#  Finished cells are saved in OUT_DIR/cells and pools in OUT_DIR/pools, so an interrupted run
#  resumes where it stopped.
#
#  Run once with QUICK <- TRUE (minutes) to check that everything works, then with FALSE.
################################################################################

rm(list = ls())

## ---- 0. Settings -------------------------------------------------------------------------

ROOT     <- "C:/Users/sfaridani6/Documents/GitHub/When_is_p_hacking_detectable"  #"C:/Users/stefa/OneDrive/Documents/GitHub/When_is_p_hacking_detectable"
CODE_DIR <- file.path(ROOT, "code")

QUICK <- FALSE   # TRUE: small run that checks the pipeline; FALSE: the settings below

CFG <- list(
  name           = "kudrin_style_var_dial",   # names the output folder and files
  S_values       = c(5000),        # one figure per S (his Figs. 2 and 3)
  taus           = seq(0, 1, by = 0.1),  # share of results that are p-hacked
  reps           = 1000,                 # Monte Carlo repetitions per (S, method, panel, tau)
  alpha          = 0.05,                 # test level
  
  # Rows, top to bottom: any of "kudrin_threshold", "kudrin_minimum", "threshold", "max",
  # "smooth_pub"
  methods        = c("max", "threshold"),
  row_titles     = c(max = "Maximum of 8 specifications", threshold = "Threshold: second draw if t <= 1.96"),                 # optional, e.g. c(max = "Best of 8 specifications")
  
  # Panels, left to right: distributions of true effects (see the header)
  panels         = list(
    list(id = "N03", pi0 = list(shape = "normal", h_center = 2, sigma_h = 0.3)),
    list(id = "N05", pi0 = list(shape = "normal", h_center = 2, sigma_h = 0.5)),
    list(id = "N07", pi0 = list(shape = "normal", h_center = 2, sigma_h = 0.7))
       ),

  
  # Other examples:
  #   list(pi0 = list(shape = "point", h_center = 1)), list(pi0 = "kudrin_pi_hat"),
  #   list(pi0 = "poisson(2)", id = "pois2"), list(pi0 = function(n) rexp(n), label = "h ~ Exp(1)")
  
  n_pool         = 1e6,                  # studies per pool (EKW and Kudrin use one million)
  
  # Kudrin's covariate searches, passed to kudrin_draw()
  K              = 3,
  N              = 200,
  gamma_max      = 0.8,
  dgp            = "regression",
  search_alpha   = 0.05,                 # significance threshold of the search
  
  # Your run_sims() methods
  nu             = 9999999,              # degrees of freedom of the Student-t noise (run_sims' nu)
  cv             = 1.96,                 # "threshold": second draw when t1 <= cv
  max_draws      = 8,                    # "max": number of t-statistics (8 in run_sims)
  smooth_pub     = list(sw2 = NULL, theta = NULL),   # "smooth_pub": r_smooth_pub()'s sw2, theta
  
  # Kudrin tests: your kudrin_tests()
  kudrin_boots   = 199,                  # T_inf bootstrap draws
  eig_floor      = 1e-8,                 # hybrid metric floor; extra values add sensitivity columns
  
  # Projection test: your run_test_weight_fast(), with the paper's settings
  run_projection = TRUE,
  num_coeffs     = 41,                   # J = 40
  J0             = 20,                   # weights w_j = min(1, (1/p)^(j - J0))
  p_weight       = 2,
  shift          = 1.96,
  L              = 6.5,
  numgrid        = 3000,
  proj_boots     = 100,
  proj_mode      = "full",               # "full" as in the paper, or "warp" (pooled critical value)
  epsilon_U_pre  = NULL,                 # NULL: compute_epsilons(L, numgrid) once and cache it
  # (slow); or 0.00063, the value in the paper
  
  # EKW tests (Tests.R); without them the figures have no EKW line
  run_ekw        = TRUE,
  ekw_J          = 30,
  
  # Computing
  n_cores        = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
  seed           = 20261023,
  resume         = TRUE,
  run_selfchecks = TRUE
)

OUT_DIR <- file.path(ROOT, "simulations", CFG$name)

if (QUICK) {
  CFG$reps          <- 100
  CFG$n_pool        <- 2e5
  CFG$kudrin_boots  <- 99
  CFG$proj_mode     <- "warp"
  CFG$epsilon_U_pre <- 0.00063
  OUT_DIR           <- paste0(OUT_DIR, "_quick")
}


## ---- 1. Methods -------------------------------------------------------------------------
## family: which pool the method draws from; methods of one family share the pool and its
## honest column (common random numbers). column: the pool column the method reports.

METHODS <- list(
  kudrin_threshold = list(family = "kudrin", column = "kudrin_threshold",
                          title  = sprintf("Covariate search, thresholding (K = %d)", as.integer(CFG$K))),
  kudrin_minimum   = list(family = "kudrin", column = "kudrin_minimum",
                          title  = sprintf("Covariate search, minimum (K = %d)", as.integer(CFG$K))),
  threshold        = list(family = "iid", column = "threshold",
                          title  = sprintf("Threshold: second draw if t <= %s", format(CFG$cv))),
  max              = list(family = "iid", column = "max",
                          title  = sprintf("Maximum of %d t-statistics", as.integer(CFG$max_draws))),
  smooth_pub       = list(family = "iid", column = "smooth_pub",
                          title  = "Smooth publication bias")
)
bad <- setdiff(CFG$methods, names(METHODS))
if (length(bad) > 0) {
  stop("Unknown p-hacking method(s): ", paste(bad, collapse = ", "), ". Choose from ",
       paste(names(METHODS), collapse = ", "), ".")
}
ROWS <- lapply(CFG$methods, function(m) {
  rw <- c(list(id = m), METHODS[[m]])
  if (!is.null(CFG$row_titles) && m %in% names(CFG$row_titles)) rw$title <- CFG$row_titles[[m]]
  rw
})
names(ROWS) <- CFG$methods
NEEDS_KUDRIN <- any(vapply(ROWS, function(rw) rw$family == "kudrin", logical(1)))
NEEDS_SMOOTH <- "smooth_pub" %in% CFG$methods


## ---- 2. Your code ------------------------------------------------------------------------

CODE_FILES  <- file.path(CODE_DIR, c(if (CFG$run_ekw) "Tests.R",
                                     "reproduce_kudrin_simulations_functions.R",
                                     "functions_all_weighted.R",
                                     "kudrin_tests.R",
                                     if (NEEDS_SMOOTH) "smooth_pub_dgp.R"))
ENGINE_FILE <- file.path(CODE_DIR, "kudrin_figures_engine.R")
missing_files <- c(CODE_FILES, ENGINE_FILE)[!file.exists(c(CODE_FILES, ENGINE_FILE))]
if (length(missing_files) > 0) stop("Not found: ", paste(missing_files, collapse = ", "))
for (f in CODE_FILES) source(f)
source(ENGINE_FILE)
stopifnot(exists("kudrin_draw"), exists("kudrin_tests"), exists("run_test_weight_fast"),
          exists("create_basis_stable"), exists("compute_epsilons"), exists("plot_figure"))
if (NEEDS_SMOOTH) stopifnot(exists("r_smooth_pub"))
dir.create(file.path(OUT_DIR, "cells"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "pools"), recursive = TRUE, showWarnings = FALSE)

if (CFG$run_selfchecks) {
  if (NEEDS_KUDRIN) kudrin_selfcheck()    # DGP: vectorised t-statistics equal lm()
  kudrin_tests_selfcheck()                # Kudrin tests: Si, I(a), E[G_tilde], QP, full runs
}


## ---- 3. Distributions of true effects -----------------------------------------------------

PI0_ALIASES <- list(                     # run_sims() pi0_shape names
  null            = list(shape = "point",   h_center = 0),
  onehalf         = list(shape = "point",   h_center = 0.5),
  two             = list(shape = "point",   h_center = 2),
  standard_normal = list(shape = "normal",  h_center = 0, sigma_h = 1),
  normal_11       = list(shape = "normal",  h_center = 1, sigma_h = 1),
  hard            = list(shape = "normal",  h_center = 2, sigma_h = 0.7),
  unif_half       = list(shape = "uniform", h_center = 0, sigma_h = 1),
  "poisson(2)"    = list(shape = "poisson", h_center = 2),
  "chi2(2)"       = list(shape = "chi2",    h_center = 2),
  fitted          = list(shape = "fitted"),
  kudrin_gamma    = list(shape = "gamma",   gamma_shape = 1.4134, gamma_rate = 1.2980),
  kudrin_pi_hat   = list(shape = "gamma",   gamma_shape = 1.4134, gamma_scale = 1.2980,
                         label = quote(italic(h) %~% hat(Pi)))
)

as_pi0 <- function(x) {
  if (is.function(x)) return(list(shape = "custom", fun = x))
  if (is.character(x) && length(x) == 1) {
    if (is.null(PI0_ALIASES[[x]])) {
      stop("Unknown pi0 name: ", x, ". Known names: ", paste(names(PI0_ALIASES), collapse = ", "))
    }
    return(PI0_ALIASES[[x]])
  }
  if (is.list(x) && !is.null(x$shape)) return(x)
  stop("A pi0 must be a run_sims() pi0_shape name, a list with a shape, or a function(n).")
}

# n draws of h, with run_sims()' formulas
draw_pi0 <- function(spec, n) {
  s <- spec$shape
  if (s == "custom") {
    h <- spec$fun(n)
  } else if (s == "point") {
    h <- rep(spec$h_center, n)
  } else if (s == "normal") {
    h <- rnorm(n) * spec$sigma_h + spec$h_center
  } else if (s == "uniform") {
    h <- (runif(n) - 0.5) * spec$sigma_h + spec$h_center
  } else if (s == "chi2") {
    h <- rchisq(n, spec$h_center)
  } else if (s == "poisson") {
    h <- rpois(n, spec$h_center)
  } else if (s == "gamma") {
    if (is.null(spec$gamma_scale)) {
      h <- rgamma(n, shape = spec$gamma_shape, rate = spec$gamma_rate)
    } else {
      h <- rgamma(n, shape = spec$gamma_shape, scale = spec$gamma_scale)
    }
  } else if (s == "double_normal") {
    h  <- rnorm(n) * spec$sigma_h + spec$h_center
    rr <- runif(n)
    h[rr >= 0.5] <- h[rr >= 0.5] - 2
  } else if (s == "fitted") {             # fitted by zcurve, Faridani (2026)
    h <- sample(c(0, 1, 2, 3, 4, 5, 6), size = n, replace = TRUE,
                prob = c(0.02, 0.47, 0.01, 0.27, 0.14, 0.00, 0.09))
  } else {
    stop("Unknown pi0 shape: ", s)
  }
  h <- as.numeric(h)
  if (length(h) != n || any(!is.finite(h))) stop("pi0 '", s, "' did not return ", n, " finite draws.")
  h
}

# Panel title (plotmath)
pi0_label <- function(spec) {
  if (!is.null(spec$label)) return(as.expression(spec$label))
  s  <- spec$shape
  g4 <- function(x) signif(x, 4)
  lab <- if (s == "point") {
    bquote(italic(h) == .(g4(spec$h_center)))
  } else if (s == "normal") {
    bquote(italic(h) %~% italic(N)*"("*.(g4(spec$h_center))*","~.(g4(spec$sigma_h))^2*")")
  } else if (s == "uniform") {
    bquote(italic(h) %~% italic(U)*"["*.(g4(spec$h_center - spec$sigma_h / 2))*","~
             .(g4(spec$h_center + spec$sigma_h / 2))*"]")
  } else if (s == "chi2") {
    bquote(italic(h) %~% chi[.(g4(spec$h_center))]^2)
  } else if (s == "poisson") {
    bquote(italic(h) %~% "Poisson("*.(g4(spec$h_center))*")")
  } else if (s == "gamma" && is.null(spec$gamma_scale)) {
    bquote(italic(h) %~% "Gamma("*.(g4(spec$gamma_shape))*", rate"~.(g4(spec$gamma_rate))*")")
  } else if (s == "gamma") {
    bquote(italic(h) %~% "Gamma("*.(g4(spec$gamma_shape))*", scale"~.(g4(spec$gamma_scale))*")")
  } else if (s == "double_normal") {
    bquote(italic(h) %~% "mixture of"~italic(N)*"("*.(g4(spec$h_center))*","~.(g4(spec$sigma_h))^2*")"~
             "and"~italic(N)*"("*.(g4(spec$h_center - 2))*","~.(g4(spec$sigma_h))^2*")")
  } else if (s == "fitted") {
    quote(italic(h) %~% "z-curve fit")
  } else {
    quote(italic(h) %~% "custom")
  }
  as.expression(lab)
}

PANELS <- lapply(seq_along(CFG$panels), function(k) {
  p    <- CFG$panels[[k]]
  spec <- as_pi0(p$pi0)
  list(id    = if (is.null(p$id)) sprintf("pi%d", k) else p$id,
       spec  = spec,
       title = if (is.null(p$label)) pi0_label(spec) else as.expression(p$label))
})
PANEL_IDS <- vapply(PANELS, function(pn) pn$id, "")
if (anyDuplicated(PANEL_IDS) > 0) stop("Panel ids must be distinct: ", paste(PANEL_IDS, collapse = ", "))
names(PANELS) <- PANEL_IDS

# smooth_pub: r_smooth_pub() draws h ~ N(0, tau2) itself, so tau2 must come from the panel
smooth_pub_settings <- function(spec) {
  if (!identical(spec$shape, "normal") || !isTRUE(all.equal(spec$h_center, 0))) {
    stop("smooth_pub draws its own true effects h ~ N(0, tau2), as pi0_shape = \"smooth_pub\" does ",
         "in run_sims(). Give every panel pi0 = list(shape = \"normal\", h_center = 0, ",
         "sigma_h = sqrt(tau2)), or \"standard_normal\" for tau2 = 1.")
  }
  if (is.null(CFG$smooth_pub$sw2) || is.null(CFG$smooth_pub$theta)) {
    stop("Set CFG$smooth_pub$sw2 and CFG$smooth_pub$theta (the values you pass to r_smooth_pub()).")
  }
  list(tau2 = spec$sigma_h^2, sw2 = CFG$smooth_pub$sw2, theta = CFG$smooth_pub$theta)
}

if (NEEDS_SMOOTH) {
  for (pn in PANELS) {
    sp  <- smooth_pub_settings(pn$spec)          # stops here if a panel does not fit smooth_pub
    msg <- tryCatch(sprintf("smooth_pub, panel %s: tau2 = %g, sw2 = %g, theta = %g, P(R) = %.3f, detectable = %s",
                            pn$id, sp$tau2, sp$sw2, sp$theta, smooth_pub_prob(sp$tau2, sp$sw2, sp$theta),
                            sp$theta > smooth_pub_threshold(sp$tau2, sp$sw2)),
                    error = function(e) paste0("smooth_pub, panel ", pn$id, ": ", conditionMessage(e)))
    message(msg)
  }
}


## ---- 4. Pools of honest and p-hacked t-statistics ----------------------------------------

FAMILY_COLS <- lapply(split(vapply(ROWS, function(rw) rw$column, ""),
                            vapply(ROWS, function(rw) rw$family, "")), unique)
pool_id <- function(family, panel_id) paste(family, panel_id, sep = "__")

# Your run_sims() DGPs for "threshold", "max" and "smooth_pub" (signed t-statistics)
iid_pool <- function(h, spec, cols) {
  n   <- length(h)
  ts1 <- rt(n, CFG$nu) + h
  ts2 <- rt(n, CFG$nu) + h
  out <- cbind(honest = ts1)
  if ("threshold" %in% cols) {            # run_sims: ts1 if ts1 > cv, else max(ts1, ts2)
    out <- cbind(out, threshold = ifelse(ts1 > CFG$cv, ts1, pmax(ts1, ts2)))
  }
  if ("max" %in% cols) {                  # run_sims: pmax(ts1, ..., ts8)
    mx <- pmax(ts1, ts2)
    for (j in seq_len(max(0, CFG$max_draws - 2))) mx <- pmax(mx, rt(n, CFG$nu) + h)
    out <- cbind(out, max = mx)
  }
  if ("smooth_pub" %in% cols) {           # run_sims: r_smooth_pub(n, tau2, sw2, theta, nu)
    sp <- smooth_pub_settings(spec)
    tp <- as.numeric(r_smooth_pub(n, tau2 = sp$tau2, sw2 = sp$sw2, theta = sp$theta, nu = CFG$nu))
    if (length(tp) != n) stop("r_smooth_pub() returned ", length(tp), " t-statistics instead of ", n, ".")
    out <- cbind(out, smooth_pub = tp)
  }
  out
}

# One pool per (family, panel), cached in OUT_DIR/pools. Its seed depends only on the panel's
# position and the family, so adding a method or a panel does not change the other pools.
build_pool <- function(family, k) {
  pn   <- PANELS[[k]]
  cols <- FAMILY_COLS[[family]]
  key  <- paste(family, paste(cols, collapse = ","), paste(deparse(pn$spec), collapse = " "),
                CFG$n_pool, CFG$seed, k,
                if (family == "kudrin") {
                  paste(CFG$K, CFG$N, CFG$gamma_max, CFG$search_alpha, CFG$dgp)
                } else {
                  paste(CFG$nu, CFG$cv, CFG$max_draws, format(CFG$smooth_pub$sw2), format(CFG$smooth_pub$theta))
                })
  file <- file.path(OUT_DIR, "pools", paste0(safe_id(pool_id(family, pn$id)), ".rds"))
  if (file.exists(file)) {
    pool <- readRDS(file)
    if (identical(attr(pool, "key"), key)) return(pool)
  }
  set.seed(CFG$seed + 1000 * k + (if (family == "kudrin") 1 else 2))
  message(sprintf("Simulating the %s pool for panel %s (%s studies) ...", family, pn$id,
                  format(CFG$n_pool, big.mark = ",", scientific = FALSE)))
  secs <- system.time({
    h <- draw_pi0(pn$spec, CFG$n_pool)
    if (family == "kudrin") {
      d    <- kudrin_draw(h, K = CFG$K, N = CFG$N, gamma_max = CFG$gamma_max,
                          alpha = CFG$search_alpha, dgp = CFG$dgp)
      pool <- cbind(honest = d[, "honest"], kudrin_threshold = d[, "threshold"],
                    kudrin_minimum = d[, "minimum"])
    } else {
      pool <- iid_pool(h, pn$spec, cols)
    }
  })[["elapsed"]]
  pool <- abs(pool[, c("honest", cols), drop = FALSE])   # all tests use |t| only
  attr(pool, "key") <- key
  saveRDS(pool, file)
  message(sprintf("  done in %.0f s", secs))
  pool
}

# Shares of p = 2 Phi(-|t|) <= 0.05 in every pool column, as a sanity check
describe_pools <- function() {
  for (id in names(POOLS)) {
    sh <- colMeans(2 * pnorm(-POOLS[[id]]) <= 0.05)
    message(sprintf("Pool %-28s share of p <= .05: %s", id,
                    paste(sprintf("%s %.4f", names(sh), sh), collapse = ", ")))
  }
}


## ---- 5. Run ------------------------------------------------------------------------------

POOLS     <- list()
POOL_KEYS <- list()
for (k in seq_along(PANELS)) {
  for (fam in names(FAMILY_COLS)) {
    id              <- pool_id(fam, PANELS[[k]]$id)
    POOLS[[id]]     <- build_pool(fam, k)
    POOL_KEYS[[id]] <- attr(POOLS[[id]], "key")
  }
}
describe_pools()

PROJ <- if (CFG$run_projection) proj_constants(file.path(OUT_DIR, "proj_constants.rds")) else NULL
if (!is.null(PROJ)) {
  message(sprintf("Projection test: epsilon_U = %.5f per sqrt(n) (paper: 0.00063).", PROJ$epsilon_U_pre))
}

# Methods of one family draw from the same pool; every cell's seed is shared across methods
CELLS <- make_cells(ROWS, PANELS, pool_of = function(rw, pn) pool_id(rw$family, pn$id))
pilot_estimate(CELLS)

cl  <- start_cluster()
res <- tryCatch({
  acc <- NULL
  for (S in CFG$S_values) {
    acc <- rbind(acc, run_cells(Filter(function(cc) cc$S == S, CELLS), cl))
    saveRDS(acc, file.path(OUT_DIR, paste0("results_", CFG$name, ".rds")))
    plot_figure(acc, S, ROWS, PANELS, file.path(OUT_DIR, sprintf("%s_S%d", CFG$name, as.integer(S))))
    message(sprintf("Figure for S = %d written to %s", as.integer(S), OUT_DIR))
  }
  acc
}, finally = {
  if (!is.null(cl)) parallel::stopCluster(cl)
})

saveRDS(res, file.path(OUT_DIR, paste0("results_", CFG$name, ".rds")))
write.csv(res, file.path(OUT_DIR, paste0("results_", CFG$name, ".csv")), row.names = FALSE)
tau0_check(res)
if (length(CFG$eig_floor) > 1) {
  message("Hybrid sensitivity to eig_floor (rejection rates):")
  print(res[, c("S", "row", "panel", "tau", "T_H_hybrid", sapply(CFG$eig_floor[-1], hyb_name))], digits = 3)
}
message("Figures and tables written to ", OUT_DIR)


## ---- 6. LaTeX tables ----------------------------------------------------------------------
## Paste at the end of kudrin_style_figures.R. One table per S, laid out like the figure: a
## block of rows per p-hacking method (the figure's rows), tau down each block, and a column
## group per panel (the figure's panels) with one column per test. Written to
## OUT_DIR/NAME_S1000.tex, OUT_DIR/NAME_S5000.tex; \input{} them in the paper, which needs
## \usepackage{booktabs}. Rates are among the repetitions in which a test could be computed;
## cells where Kudrin's tests (dagger) or the projection test (double dagger) failed in some
## repetitions are marked.
## Panel headers come from the panel's pi0; to set one yourself, add latex = "..." to that
## panel in CFG$panels, e.g. list(id = "N05", pi0 = ..., latex = "$h \\sim N(2, 0.5^2)$").
## Method headers come from the row titles; to set them yourself, add to CFG, e.g.
## CFG$row_latex <- c(max = "Maximum of 8 specifications").
## To rebuild the tables from saved results without rerunning the simulations, run
## sections 0-3, then res <- readRDS(file.path(OUT_DIR, paste0("results_", CFG$name, ".rds")))
## and this section.

TABLE_OPTS <- list(
  taus      = NULL,     # NULL: every tau in the results; or a subset such as seq(0, 1, by = 0.2)
  digits    = 3,        # decimals of the rejection rates
  bold_best = FALSE,    # TRUE: bold the highest rejection rate in each panel (tau > 0 only)
  fragment  = FALSE     # TRUE: write only the tabular, without the table environment and notes
)

TEST_LATEX <- c(proj = "Proj.", T_inf = "$T_\\infty$", T_H = "$T_H$",
                T_H_hybrid = "$T_H^{\\mathrm{hyb}}$", EKW_max = "EKW")

# Plain-text titles: escape the characters LaTeX treats as special (math in $...$ is kept)
latex_text <- function(x) {
  x <- gsub("([&%#_])", "\\\\\\1", x)
  x <- gsub("<=", "$\\\\le$", x)
  gsub(">=", "$\\\\ge$", x)
}

# Panel header from a pi0 spec (as pi0_label(), in LaTeX)
pi0_latex <- function(spec) {
  g <- function(x) format(signif(x, 4))
  s <- spec$shape
  if (identical(spec, PI0_ALIASES$kudrin_pi_hat)) return("$h \\sim \\hat{\\Pi}$")
  if (s == "point")   return(sprintf("$h = %s$", g(spec$h_center)))
  if (s == "normal")  return(sprintf("$h \\sim N(%s, %s^2)$", g(spec$h_center), g(spec$sigma_h)))
  if (s == "uniform") return(sprintf("$h \\sim U[%s, %s]$", g(spec$h_center - spec$sigma_h / 2),
                                     g(spec$h_center + spec$sigma_h / 2)))
  if (s == "chi2")    return(sprintf("$h \\sim \\chi^2_{%s}$", g(spec$h_center)))
  if (s == "poisson") return(sprintf("$h \\sim \\mathrm{Poisson}(%s)$", g(spec$h_center)))
  if (s == "gamma") {
    if (is.null(spec$gamma_scale)) {
      return(sprintf("$h \\sim \\Gamma(%s, \\mathrm{rate}\\ %s)$", g(spec$gamma_shape), g(spec$gamma_rate)))
    }
    return(sprintf("$h \\sim \\Gamma(%s, \\mathrm{scale}\\ %s)$", g(spec$gamma_shape), g(spec$gamma_scale)))
  }
  if (s == "double_normal") {
    return(sprintf("$h \\sim \\frac{1}{2} N(%s, %s^2) + \\frac{1}{2} N(%s, %s^2)$", g(spec$h_center),
                   g(spec$sigma_h), g(spec$h_center - 2), g(spec$sigma_h)))
  }
  if (s == "fitted") return("$h \\sim$ z-curve fit")
  "$h \\sim$ custom"
}

write_latex_tables <- function(res, opts = TABLE_OPTS) {
  tests <- intersect(c("proj", "T_inf", "T_H", "T_H_hybrid", "EKW_max"), names(res))
  kt    <- intersect(c("T_inf", "T_H", "T_H_hybrid"), tests)
  nt    <- length(tests)
  np    <- length(PANELS)
  ncols <- 1 + nt * np
  plab  <- vapply(seq_along(PANELS), function(k) {
    l <- CFG$panels[[k]]$latex
    if (is.null(l)) pi0_latex(PANELS[[k]]$spec) else l
  }, "")
  rlab  <- vapply(ROWS, function(rw) {
    if (!is.null(CFG$row_latex) && rw$id %in% names(CFG$row_latex)) CFG$row_latex[[rw$id]] else latex_text(rw$title)
  }, "")
  head <- c(paste0(" & ", paste(sprintf("\\multicolumn{%d}{c}{%s}", nt, plab), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + (seq_len(np) - 1) * nt, 1 + seq_len(np) * nt),
                  collapse = " "),
            paste0("$\\tau$ & ", paste(rep(TEST_LATEX[tests], np), collapse = " & "), " \\\\"))
  files <- character(0)
  
  for (S in sort(unique(res$S))) {
    rs   <- res[res$S == S, , drop = FALSE]
    taus <- sort(unique(rs$tau))
    if (!is.null(opts$taus)) taus <- taus[round(taus, 8) %in% round(opts$taus, 8)]
    tdig <- if (all(abs(taus * 10 - round(taus * 10)) < 1e-8)) 1 else 2
    dag_k <- FALSE
    dag_p <- FALSE
    body  <- character(0)
    for (i in seq_along(ROWS)) {
      if (i > 1) body <- c(body, "\\midrule")
      body <- c(body, sprintf("\\multicolumn{%d}{l}{\\textit{%s}} \\\\", ncols, rlab[i]))
      for (tau in taus) {
        entries <- character(0)
        for (k in seq_along(PANELS)) {
          r <- rs[rs$row == ROWS[[i]]$id & rs$panel == PANELS[[k]]$id & abs(rs$tau - tau) < 1e-8, ,
                  drop = FALSE]
          if (nrow(r) == 0) {
            entries <- c(entries, rep("", nt))
            next
          }
          v  <- vapply(tests, function(tn) as.numeric(r[[tn]][1]), 0)
          tx <- ifelse(is.finite(v), formatC(v, format = "f", digits = opts$digits), "--")
          if (isTRUE(opts$bold_best) && tau > 0 && any(is.finite(v))) {
            b <- is.finite(v) & v >= max(v, na.rm = TRUE) - 1e-12
            tx[b] <- paste0("\\textbf{", tx[b], "}")
          }
          if (length(kt) > 0 && isTRUE(r$kudrin_failed[1] > 0)) {
            j <- match(kt, tests)
            tx[j] <- paste0(tx[j], "$^\\dagger$")
            dag_k <- TRUE
          }
          if ("proj" %in% tests && isTRUE(r$proj_failed[1] > 0)) {
            j <- match("proj", tests)
            tx[j] <- paste0(tx[j], "$^\\ddagger$")
            dag_p <- TRUE
          }
          entries <- c(entries, tx)
        }
        body <- c(body, paste0(formatC(tau, format = "f", digits = tdig), " & ",
                               paste(entries, collapse = " & "), " \\\\"))
      }
    }
    
    tab <- c(sprintf("\\begin{tabular}{l%s}", strrep(strrep("c", nt), np)),
             "\\toprule", head, "\\midrule", body, "\\bottomrule", "\\end{tabular}")
    S_tex <- gsub(",", "{,}", format(S, big.mark = ",", scientific = FALSE), fixed = TRUE)
    notes <- paste0(
      sprintf("Rejection rates of %s\\%%-level tests over %s Monte Carlo repetitions; each meta-sample has ",
              format(100 * CFG$alpha), paste(format(unique(rs$reps), big.mark = ","), collapse = "/")),
      sprintf("$S = %s$ $t$-statistics. ", S_tex),
      "Blocks of rows: p-hacking method. $\\tau$: share of results that are p-hacked, so the $\\tau = 0$ ",
      "rows are sizes. Column groups: distribution of true effects $h$. ",
      "Proj.: projection test; $T_\\infty$, $T_H$, $T_H^{\\mathrm{hyb}}$: Kudrin (2025)",
      if ("EKW_max" %in% tests) "; EKW: highest rejection rate among the EKW (2022) tests" else "", ".",
      if (dag_k) " $^\\dagger$Kudrin's tests could not be computed in some repetitions; the rate is among the others." else "",
      if (dag_p) " $^\\ddagger$The projection test could not be computed in some repetitions; the rate is among the others." else "")
    caption <- sprintf("Rejection rates, $S = %s$%s", S_tex,
                       if (S == 1000) " (as in Kudrin's Fig.~2)" else if (S == 5000) " (as in Kudrin's Fig.~3)" else "")
    out <- if (isTRUE(opts$fragment)) tab else c(
      "\\begin{table}[htbp]", "\\centering",
      sprintf("\\caption{%s}", caption),
      sprintf("\\label{tab:%s_S%d}", CFG$name, as.integer(S)),
      "\\footnotesize", "\\setlength{\\tabcolsep}{4pt}",
      tab,
      "\\par\\smallskip",
      sprintf("\\begin{minipage}{\\linewidth}\\footnotesize\\textit{Notes:} %s\\end{minipage}", notes),
      "\\end{table}")
    f <- file.path(OUT_DIR, sprintf("%s_S%d.tex", CFG$name, as.integer(S)))
    writeLines(out, f)
    files <- c(files, f)
  }
  invisible(files)
}

TABLE_FILES <- write_latex_tables(res)
message("LaTeX tables written: ", paste(basename(TABLE_FILES), collapse = ", "))