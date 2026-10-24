################################################################################
#  kudrin_figures_engine.R
#
#  Shared machinery for
#    reproduce_kudrin_fig2_fig3.R   Kudrin (2025) Figs. 2 and 3, both rows, plus the projection test
#    kudrin_style_figures.R         new versions of those figures for any Pi0 and p-hacking method
#
#  The method is that of your two replication scripts (and of EKW and Kudrin). For each panel,
#  a pool of studies is simulated once; each study has an honest t-statistic and one p-hacked
#  t-statistic per p-hacking method. A meta-sample draws S studies from the pool with replacement
#  and reports, for each, the p-hacked t-statistic with probability tau and the honest one
#  otherwise. Every meta-sample goes to your kudrin_tests() and run_test_weight_fast(), called
#  with exactly the arguments of the replication scripts. Nothing here changes how a test works:
#  draw_meta(), ekw_rejections(), one_rep(), summarise_cell() and proj_constants() are the
#  replication scripts' functions, except that a meta-sample now names the pool column to report.
#
#  What is new (orchestration and graphics only):
#    - A cell is one (S, row, panel, tau); a row is a p-hacking method, that is, a pool column.
#    - Every cell has its own RNG seed, which depends on S, the panel and tau but not on the row.
#      Rows that share a pool therefore test the same meta-samples at every tau (common random
#      numbers); at tau = 0 they agree exactly. A resumed run gives the same numbers as an
#      uninterrupted one (on the same number of cores). The draws differ from those of your
#      earlier, separate runs.
#    - New colours, markers, line types and axes, and the legend in a strip below the panels.
#
#  Source after your code files. Uses these objects from the global environment:
#    CFG, OUT_DIR, POOLS (named list of pool matrices), POOL_KEYS (same names), PROJ.
################################################################################


## ---- Projection-test constants (your create_basis_stable / compute_epsilons) ---------------

proj_constants <- function(file) {
  key <- paste(CFG$num_coeffs, CFG$L, CFG$numgrid,
               if (is.null(CFG$epsilon_U_pre)) "computed" else CFG$epsilon_U_pre)
  if (file.exists(file)) {
    pc <- readRDS(file)
    if (identical(pc$key, key)) return(pc)
  }
  eps <- CFG$epsilon_U_pre
  if (is.null(eps)) {
    message("Computing epsilon_U(L = ", CFG$L, ", numgrid = ", CFG$numgrid, ") once; this is slow ...")
    eps <- compute_epsilons(CFG$L, nx = CFG$numgrid)$epsilon_U
  }
  pc <- list(key = key, epsilon_U_pre = eps,
             U = create_basis_stable(CFG$num_coeffs, L = CFG$L, numgrid = CFG$numgrid, sigma_Y = 1))
  saveRDS(pc, file)
  pc
}


## ---- One meta-sample ----------------------------------------------------------------------

# S studies drawn with replacement; each reports column `column` with probability tau and the
# honest t-statistic otherwise. Same random draws, in the same order, as the replication scripts.
draw_meta <- function(pool, column, tau, S) {
  idx  <- sample.int(nrow(pool), S, replace = TRUE)
  hack <- runif(S) < tau
  ifelse(hack, pool[idx, column], pool[idx, "honest"])
}

hyb_name <- function(f) paste0("T_H_hybrid_floor", formatC(f, format = "e", digits = 0))

# EKW (2022) tests on (0, 0.15], called as in your run_sims(); solver failures count as
# non-rejections (EKW 2025, footnote 19). The robust caliper test is not in Tests.R.
ekw_rejections <- function(tstats) {
  ps   <- 2 * pnorm(-abs(tstats))
  pmn  <- min(ps[ps > 0])
  pmx  <- 0.15
  J    <- CFG$ekw_J
  safe <- function(expr) {
    v <- tryCatch(suppressWarnings(expr), error = function(e) NA_real_)
    if (length(v) != 1 || !is.numeric(v)) NA_real_ else as.numeric(v)
  }
  pv <- c(EKW_Binomial      = safe(Binomial(ps, 0.04, 0.05, "c")),
          EKW_Discontinuity = safe(Discontinuity_test(ps, 0.05)),
          EKW_Fisher        = safe(Fisher(ps, pmn, pmx)),
          EKW_LCM           = safe(LCM(ps, pmn, pmx)),
          EKW_CS1           = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 1, 0)),
          EKW_CSUB          = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 0, 1)),
          EKW_CS2B          = safe(CoxShi(ps, seq_along(ps), pmn, pmx, J, 2, 1)))
  setNames(as.numeric(!is.na(pv) & pv <= CFG$alpha), names(pv))
}

# All tests on one meta-sample of cell `cell`. Returns a named numeric vector (1 = reject).
one_rep <- function(r, cell) {
  tstats <- draw_meta(POOLS[[cell$pool]], cell$column, cell$tau, cell$S)
  out    <- c(T_inf = NA_real_, T_H = NA_real_, T_H_hybrid = NA_real_, hyb_restricted = NA_real_)

  kt <- tryCatch(kudrin_tests(tstats, boots = CFG$kudrin_boots, alpha = CFG$alpha,
                              eig_floor = CFG$eig_floor[1]),
                 error = function(e) NULL)
  if (!is.null(kt)) {
    out[c("T_inf", "T_H", "T_H_hybrid")] <- as.numeric(kt$pval[c("T_inf", "T_H", "T_H_hybrid")] <= CFG$alpha)
    out[["hyb_restricted"]] <- as.numeric(kt$restricted)
  }
  for (f in CFG$eig_floor[-1]) {             # hybrid only; T_H and the hybrid use no bootstrap
    k2 <- tryCatch(kudrin_tests(tstats, boots = 1, alpha = CFG$alpha, eig_floor = f),
                   error = function(e) NULL)
    out[hyb_name(f)] <- if (is.null(k2)) NA_real_ else as.numeric(k2$pval[["T_H_hybrid"]] <= CFG$alpha)
  }

  if (CFG$run_projection) {
    pr <- tryCatch(run_test_weight_fast(data.frame(t = tstats, title = as.character(seq_along(tstats))),
                                        numcoeffs = CFG$num_coeffs, sigma_Y = 1, shift = CFG$shift,
                                        L = CFG$L, numgrid = CFG$numgrid,
                                        boots = if (CFG$proj_mode == "warp") 1 else CFG$proj_boots,
                                        U = PROJ$U, epsilon_U_pre = PROJ$epsilon_U_pre,
                                        p = CFG$p_weight, J0 = CFG$J0),
                   error = function(e) NULL)
    proj <- c(proj = NA_real_, proj_resid = NA_real_, proj_boot1 = NA_real_, proj_epsU = NA_real_)
    if (!is.null(pr)) {
      proj <- c(proj       = as.numeric(pr$pval <= CFG$alpha),
                proj_resid = as.numeric(pr$resid),
                proj_boot1 = as.numeric(pr$resids_boot[1]),
                proj_epsU  = as.numeric(pr$epsilon_U))
    }
    out <- c(out, proj)
  }

  if (CFG$run_ekw) out <- c(out, ekw_rejections(tstats))
  out
}


## ---- Cells --------------------------------------------------------------------------------

safe_id <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)

# rows:   list of list(id, title, column); the column of the pool that the row reports
# panels: list of list(id, title)
# pool_of(row, panel): name of the pool in POOLS that the (row, panel) cell draws from
# The seed depends on S, the panel and tau, not on the row (common random numbers across rows).
make_cells <- function(rows, panels, pool_of) {
  cells <- list()
  for (fi in seq_along(CFG$S_values)) {
    for (rw in rows) {
      for (k in seq_along(panels)) {
        for (tau in CFG$taus) {
          S  <- CFG$S_values[fi]
          pn <- panels[[k]]
          cells[[length(cells) + 1]] <- list(
            S = S, row = rw$id, panel = pn$id, tau = tau,
            pool = pool_of(rw, pn), column = rw$column,
            seed = CFG$seed + 1000000 * fi + 10000 * k + round(1000 * tau),
            file = file.path(OUT_DIR, "cells",
                             sprintf("cell_S%d_%s_%s_tau%04d.rds", as.integer(S), safe_id(rw$id),
                                     safe_id(pn$id), as.integer(round(1000 * tau)))))
        }
      }
    }
  }
  cells
}

cell_key <- function(cell) {
  paste(POOL_KEYS[[cell$pool]], cell$column, cell$S, cell$tau, cell$seed, CFG$reps,
        CFG$kudrin_boots, paste(CFG$eig_floor, collapse = ","), CFG$run_projection, CFG$proj_mode,
        CFG$proj_boots, CFG$num_coeffs, CFG$J0, CFG$p_weight, CFG$shift, CFG$L, CFG$numgrid,
        if (is.null(PROJ)) "no_proj" else PROJ$epsilon_U_pre, CFG$run_ekw, CFG$ekw_J)
}

# PSOCK cluster, so it also runs on Windows; each worker sources your files.
start_cluster <- function() {
  if (CFG$n_cores <= 1) return(NULL)
  cl <- parallel::makePSOCKcluster(CFG$n_cores)
  parallel::clusterExport(cl, "CODE_FILES", envir = globalenv())
  invisible(parallel::clusterEvalQ(cl, { for (f in CODE_FILES) source(f); NULL }))
  parallel::clusterExport(cl, c("CFG", "POOLS", "PROJ", "draw_meta", "hyb_name", "ekw_rejections",
                                "one_rep"), envir = globalenv())
  cl
}

run_cell <- function(cell, cl) {
  if (CFG$resume && file.exists(cell$file)) {
    done <- readRDS(cell$file)
    if (identical(done$key, cell_key(cell))) return(done)
  }
  t0 <- Sys.time()
  if (is.null(cl)) {
    set.seed(cell$seed)
    L <- lapply(seq_len(CFG$reps), one_rep, cell = cell)
  } else {
    parallel::clusterSetRNGStream(cl, cell$seed)
    L <- parallel::parLapply(cl, seq_len(CFG$reps), one_rep, cell = cell)
  }
  nm   <- names(L[[1]])
  M    <- do.call(rbind, lapply(L, function(v) v[nm]))
  done <- list(key = cell_key(cell), S = cell$S, row = cell$row, panel = cell$panel, tau = cell$tau,
               M = M, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(done, cell$file)
  done
}


## ---- Rejection rates ----------------------------------------------------------------------

summarise_cell <- function(done) {
  M    <- done$M
  rate <- function(col) mean(M[, col], na.rm = TRUE)
  s <- data.frame(S = done$S, row = done$row, panel = done$panel, tau = done$tau, reps = nrow(M),
                  T_inf = rate("T_inf"), T_H = rate("T_H"), T_H_hybrid = rate("T_H_hybrid"),
                  hybrid_restricted = rate("hyb_restricted"),
                  kudrin_failed = mean(is.na(M[, "T_inf"])), secs = done$secs,
                  stringsAsFactors = FALSE)
  for (f in CFG$eig_floor[-1]) s[[hyb_name(f)]] <- rate(hyb_name(f))
  if (CFG$run_projection) {
    if (CFG$proj_mode == "warp") {
      ok     <- is.finite(M[, "proj_resid"]) & is.finite(M[, "proj_boot1"])
      crit   <- quantile(M[ok, "proj_boot1"], 1 - CFG$alpha, names = FALSE)
      s$proj <- mean(M[ok, "proj_resid"] > crit + M[ok, "proj_epsU"])
    } else {
      s$proj <- rate("proj")
    }
    s$proj_failed <- mean(is.na(M[, "proj_resid"]))
  }
  if (CFG$run_ekw) {
    ek <- grep("^EKW_", colnames(M), value = TRUE)
    for (e in ek) s[[e]] <- rate(e)
    s$EKW_max <- max(unlist(s[ek]))
  }
  s
}

run_cells <- function(cells, cl) {
  acc <- NULL
  for (cell in cells) {
    s   <- summarise_cell(run_cell(cell, cl))
    acc <- rbind(acc, s)
    message(sprintf("S = %d  %-16s %-10s tau = %.2f:  T_inf %.3f   T_H %.3f   hybrid %.3f%s   (%.0f s)",
                    as.integer(cell$S), cell$row, cell$panel, cell$tau, s$T_inf, s$T_H, s$T_H_hybrid,
                    if (CFG$run_projection) sprintf("   projection %.3f", s$proj) else "", s$secs))
  }
  acc
}

# Run time of one meta-sample at each S, and the implied total. Also prints the sup-norm bias
# bound of kudrin_tests(), as the replication scripts did.
pilot_estimate <- function(cells) {
  Ss    <- sort(unique(vapply(cells, function(cc) cc$S, numeric(1))))
  first <- cells[[1]]
  set.seed(CFG$seed + 1)
  invisible(one_rep(1, first))               # warm-up: builds kudrin_tests' cached bounds
  per <- vapply(Ss, function(S) {
    cc <- first
    cc$S <- S
    cc$tau <- 1
    system.time(one_rep(1, cc))[["elapsed"]]
  }, numeric(1))
  n_by_S <- vapply(Ss, function(S) sum(vapply(cells, function(cc) cc$S == S, logical(1))), numeric(1))
  bbar   <- kudrin_tests(draw_meta(POOLS[[first$pool]], "honest", 0, min(Ss)), boots = 1)$bias_bound_inf
  hours  <- sum(per * n_by_S) * CFG$reps / 3600 / max(1, CFG$n_cores)
  message(sprintf("Sup-norm bias bound %.3e (expected about 6.55e-04).", bbar))
  message("One meta-sample takes ", paste(sprintf("%.1f s at S = %d", per, as.integer(Ss)), collapse = ", "),
          sprintf(", so %d cells x %d repetitions should take about %.1f hours on %d core(s).",
                  length(cells), CFG$reps, hours, CFG$n_cores))
  invisible(per)
}

# tau = 0 tests honest meta-samples only. Rows that share a pool draw the same meta-samples, so
# their entries agree exactly (when both cells come from runs on the same number of cores).
tau0_check <- function(res) {
  cols <- intersect(c("S", "panel", "row", "T_inf", "T_H", "T_H_hybrid", "proj"), names(res))
  z    <- res[abs(res$tau) < 1e-8, cols, drop = FALSE]
  z    <- z[order(z$S, z$panel, z$row), , drop = FALSE]
  message("tau = 0 (honest meta-samples only). Rows that share a pool should agree exactly:")
  print(z, digits = 3, row.names = FALSE)
  invisible(z)
}


## ---- Comparison with published curves ---------------------------------------------------
## pub: list(tau, tau_ekw, reps, rows = list(<row id> = list(<panel id> = list(T_inf, T_H,
## T_H_hybrid, EKW_max)))). Your rates minus his on his tau grid; the flag counts taus where the
## gap exceeds 3 standard errors (both simulations plus 0.002 of slack).

compare_with_published <- function(res, S, pub) {
  out <- list()
  for (rw in intersect(unique(res$row), names(pub$rows))) {
    for (pn in intersect(unique(res$panel), names(pub$rows[[rw]]))) {
      r <- res[res$S == S & res$row == rw & res$panel == pn, , drop = FALSE]
      i <- match(round(r$tau, 8), round(pub$tau, 8))
      r <- r[!is.na(i), , drop = FALSE]
      i <- i[!is.na(i)]
      if (nrow(r) == 0) next
      for (tn in c("T_inf", "T_H", "T_H_hybrid")) {
        ours <- r[[tn]]
        his  <- pub$rows[[rw]][[pn]][[tn]][i]
        se   <- sqrt(ours * (1 - ours) / r$reps + his * (1 - his) / pub$reps + 0.002^2)
        out[[length(out) + 1]] <- data.frame(
          S = S, row = rw, panel = pn, test = tn, n_tau = nrow(r),
          mean_diff = mean(ours - his), max_abs_diff = max(abs(ours - his)),
          tau_of_max = r$tau[which.max(abs(ours - his))],
          n_beyond_3se = sum(abs(ours - his) > 3 * se), stringsAsFactors = FALSE)
      }
    }
  }
  do.call(rbind, out)
}

# His published curves as a results table (to redraw his panels with this code)
published_as_res <- function(pub, S) {
  out <- list()
  for (rw in names(pub$rows)) {
    for (pn in names(pub$rows[[rw]])) {
      p <- pub$rows[[rw]][[pn]]
      out[[length(out) + 1]] <- data.frame(S = S, row = rw, panel = pn, tau = pub$tau, reps = pub$reps,
                                           T_inf = p$T_inf, T_H = p$T_H, T_H_hybrid = p$T_H_hybrid,
                                           stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

pub_panel <- function(pub, row, panel) {
  if (is.null(pub) || is.null(pub$rows[[row]]) || is.null(pub$rows[[row]][[panel]])) return(NULL)
  c(list(tau = pub$tau, tau_ekw = pub$tau_ekw), pub$rows[[row]][[panel]])
}


## ---- Figures ------------------------------------------------------------------------------
## Okabe-Ito colours, with markers, line types and axes that differ from Kudrin's figures:
## T_inf blue, solid, open triangles; T_H orange, dashed, plus signs; T_H^hybrid reddish
## purple, dot-dash, open diamonds; projection test black, solid, filled circles; EKW maximum
## grey, long dashes; 5% level black, dotted. Light grid, ticks outside on the left and bottom
## axes, and one legend in a strip below the panels.

STYLE <- list(
  sig        = list(col = "#000000", lty = "dotted",   lwd = 0.8),
  EKW_max    = list(col = "#8C8C8C", lty = "longdash", lwd = 1.1),
  T_inf      = list(col = "#0072B2", lty = "solid",    lwd = 1.1, pch = 2,  cex = 0.62, plwd = 0.9),
  T_H        = list(col = "#E69F00", lty = "dashed",   lwd = 1.1, pch = 3,  cex = 0.62, plwd = 0.9),
  T_H_hybrid = list(col = "#CC79A7", lty = "dotdash",  lwd = 1.1, pch = 5,  cex = 0.62, plwd = 0.9),
  proj       = list(col = "#000000", lty = "solid",    lwd = 1.6, pch = 16, cex = 0.55, plwd = 1.0)
)
AXIS_COL  <- "grey25"
GRID_COL  <- "grey90"
LEG_LABEL <- list(sig        = "5% level",
                  EKW_max    = "Max power of EKW tests",
                  T_inf      = expression(italic(T)[infinity]),
                  T_H        = expression(italic(T)[italic(H)]),
                  T_H_hybrid = expression(italic(T)[italic(H)]^hybrid),
                  proj       = "Projection test")
LEG_ORDER <- c("sig", "EKW_max", "T_inf", "T_H", "T_H_hybrid", "proj")

lighten <- function(col, f = 0.65) {     # opaque tint (works on every device)
  m <- col2rgb(col) / 255
  rgb(t(m + (1 - m) * f))
}

draw_series <- function(x, y, st, faint = FALSE) {
  ok   <- is.finite(y)
  col  <- if (faint) lighten(st$col) else st$col
  lwd  <- if (faint) 0.8 * st$lwd else st$lwd
  lines(x[ok], y[ok], col = col, lty = st$lty, lwd = lwd)
  if (!is.null(st$pch)) {
    plwd <- if (faint) 0.8 * st$plwd else st$plwd
    points(x[ok], y[ok], col = col, pch = st$pch, cex = st$cex, lwd = plwd, xpd = TRUE)
  }
  invisible(NULL)
}

# One panel. r: results for one (S, row, panel); pub: that panel's published curves or NULL.
# The EKW line comes from r when the EKW tests were run, otherwise from pub when available.
draw_panel <- function(r, title, pub = NULL, overlay_published = FALSE, add_projection = TRUE) {
  r <- r[order(r$tau), , drop = FALSE]
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  abline(h = c(0.25, 0.5, 0.75), v = c(0.25, 0.5, 0.75), col = GRID_COL, lwd = 0.5)
  at   <- c(0, 0.25, 0.5, 0.75, 1)
  labs <- c("0", ".25", ".5", ".75", "1")
  for (side in 1:2) {
    axis(side, at = at, labels = labs, lwd = 0, lwd.ticks = 0.5, col.ticks = AXIS_COL,
         col.axis = AXIS_COL, cex.axis = 0.8)
  }
  box(lwd = 0.6, col = AXIS_COL)
  mtext(expression("Extent of selective reporting ("*tau*")"), side = 1, line = 1.35, cex = 0.78)
  mtext("Rejection rate", side = 2, line = 1.9, cex = 0.78, las = 0)
  mtext(title, side = 3, line = 0.25, cex = 0.8)

  abline(h = CFG$alpha, col = STYLE$sig$col, lty = STYLE$sig$lty, lwd = STYLE$sig$lwd)
  if ("EKW_max" %in% names(r)) {
    draw_series(r$tau, r$EKW_max, STYLE$EKW_max)
  } else if (!is.null(pub) && !is.null(pub$EKW_max)) {
    draw_series(pub$tau_ekw, pub$EKW_max, STYLE$EKW_max)
  }
  tests <- c("T_inf", "T_H", "T_H_hybrid")
  if (overlay_published && !is.null(pub)) {
    for (tn in tests) draw_series(pub$tau, pub[[tn]], STYLE[[tn]], faint = TRUE)
  }
  for (tn in tests) {
    if (tn %in% names(r)) draw_series(r$tau, r[[tn]], STYLE[[tn]])
  }
  if (add_projection && ("proj" %in% names(r))) draw_series(r$tau, r$proj, STYLE$proj)
  invisible(NULL)
}

open_device <- function(file, ext, height, width = 6.6) {
  if (ext == "pdf") {
    pdf(file, width = width, height = height, pointsize = 8)
  } else if (isTRUE(capabilities("cairo"))) {
    png(file, width = width, height = height, units = "in", res = 300, pointsize = 8, type = "cairo")
  } else {
    png(file, width = width, height = height, units = "in", res = 300, pointsize = 8)
  }
}

# Legend in a strip across the bottom of the device (drawn last, over the whole page)
draw_legend_strip <- function(keys, ncol) {
  sty  <- STYLE[keys]
  labs <- do.call(c, lapply(keys, function(k) as.expression(LEG_LABEL[[k]])))
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend("bottom", inset = c(0, 0.004), legend = labs, ncol = ncol, bty = "n", xpd = NA,
         col    = vapply(sty, function(s) s$col, ""),
         lty    = vapply(sty, function(s) s$lty, ""),
         lwd    = vapply(sty, function(s) s$lwd, 0),
         pch    = vapply(sty, function(s) if (is.null(s$pch)) NA_real_ else s$pch, 0),
         pt.cex = vapply(sty, function(s) if (is.null(s$cex)) 1 else s$cex, 0),
         pt.lwd = vapply(sty, function(s) if (is.null(s$plwd)) 1 else s$plwd, 0),
         cex = 0.8, seg.len = 2.8, x.intersp = 0.6, y.intersp = 1.05)
}

# One figure: rows x panels at sample size S. rows, panels: lists of list(id, title).
# pub: published curves for this S (or NULL). note: a line of text above the panels.
plot_figure <- function(res, S, rows, panels, stem, pub = NULL, overlay_published = FALSE,
                        add_projection = TRUE, note = NULL) {
  nr   <- length(rows)
  np   <- length(panels)
  mid  <- ceiling(np / 2)
  keys <- c("sig",
            if ("EKW_max" %in% names(res) || !is.null(pub)) "EKW_max",
            "T_inf", "T_H", "T_H_hybrid",
            if (add_projection && "proj" %in% names(res)) "proj")
  keys <- intersect(LEG_ORDER, keys)
  ncol <- if (np >= 2) 3 else 2
  leg_lines <- 0.9 + 1.1 * ceiling(length(keys) / ncol)
  top_lines <- if (is.null(note)) 0 else 1.3
  line_in   <- 8 * 1.2 / 72                         # height of one margin line at pointsize 8
  height    <- 2.05 * nr + (leg_lines + top_lines) * line_in
  width     <- max(3.3, 2.2 * np)
  for (ext in c("pdf", "png")) {
    open_device(paste0(stem, ".", ext), ext, height = height, width = width)
    par(mfrow = c(nr, np))                          # separate call: mfrow resets cex
    par(cex = 1, mar = c(2.6, 3.0, 2.5, 0.8), oma = c(leg_lines, 0, top_lines, 0),
        mgp = c(1.5, 0.35, 0), tcl = -0.22, las = 1, family = "sans")
    for (i in seq_len(nr)) {
      for (k in seq_len(np)) {
        r <- res[res$S == S & res$row == rows[[i]]$id & res$panel == panels[[k]]$id, , drop = FALSE]
        draw_panel(r, panels[[k]]$title, pub_panel(pub, rows[[i]]$id, panels[[k]]$id),
                   overlay_published = overlay_published, add_projection = add_projection)
        if (k == mid) mtext(rows[[i]]$title, side = 3, line = 1.2, font = 2, cex = 1.05)
      }
    }
    if (!is.null(note)) mtext(note, side = 3, outer = TRUE, line = 0.3, cex = 0.78)
    draw_legend_strip(keys, ncol)
    dev.off()
  }
  invisible(paste0(stem, c(".pdf", ".png")))
}
