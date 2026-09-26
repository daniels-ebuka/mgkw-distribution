# ============================================================================
# Application: fitting MGKw to two real datasets and comparing with
# competing distributions.
#
# Datasets
# --------
# 1. COVID-19 mortality rate, France (1 Jan – 20 Feb 2021), n = 51
#    Source: Table 3 of the manuscript.
# 2. Unit capacity factors estimated by the P3 algorithm, n = 22
#    Source: Table 7 of the manuscript; Caramanis, M., Stremel, J., Fleck, W.
#    & Daniel, S. (1983). Probabilistic production costing: an investigation
#    of alternative algorithms. International Journal of Electrical Power &
#    Energy Systems, 5(2), 75-86.
#
# Competing models
# ----------------
#   Kw      — Kumaraswamy (unit interval, 2 params)
#   IKw     — Inverted Kumaraswamy   (x > 0, 2 params)
#   GIKw    — Generalized Inverted Kumaraswamy (x > 0, 3 params)
#   MOEIKw  — Marshall–Olkin Extended IKw      (x > 0, 3 params)
#
# Note on information criteria:
#   AIC/BIC/HQIC use the standard formulas (-2*logLik + penalty). Densities on
#   (0,1) can exceed 1, so the maximised log-likelihood is positive here and
#   the criteria are negative, as in Tables 5 and 9 of the manuscript.
#
# Note on goodness-of-fit:
#   CVM and AD are the Chen & Balakrishnan (1995) modified statistics W* and
#   A*, as computed by AdequacyModel::goodness.fit. These reproduce the
#   MGKw values in Table 6 exactly.
#
# Note on competitor support:
#   IKw, GIKw and MOEIKw are defined on x > 0, not (0,1).  Although the
#   datasets happen to lie in (0,1), these models are fit using the
#   correct positive-support formulations (see competing-distributions.R).
# ============================================================================

source(file.path("R", "mgkw-distribution.R"))
source(file.path("R", "competing-distributions.R"))

# ---- Datasets ---------------------------------------------------------------

covid_france <- c(
  0.0995, 0.0525, 0.0615, 0.0455, 0.1474, 0.3373, 0.1087, 0.1055, 0.2235,
  0.0633, 0.0565, 0.2577, 0.1345, 0.0843, 0.1023, 0.2296, 0.0691, 0.0505,
  0.1434, 0.2326, 0.1089, 0.1206, 0.2242, 0.0786, 0.0587, 0.1516, 0.2070,
  0.1170, 0.1141, 0.2705, 0.0793, 0.0635, 0.1474, 0.2345, 0.1131, 0.1129,
  0.2054, 0.0600, 0.0534, 0.1422, 0.2235, 0.0908, 0.1092, 0.1958, 0.0580,
  0.0502, 0.1229, 0.1738, 0.0917, 0.0787, 0.1654
)  # Table 3: COVID-19 mortality rate, France, 1 Jan – 20 Feb 2021, n = 51

p3_algorithm <- c(
  0.010, 0.014, 0.019, 0.026, 0.036, 0.044, 0.056, 0.067, 0.078, 0.097,
  0.118, 0.118, 0.207, 0.334, 0.399, 0.503, 0.557, 0.716, 0.759, 0.800,
  0.853, 0.874
)  # Table 7: unit capacity factors, P3 algorithm, n = 22

# ---- Information criteria (standard definitions) ----------------------------

ic_table <- function(loglik, k, n) {
  aic  <- -2 * loglik + 2 * k
  caic <- aic + (2 * k * (k + 1)) / (n - k - 1)
  bic  <- -2 * loglik + k * log(n)
  hqic <- -2 * loglik + 2 * k * log(log(n))
  c(AIC = aic, CAIC = caic, BIC = bic, HQIC = hqic)
}

# ---- Goodness-of-fit statistics ---------------------------------------------

# Chen & Balakrishnan (1995) W* (Cramér–von Mises) and A* (Anderson–Darling),
# plus the Kolmogorov–Smirnov statistic.
gof_stats <- function(x, pfun, ...) {
  n  <- length(x)
  i  <- seq_len(n)
  v  <- pfun(sort(x), ...)
  v  <- pmax(pmin(v, 1 - 1e-12), 1e-12)
  y  <- stats::qnorm(v)
  u  <- stats::pnorm((y - mean(y)) / stats::sd(y))

  W   <- sum((u - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)
  A   <- -n - mean((2 * i - 1) * log(u) + (2 * n + 1 - 2 * i) * log(1 - u))
  cvm <- W * (1 + 0.5 / n)
  ad  <- A * (1 + 0.75 / n + 2.25 / n^2)

  # Both datasets contain tied values, which triggers a harmless ks.test
  # warning; the p-value is approximate in any case because parameters
  # were estimated from the same data.
  ks <- suppressWarnings(stats::ks.test(x, pfun, ...))

  c(CVM = unname(cvm), AD = unname(ad),
    KS  = unname(ks$statistic), KS_p = ks$p.value)
}

# ---- Standard errors from Hessian -------------------------------------------

se_from_hessian <- function(hessian) {
  tryCatch(sqrt(diag(solve(hessian))), error = function(e) rep(NA_real_, nrow(hessian)))
}

# ---- Multi-start optimiser for the competing models -------------------------
# The three-parameter likelihoods have local optima: from a single start of
# (1, 1, 1) the GIKw fit to the P3 data stops at logLik = 4.43, whereas the
# global maximum found here is about 6.54. Every competitor is therefore
# fitted from a grid of starting values and the best solution is kept.

fit_positive <- function(x, fn, p) {
  grid <- if (p == 2) {
    expand.grid(c(0.1, 0.5, 1, 5, 20), c(0.1, 0.5, 1, 5, 20))
  } else {
    expand.grid(c(0.05, 0.1, 0.5, 1, 5, 20), c(0.3, 1, 3, 10),
                c(0.3, 0.5, 1, 2, 5))
  }
  best <- NULL
  for (i in seq_len(nrow(grid))) {
    f <- tryCatch(
      stats::optim(par = unlist(grid[i, ]), fn = fn, x = x,
                   method = "L-BFGS-B", lower = rep(1e-8, p),
                   control = list(maxit = 10000)),
      error = function(e) NULL)
    if (!is.null(f) && (is.null(best) || f$value < best$value)) best <- f
  }
  best$par     <- unname(best$par)
  best$hessian <- stats::optimHess(best$par, fn, x = x)
  best
}

# ---- Fit all models to one dataset ------------------------------------------

fit_all <- function(x, label) {
  n <- length(x)
  cat(sprintf("\n=== %s  (n = %d) ===\n", label, n))

  # --- MGKw ---
  fit_mg  <- fit_mgkw(x, start = c(0.5, 0.5))
  ll_mg   <- -fit_mg$value
  par_mg  <- fit_mg$par
  se_mg   <- se_from_hessian(fit_mg$hessian)
  cat(sprintf("MGKw MLE:  a = %.4f (SE %.4f),  b = %.4f (SE %.4f),  logLik = %.4f\n",
              par_mg[1], se_mg[1], par_mg[2], se_mg[2], ll_mg))

  # --- Kumaraswamy (unit-interval) ---
  fit_kw  <- fit_positive(x, nll_kw, 2)
  ll_kw   <- -fit_kw$value

  # --- Inverted Kumaraswamy (x > 0) ---
  fit_ikw <- fit_positive(x, nll_ikw, 2)
  ll_ikw  <- -fit_ikw$value

  # --- Generalized Inverted Kumaraswamy (x > 0) ---
  fit_gik <- fit_positive(x, nll_gikw, 3)
  ll_gik  <- -fit_gik$value

  # --- Marshall-Olkin Extended IKw (x > 0) ---
  fit_moe <- fit_positive(x, nll_moeikw, 3)
  ll_moe  <- -fit_moe$value

  # --- Assemble results table ---
  models <- list(
    MGKw   = list(ll = ll_mg,  k = 2, par = par_mg,       pfun = pmgkw),
    Kw     = list(ll = ll_kw,  k = 2, par = fit_kw$par,   pfun = pkw),
    IKw    = list(ll = ll_ikw, k = 2, par = fit_ikw$par,  pfun = pikw),
    GIKw   = list(ll = ll_gik, k = 3, par = fit_gik$par,  pfun = pgikw),
    MOEIKw = list(ll = ll_moe, k = 3, par = fit_moe$par,  pfun = pmoeikw)
  )

  rows <- lapply(names(models), function(nm) {
    m   <- models[[nm]]
    ic  <- ic_table(m$ll, m$k, n)
    gof <- tryCatch(
      do.call(gof_stats, c(list(x = x, pfun = m$pfun), as.list(m$par))),
      error = function(e) c(CVM = NA, AD = NA, KS = NA, KS_p = NA)
    )
    data.frame(
      Model  = nm,
      LogLik = round(m$ll,        4),
      AIC    = round(ic["AIC"],   4),
      CAIC   = round(ic["CAIC"],  4),
      BIC    = round(ic["BIC"],   4),
      HQIC   = round(ic["HQIC"],  4),
      CVM    = round(gof["CVM"],  4),
      AD     = round(gof["AD"],   4),
      KS     = round(gof["KS"],   4),
      KS_p   = round(gof["KS_p"], 4),
      stringsAsFactors = FALSE, row.names = NULL
    )
  })

  tab <- do.call(rbind, rows)
  print(tab, row.names = FALSE)

  dir.create("output", showWarnings = FALSE)
  out_file <- file.path(
    "output",
    paste0("application_", gsub("[^A-Za-z0-9]+", "_", label), ".csv")
  )
  write.csv(tab, out_file, row.names = FALSE)
  cat(sprintf("Results saved to %s\n", out_file))

  invisible(tab)
}

# ---- Run --------------------------------------------------------------------

fit_all(covid_france, "COVID-19 mortality rate, France")
fit_all(p3_algorithm, "Unit capacity factors, P3 algorithm")
