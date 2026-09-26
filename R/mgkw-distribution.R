# ============================================================================
# MG-Kumaraswamy (MGKw) distribution: core functions
#
# Obtained by applying the MG transmutation map (Kumar et al., 2017) to the
# Kumaraswamy baseline. Support: x in (0, 1). Shape parameters a, b > 0.
# No additional parameters are introduced by the transmutation.
#
#   F(x) = exp( -(1-x^a)^b / (1 - (1-x^a)^b) )                      eq. (5)
#   f(x) = a b x^(a-1) (1-x^a)^(b-1) / (1-(1-x^a)^b)^2 * exp(...)   eq. (6)
#   Q(u) = [ 1 - ( ln(u) / (ln(u)-1) )^(1/b) ]^(1/a)                eq. (10)
#
# Author: Daniels Ebuka Marvelous
# ============================================================================

.check_ab <- function(a, b) {
  if (any(a <= 0) || any(b <= 0)) stop("Parameters a and b must be positive.")
}

#' Density of the MGKw distribution
dmgkw <- function(x, a, b, log = FALSE) {
  .check_ab(a, b)
  out <- rep(-Inf, length(x))
  ok  <- x > 0 & x < 1
  if (any(ok)) {
    xo  <- x[ok]
    # log1p/expm1 keep 1 - t accurate when x^a is tiny (else 1 - t rounds
    # to 0 and the density evaluates to Inf - Inf = NaN).
    l1  <- log1p(-xo^a)           # log(1 - x^a)
    t   <- exp(b * l1)            # (1 - x^a)^b, in (0,1)
    omt <- -expm1(b * l1)         # 1 - t
    ld  <- log(a) + log(b) + (a - 1) * log(xo) + (b - 1) * l1 -
           2 * log(omt) - t / omt
    out[ok] <- ld
  }
  if (log) out else exp(out)
}

#' Distribution function of the MGKw distribution
pmgkw <- function(q, a, b, lower.tail = TRUE) {
  .check_ab(a, b)
  out <- ifelse(q <= 0, 0, ifelse(q >= 1, 1, NA_real_))
  ok  <- q > 0 & q < 1
  if (any(ok)) {
    l1 <- b * log1p(-q[ok]^a)
    out[ok] <- exp(-exp(l1) / -expm1(l1))
  }
  if (lower.tail) out else 1 - out
}

#' Quantile function of the MGKw distribution (manuscript eq. 10)
qmgkw <- function(p, a, b) {
  .check_ab(a, b)
  if (any(p <= 0 | p >= 1)) stop("p must lie strictly in (0, 1).")
  lp <- log(p)
  (1 - (lp / (lp - 1))^(1 / b))^(1 / a)
}

#' Random generation by inversion
rmgkw <- function(n, a, b) {
  qmgkw(runif(n), a, b)
}

#' Survival function                                              eq. (9)
smgkw <- function(x, a, b) 1 - pmgkw(x, a, b)

#' Hazard rate function
hmgkw <- function(x, a, b) dmgkw(x, a, b) / smgkw(x, a, b)

#' Cumulative hazard function
Hmgkw <- function(x, a, b) -log(smgkw(x, a, b))

#' Reversed hazard rate function
rhmgkw <- function(x, a, b) dmgkw(x, a, b) / pmgkw(x, a, b)

#' r-th crude moment, by numerical integration of x^r f(x)
#' Cross-check for the series representation in eq. (11).
moment_mgkw <- function(a, b, r) {
  integrand <- function(x) x^r * dmgkw(x, a, b)
  stats::integrate(integrand, lower = 0, upper = 1,
                   rel.tol = 1e-10, subdivisions = 500L)$value
}

#' Moment-based summary measures (SD, CV, skewness, kurtosis)
moment_measures <- function(a, b) {
  m <- vapply(1:4, function(r) moment_mgkw(a, b, r), numeric(1))
  sd_  <- sqrt(m[2] - m[1]^2)
  cv   <- sqrt(m[2] / m[1]^2 - 1)
  cs   <- (m[3] - 3 * m[1] * m[2] + 2 * m[1]^3) / (m[2] - m[1]^2)^(3 / 2)
  ck   <- (m[4] - 4 * m[1] * m[3] + 6 * m[2] * m[1]^2 - 3 * m[1]^4) /
          (m[2] - m[1]^2)^2
  data.frame(mu1 = m[1], mu2 = m[2], mu3 = m[3], mu4 = m[4],
             SD = sd_, CV = cv, CS = cs, CK = ck)
}

#' Negative log-likelihood for a sample                           eq. (38)
nll_mgkw <- function(par, x) {
  a <- par[1]; b <- par[2]
  if (a <= 0 || b <= 0) return(1e10)
  ll <- sum(dmgkw(x, a, b, log = TRUE))
  if (!is.finite(ll)) return(1e10)
  -ll
}

#' Maximum likelihood estimation, BFGS as used in the manuscript
fit_mgkw <- function(x, start = c(1, 1)) {
  stats::optim(par = start, fn = nll_mgkw, x = x, method = "BFGS",
               hessian = TRUE)
}
