# ============================================================================
# Competing distributions used in the manuscript.
#
# NOTE ON SUPPORT:
#   InvKw, GIKw, and MOEIKw all have support x > 0, NOT (0, 1).
#   The datasets used in the paper happen to fall in (0, 1), but these
#   distributions are fit as positive-support distributions.  The original
#   (now removed) implementation used a unit-interval formula for InvKw:
#       (1 - (1-q)^a / q^a)^b
#   which turns negative for q < ~0.5 and produces NaN with non-integer b.
#   The correct formulas (from the cited references) are given below.
# ============================================================================

# ---- Kumaraswamy  (unit-interval baseline, 2 params) ------------------------

dkw <- function(x, a, b, log = FALSE) {
  if (a <= 0 || b <= 0) stop("a and b must be positive.")
  out <- rep(if (log) -Inf else 0, length(x))
  ok  <- is.finite(x) & x > 0 & x < 1
  if (any(ok)) {
    ld <- log(a) + log(b) + (a - 1)*log(x[ok]) + (b - 1)*log1p(-x[ok]^a)
    out[ok] <- if (log) ld else exp(ld)
  }
  out
}
pkw <- function(q, a, b) {
  if (a <= 0 || b <= 0) stop("a and b must be positive.")
  ifelse(q <= 0, 0, ifelse(q >= 1, 1, 1 - (1 - q^a)^b))
}
nll_kw <- function(par, x) {
  if (length(par) != 2L || any(!is.finite(par)) || any(par <= 0)) return(1e100)
  ll <- sum(dkw(x, par[1], par[2], log = TRUE))
  if (is.finite(ll)) -ll else 1e100
}

# ---- Inverted Kumaraswamy  (IKw, x > 0, 2 params) --------------------------
# Reference: G(x) = [1 - (1+x)^(-a)]^b,  x > 0.
# Using log1p for numerical stability: (1+x)^(-a) = exp(-a * log1p(x)).

dikw <- function(x, a, b, log = FALSE) {
  if (a <= 0 || b <= 0) stop("a and b must be positive.")
  out <- rep(if (log) -Inf else 0, length(x))
  ok  <- is.finite(x) & x >= 0
  if (any(ok)) {
    z  <- 1 - exp(-a * log1p(x[ok]))   # [1 - (1+x)^(-a)]
    ld <- log(a) + log(b) - (a + 1)*log1p(x[ok]) + (b - 1)*log(z)
    out[ok] <- if (log) ld else exp(ld)
  }
  out
}
pikw <- function(q, a, b) {
  if (a <= 0 || b <= 0) stop("a and b must be positive.")
  out <- rep(NA_real_, length(q))
  out[q <= 0] <- 0
  ok <- is.finite(q) & q > 0
  if (any(ok)) out[ok] <- (1 - exp(-a * log1p(q[ok])))^b
  out[is.infinite(q) & q > 0] <- 1
  out
}
nll_ikw <- function(par, x) {
  if (length(par) != 2L || any(!is.finite(par)) || any(par <= 0)) return(1e100)
  ll <- sum(dikw(x, par[1], par[2], log = TRUE))
  if (is.finite(ll)) -ll else 1e100
}

# ---- Generalized Inverted Kumaraswamy  (GIKw, x > 0, 3 params) -------------
# F(x) = [1 - (1 + x^gamma)^(-alpha)]^beta,  x > 0.

dgikw <- function(x, alpha, beta, gamma, log = FALSE) {
  if (any(c(alpha, beta, gamma) <= 0) || any(!is.finite(c(alpha, beta, gamma))))
    stop("Parameters must be positive and finite.")
  out <- rep(if (log) -Inf else 0, length(x))
  ok  <- is.finite(x) & x > 0
  if (any(ok)) {
    xx  <- x[ok]
    xg  <- xx^gamma
    z   <- 1 - (1 + xg)^(-alpha)         # [1 - (1+x^g)^(-a)]
    ld  <- log(alpha) + log(beta) + log(gamma) +
           (gamma - 1)*log(xx) -
           (alpha + 1)*log1p(xg) +
           (beta - 1)*log(z)
    out[ok] <- if (log) ld else exp(ld)
  }
  out
}
pgikw <- function(q, alpha, beta, gamma) {
  if (any(c(alpha, beta, gamma) <= 0) || any(!is.finite(c(alpha, beta, gamma))))
    stop("Parameters must be positive and finite.")
  out <- rep(NA_real_, length(q)); out[q <= 0] <- 0
  ok  <- is.finite(q) & q > 0
  if (any(ok)) out[ok] <- (1 - (1 + q[ok]^gamma)^(-alpha))^beta
  out[is.infinite(q) & q > 0] <- 1
  out
}
nll_gikw <- function(par, x) {
  if (length(par) != 3L || any(!is.finite(par)) || any(par <= 0)) return(1e100)
  ll <- sum(dgikw(x, par[1], par[2], par[3], log = TRUE))
  if (is.finite(ll)) -ll else 1e100
}

# ---- Marshall-Olkin Extended Inverted Kumaraswamy  (MOEIKw, x > 0, 3 params)
# F(x) = G(x) / [kappa + (1 - kappa) * G(x)],  where G is IKw.

dmoeikw <- function(x, alpha, beta, kappa, log = FALSE) {
  if (any(c(alpha, beta, kappa) <= 0) || any(!is.finite(c(alpha, beta, kappa))))
    stop("Parameters must be positive and finite.")
  out <- rep(if (log) -Inf else 0, length(x))
  ok  <- is.finite(x) & x > 0
  if (any(ok)) {
    G  <- pikw(x[ok], alpha, beta)
    gd <- dikw(x[ok], alpha, beta, log = TRUE)
    ld <- log(kappa) + gd - 2*log(kappa + (1 - kappa)*G)
    out[ok] <- if (log) ld else exp(ld)
  }
  out
}
pmoeikw <- function(q, alpha, beta, kappa) {
  if (any(c(alpha, beta, kappa) <= 0) || any(!is.finite(c(alpha, beta, kappa))))
    stop("Parameters must be positive and finite.")
  G <- pikw(q, alpha, beta)
  G / (kappa + (1 - kappa)*G)
}
nll_moeikw <- function(par, x) {
  if (length(par) != 3L || any(!is.finite(par)) || any(par <= 0)) return(1e100)
  ll <- sum(dmoeikw(x, par[1], par[2], par[3], log = TRUE))
  if (is.finite(ll)) -ll else 1e100
}
