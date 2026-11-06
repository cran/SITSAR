#' The G0I distribution for SAR intensities
#'
#' Density, distribution function, quantile function and random generation
#' for the \eqn{\mathcal{G}_I^0(\alpha, \gamma, L)} distribution, the
#' multiplicative model for multilook SAR intensities with reciprocal gamma
#' backscatter and gamma speckle (Frery et al., 1997).
#'
#' The density is
#' \deqn{f(z) = \frac{L^L \Gamma(L - \alpha)}{\gamma^{\alpha} \Gamma(-\alpha)
#'   \Gamma(L)} \frac{z^{L - 1}}{(\gamma + L z)^{L - \alpha}}, \quad z > 0,}
#' with roughness \eqn{\alpha < 0}, scale \eqn{\gamma > 0} and number of
#' looks \eqn{L > 0}. If \eqn{Z \sim \mathcal{G}_I^0(\alpha, \gamma, L)} then
#' \eqn{-\alpha Z / \gamma} follows a Snedecor F distribution with
#' \eqn{2L} and \eqn{-2\alpha} degrees of freedom, which is used to compute
#' the distribution and quantile functions exactly.
#'
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param alpha Roughness parameter (\eqn{\alpha < 0}).
#' @param gamma Scale parameter (\eqn{\gamma > 0}).
#' @param L Number of looks (\eqn{L > 0}).
#' @param log,log.p Logical; if `TRUE`, probabilities/densities are given as
#'   logarithms.
#' @param lower.tail Logical; if `TRUE` (default) probabilities are
#'   \eqn{P(X \le x)}, otherwise \eqn{P(X > x)}.
#'
#' @return `dg0i()` gives the density, `pg0i()` the distribution function,
#'   `qg0i()` the quantile function and `rg0i()` generates random deviates.
#'
#' @references
#' Frery, A. C., Muller, H.-J., Yanasse, C. C. F. and Sant'Anna, S. J. S.
#' (1997). A model for extremely heterogeneous clutter. *IEEE Transactions
#' on Geoscience and Remote Sensing*, 35(3), 648--659.
#'
#' Nascimento, A. D. C., Cintra, R. J. and Frery, A. C. (2010). Hypothesis
#' testing in speckled data with stochastic distances. *IEEE Transactions on
#' Geoscience and Remote Sensing*, 48(1), 373--385.
#'
#' @seealso [G0I] for the R6 interface (estimation and entropy).
#'
#' @examples
#' set.seed(1)
#' z <- rg0i(5, alpha = -1.5, gamma = 2, L = 8)
#' dg0i(z, alpha = -1.5, gamma = 2, L = 8)
#' pg0i(c(0.5, 1, 2), alpha = -1.5, gamma = 2, L = 8)
#' qg0i(c(0.1, 0.5, 0.9), alpha = -1.5, gamma = 2, L = 8)
#'
#' @name dg0i
#' @aliases dg0i pg0i qg0i rg0i
NULL

.check_g0i <- function(alpha, gamma, L) {
  if (!is.numeric(alpha) || any(alpha >= 0)) stop("`alpha` must be negative.", call. = FALSE)
  if (!is.numeric(gamma) || any(gamma <= 0)) stop("`gamma` must be positive.", call. = FALSE)
  if (!is.numeric(L) || any(L <= 0)) stop("`L` must be positive.", call. = FALSE)
}

#' @rdname dg0i
#' @export
dg0i <- function(x, alpha, gamma, L, log = FALSE) {
  .check_g0i(alpha, gamma, L)
  out <- rep(-Inf, length(x))
  pos <- !is.na(x) & x > 0
  xp <- x[pos]
  out[pos] <- L * log(L) + lgamma(L - alpha) - alpha * log(gamma) -
    lgamma(-alpha) - lgamma(L) + (L - 1) * log(xp) -
    (L - alpha) * log(gamma + L * xp)
  out[is.na(x)] <- NA_real_
  if (log) out else exp(out)
}

#' @rdname dg0i
#' @export
pg0i <- function(q, alpha, gamma, L, lower.tail = TRUE, log.p = FALSE) {
  .check_g0i(alpha, gamma, L)
  stats::pf(pmax(q, 0) * (-alpha) / gamma, df1 = 2 * L, df2 = -2 * alpha,
            lower.tail = lower.tail, log.p = log.p)
}

#' @rdname dg0i
#' @export
qg0i <- function(p, alpha, gamma, L, lower.tail = TRUE, log.p = FALSE) {
  .check_g0i(alpha, gamma, L)
  gamma / (-alpha) * stats::qf(p, df1 = 2 * L, df2 = -2 * alpha,
                               lower.tail = lower.tail, log.p = log.p)
}

#' @rdname dg0i
#' @export
rg0i <- function(n, alpha, gamma, L) {
  .check_g0i(alpha, gamma, L)
  gamma / (-alpha) * stats::rf(n, df1 = 2 * L, df2 = -2 * alpha)
}


#' The Gamma distribution for SAR intensities
#'
#' Density, distribution function, quantile function and random generation
#' for the multilook intensity Gamma law parametrised by its mean \eqn{\mu}
#' and number of looks \eqn{L}:
#' \deqn{f(z) = \frac{L^L}{\mu^L \Gamma(L)} z^{L-1} \exp(-L z / \mu),
#' \quad z > 0.}
#' This is the model for homogeneous areas (fully developed speckle).
#'
#' @inheritParams dg0i
#' @param mu Mean (\eqn{\mu > 0}).
#'
#' @return Density, distribution function, quantile function or random
#'   deviates.
#'
#' @seealso [GammaSAR] for the R6 interface.
#'
#' @examples
#' dgammasar(c(0.5, 1, 2), mu = 1, L = 4)
#' rgammasar(3, mu = 1, L = 4)
#'
#' @name dgammasar
#' @aliases dgammasar pgammasar qgammasar rgammasar
NULL

#' @rdname dgammasar
#' @export
dgammasar <- function(x, mu, L, log = FALSE) {
  .check_positive_vector(mu, "mu"); .check_positive_vector(L, "L")
  stats::dgamma(x, shape = L, rate = L / mu, log = log)
}

#' @rdname dgammasar
#' @export
pgammasar <- function(q, mu, L, lower.tail = TRUE, log.p = FALSE) {
  .check_positive_vector(mu, "mu"); .check_positive_vector(L, "L")
  stats::pgamma(q, shape = L, rate = L / mu, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname dgammasar
#' @export
qgammasar <- function(p, mu, L, lower.tail = TRUE, log.p = FALSE) {
  .check_positive_vector(mu, "mu"); .check_positive_vector(L, "L")
  stats::qgamma(p, shape = L, rate = L / mu, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname dgammasar
#' @export
rgammasar <- function(n, mu, L) {
  .check_positive_vector(mu, "mu"); .check_positive_vector(L, "L")
  stats::rgamma(n, shape = L, rate = L / mu)
}


#' The generalized Gamma distribution for SAR intensities
#'
#' Density, distribution function, quantile function and random generation
#' for the generalized Gamma (Stacy) law with power \eqn{\alpha > 0}, shape
#' \eqn{\phi > 0} and inverse-scale \eqn{\mu > 0}:
#' \deqn{f(z) = \frac{\alpha \mu^{\alpha\phi}}{\Gamma(\phi)}
#'   z^{\alpha\phi - 1} \exp\{-(\mu z)^\alpha\}, \quad z > 0.}
#' For \eqn{\alpha = 1} it reduces to the Gamma law with shape \eqn{\phi}
#' and rate \eqn{\mu}; in SAR applications \eqn{\phi} plays the role of the
#' number of looks.
#'
#' @inheritParams dg0i
#' @param alpha Power parameter (\eqn{\alpha > 0}).
#' @param phi Shape parameter (\eqn{\phi > 0}).
#' @param mu Inverse-scale parameter (\eqn{\mu > 0}).
#'
#' @return Density, distribution function, quantile function or random
#'   deviates.
#'
#' @references
#' Li, H.-C., Hong, W., Wu, Y.-R. and Fan, P.-Z. (2011). On the empirical-
#' statistical modeling of SAR images with generalized gamma distribution.
#' *IEEE Journal of Selected Topics in Signal Processing*, 5(3), 386--397.
#'
#' @seealso [GGammaSAR] for the R6 interface.
#'
#' @examples
#' dggamma(c(0.5, 1, 2), alpha = 1.5, phi = 2, mu = 1)
#' rggamma(3, alpha = 1.5, phi = 2, mu = 1)
#'
#' @name dggamma
#' @aliases dggamma pggamma qggamma rggamma
NULL

.check_ggamma <- function(alpha, phi, mu) {
  .check_positive_vector(alpha, "alpha")
  .check_positive_vector(phi, "phi")
  .check_positive_vector(mu, "mu")
}

#' @rdname dggamma
#' @export
dggamma <- function(x, alpha, phi, mu, log = FALSE) {
  .check_ggamma(alpha, phi, mu)
  out <- rep(-Inf, length(x))
  pos <- !is.na(x) & x > 0
  xp <- x[pos]
  out[pos] <- log(alpha) + alpha * phi * log(mu) - lgamma(phi) +
    (alpha * phi - 1) * log(xp) - (mu * xp)^alpha
  out[is.na(x)] <- NA_real_
  if (log) out else exp(out)
}

#' @rdname dggamma
#' @export
pggamma <- function(q, alpha, phi, mu, lower.tail = TRUE, log.p = FALSE) {
  .check_ggamma(alpha, phi, mu)
  stats::pgamma((mu * pmax(q, 0))^alpha, shape = phi,
                lower.tail = lower.tail, log.p = log.p)
}

#' @rdname dggamma
#' @export
qggamma <- function(p, alpha, phi, mu, lower.tail = TRUE, log.p = FALSE) {
  .check_ggamma(alpha, phi, mu)
  stats::qgamma(p, shape = phi, lower.tail = lower.tail, log.p = log.p)^(1 / alpha) / mu
}

#' @rdname dggamma
#' @export
rggamma <- function(n, alpha, phi, mu) {
  .check_ggamma(alpha, phi, mu)
  stats::rgamma(n, shape = phi)^(1 / alpha) / mu
}
