#' Shannon entropies of SAR and PolSAR models
#'
#' Closed-form Shannon entropies \eqn{H = E\{-\log f(Z)\}} for the models
#' used in SAR (intensity) and PolSAR (matrix) image analysis. All functions
#' are vectorised in the scalar parameters.
#'
#' \describe{
#'   \item{Gamma}{\eqn{H = L - \log L + \log\mu + \log\Gamma(L) +
#'     (1 - L)\psi^{(0)}(L)}.}
#'   \item{G0I}{(Ferreira and Nascimento, 2020)
#'     \eqn{H = -\log\frac{L^L\Gamma(L-\alpha)}{\Gamma(-\alpha)\gamma^\alpha
#'     \Gamma(L)} + (1 - L)[\psi^{(0)}(L) - \log L + \log\gamma -
#'     \psi^{(0)}(-\alpha)] + (L - \alpha)[\psi^{(0)}(L - \alpha) +
#'     \log\gamma - \psi^{(0)}(-\alpha)]}.}
#'   \item{Generalized Gamma}{\eqn{H = -\log\alpha + \log\Gamma(\phi) + \phi
#'     - (\phi - 1/\alpha)\psi^{(0)}(\phi) - \log\mu}.}
#'   \item{Scaled complex Wishart}{(Nascimento, Ferreira and Frery, 2023)
#'     \eqn{H = m\log|\Sigma| - m^2\log L + \log\Gamma_m(L) + mL +
#'     (m - L)\psi^{(0)}_m(L)}, with
#'     \eqn{\psi^{(0)}_m(L) = \sum_{i=0}^{m-1}\psi^{(0)}(L - i)}.}
#'   \item{G0 polarimetric}{(Nascimento, Ferreira and Frery, 2023) obtained
#'     from \eqn{E\log|Z| = m[\log(-\alpha-1) - \psi^{(0)}(-\alpha)] +
#'     \log|\Sigma| - m\log L + \psi^{(0)}_m(L)} and
#'     \eqn{E\log[L\,\mathrm{tr}(\Sigma^{-1}Z) - \alpha - 1] =
#'     \log(-\alpha-1) + \psi^{(0)}(mL - \alpha) - \psi^{(0)}(-\alpha)}.}
#' }
#'
#' @param mu Mean of the Gamma law, or inverse scale of the generalized
#'   Gamma law.
#' @param L Number of looks.
#' @param alpha Roughness (G0I, G0 polarimetric) or power (generalized
#'   Gamma) parameter.
#' @param gamma Scale parameter of the G0I law.
#' @param phi Shape parameter of the generalized Gamma law.
#' @param Sigma Hermitian positive definite covariance matrix.
#'
#' @return Numeric vector of entropies.
#'
#' @references
#' Ferreira, J. A. and Nascimento, A. D. C. (2020). Shannon entropy for the
#' G0I model: a new segmentation approach. *IEEE Journal of Selected Topics
#' in Applied Earth Observations and Remote Sensing*, 13, 2547--2553.
#' \doi{10.1109/JSTARS.2020.2997666}
#'
#' Nascimento, A. D. C., Ferreira, J. A. and Frery, A. C. (2023).
#' Unsupervised segmentation of PolSAR data with complex Wishart and G0m
#' distributions and Shannon entropy. *Statistics and Computing*, 33, 129.
#' \doi{10.1007/s11222-023-10298-0}
#'
#' @examples
#' entropy_gamma(mu = 1, L = 4)
#' entropy_g0i(alpha = -3, gamma = 2, L = 4)
#' entropy_ggamma(alpha = 1, phi = 4, mu = 4) # equals the Gamma entropy
#' S <- polsar_example_sigma("forest")
#' entropy_wishart(S, L = 4)
#' entropy_g0p(alpha = -5, Sigma = S, L = 4)
#'
#' @name entropy
NULL

#' @rdname entropy
#' @export
entropy_gamma <- function(mu, L) {
  .check_positive_vector(mu, "mu"); .check_positive_vector(L, "L")
  L - log(L) + log(mu) + lgamma(L) + (1 - L) * digamma(L)
}

#' @rdname entropy
#' @export
entropy_g0i <- function(alpha, gamma, L) {
  .check_g0i(alpha, gamma, L)
  -(L * log(L) + lgamma(L - alpha) - lgamma(-alpha) - alpha * log(gamma) -
      lgamma(L)) +
    (1 - L) * (digamma(L) - log(L) + log(gamma) - digamma(-alpha)) +
    (L - alpha) * (digamma(L - alpha) + log(gamma) - digamma(-alpha))
}

#' @rdname entropy
#' @export
entropy_ggamma <- function(alpha, phi, mu) {
  .check_ggamma(alpha, phi, mu)
  -log(alpha) + lgamma(phi) + phi - (phi - 1 / alpha) * digamma(phi) - log(mu)
}

#' @rdname entropy
#' @export
entropy_wishart <- function(Sigma, L) {
  .check_hpd(Sigma)
  m <- nrow(Sigma)
  .entropy_wishart_ld(.hlogdet(Sigma), L, m)
}

# Wishart entropy from log|Sigma| (vectorised)
.entropy_wishart_ld <- function(logdet, L, m) {
  m * logdet - m^2 * log(L) + .lmvgamma(L, m) + m * L +
    (m - L) * .mvpsigamma(L, m)
}

#' @rdname entropy
#' @export
entropy_g0p <- function(alpha, Sigma, L) {
  .check_hpd(Sigma)
  m <- nrow(Sigma)
  if (any(alpha >= -1)) stop("`alpha` must be smaller than -1.", call. = FALSE)
  .entropy_g0p_ld(alpha, .hlogdet(Sigma), L, m)
}

.entropy_g0p_ld <- function(alpha, logdet, L, m) {
  logc <- m * L * log(L) + lgamma(m * L - alpha) - L * logdet -
    .lmvgamma(L, m) - lgamma(-alpha) - alpha * log(-alpha - 1)
  e_logdet <- m * (log(-alpha - 1) - digamma(-alpha)) + logdet -
    m * log(L) + .mvpsigamma(L, m)
  e_logtr <- log(-alpha - 1) + digamma(m * L - alpha) - digamma(-alpha)
  -logc - (L - m) * e_logdet - (alpha - m * L) * e_logtr
}

# Asymptotic variance of the Wishart (and Gamma, m = 1) entropy estimator
# based on N observations; when L is estimated, the extra term comes from the
# Fisher information of L (orthogonal to Sigma).
.entropy_var_wishart <- function(L, m, N, L_estimated) {
  v <- m^3 / L
  if (L_estimated) {
    psi1 <- .mvpsigamma(L, m, deriv = 1L)
    dHdL <- m - m^2 / L + (m - L) * psi1
    v <- v + dHdL^2 / (psi1 - m / L)
  }
  v / N
}
