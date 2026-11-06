#' The scaled complex Wishart distribution
#'
#' Density and random generation for the scaled complex Wishart law
#' \eqn{\mathcal{W}_m(\Sigma, L)}, the reference model for multilook
#' PolSAR covariance (coherency) matrices in homogeneous areas:
#' \deqn{f(Z) = \frac{L^{mL} |Z|^{L-m}}{|\Sigma|^L \Gamma_m(L)}
#'   \exp\{-L\, \mathrm{tr}(\Sigma^{-1} Z)\},}
#' where \eqn{\Gamma_m(L) = \pi^{m(m-1)/2}\prod_{i=0}^{m-1}\Gamma(L-i)} and
#' \eqn{E(Z) = \Sigma}.
#'
#' Random matrices are generated with the complex Bartlett decomposition,
#' so non-integer numbers of looks \eqn{L > m - 1} are allowed.
#'
#' @param x PolSAR observations: a complex array `m x m x n`, a single
#'   Hermitian matrix or an `n x m(m+1)/2` matrix (see [as_polsar_vec()]).
#' @param n Number of observations.
#' @param Sigma Hermitian positive definite `m x m` matrix (the mean).
#' @param L Number of looks (\eqn{L > m - 1}).
#' @param log Logical; if `TRUE` the log-density is returned.
#' @param output Format of the generated sample: `"array"` (`m x m x n`)
#'   or `"vec"` (`n x m(m+1)/2`).
#'
#' @return `dwishart_sar()` returns a vector of (log-)densities;
#'   `rwishart_sar()` returns the simulated sample.
#'
#' @references
#' Goodman, N. R. (1963). Statistical analysis based on a certain
#' multivariate complex Gaussian distribution. *The Annals of Mathematical
#' Statistics*, 34(1), 152--177.
#'
#' Lee, J.-S., Grunes, M. R. and Kwok, R. (1994). Classification of
#' multi-look polarimetric SAR imagery based on complex Wishart
#' distribution. *International Journal of Remote Sensing*, 15(11),
#' 2299--2311.
#'
#' @seealso [WishartSAR] for the R6 interface.
#'
#' @examples
#' Sigma <- polsar_example_sigma("forest")
#' set.seed(1)
#' Z <- rwishart_sar(3, Sigma, L = 4)
#' dim(Z)
#' dwishart_sar(Z, Sigma, L = 4, log = TRUE)
#'
#' @name dwishart_sar
#' @aliases dwishart_sar rwishart_sar
NULL

#' @rdname dwishart_sar
#' @export
dwishart_sar <- function(x, Sigma, L, log = FALSE) {
  .check_hpd(Sigma)
  m <- nrow(Sigma)
  .check_scalar(L, "L", lower = m - 1)
  V <- .as_vec(x)
  ldz <- log(.hdet_vec(V, m))
  tr <- .htrace_vec(.hinv(Sigma), V)
  out <- m * L * log(L) - L * .hlogdet(Sigma) - .lmvgamma(L, m) +
    (L - m) * ldz - L * tr
  if (log) out else exp(out)
}

# Unscaled complex Wishart deviates (E = df * Sigma) via Bartlett.
.rcwishart_unscaled <- function(n, Sigma, df) {
  m <- nrow(Sigma)
  R <- .hsqrt(Sigma)
  df <- rep_len(df, n)
  out <- array(0 + 0i, dim = c(m, m, n))
  for (s in seq_len(n)) {
    A <- matrix(0 + 0i, m, m)
    for (i in seq_len(m)) {
      A[i, i] <- sqrt(stats::rchisq(1, df = 2 * (df[s] - i + 1)) / 2)
      if (i > 1) {
        j <- seq_len(i - 1)
        A[i, j] <- complex(real = stats::rnorm(i - 1), imaginary = stats::rnorm(i - 1)) / sqrt(2)
      }
    }
    B <- R %*% A
    W <- B %*% Conj(t(B))
    out[, , s] <- (W + Conj(t(W))) / 2
  }
  out
}

#' @rdname dwishart_sar
#' @export
rwishart_sar <- function(n, Sigma, L, output = c("array", "vec")) {
  output <- match.arg(output)
  .check_hpd(Sigma)
  m <- nrow(Sigma)
  .check_scalar(L, "L", lower = m - 1)
  out <- .rcwishart_unscaled(n, Sigma + 0i, L) / L
  if (output == "vec") as_polsar_vec(out) else out
}


#' The G0 polarimetric distribution
#'
#' Density and random generation for the polarimetric
#' \eqn{\mathcal{G}^0_m(\alpha, \Sigma, L)} law (also known as
#' \eqn{\mathcal{G}^0_P}), obtained from the multiplicative model
#' \eqn{Z = X Y} with \eqn{Y \sim \mathcal{W}_m(\Sigma, L)} and reciprocal
#' gamma texture \eqn{X \sim \Gamma^{-1}(-\alpha, -\alpha - 1)} with unit
#' mean (Freitas et al., 2005):
#' \deqn{f(Z) = \frac{L^{mL}\Gamma(mL - \alpha)|Z|^{L-m}}
#'   {|\Sigma|^L \Gamma_m(L) \Gamma(-\alpha)(-\alpha-1)^{\alpha}}
#'   \left[L\,\mathrm{tr}(\Sigma^{-1}Z) - \alpha - 1\right]^{\alpha - mL},}
#' with \eqn{\alpha < -1}. When \eqn{\alpha \to -\infty} it converges to
#' the scaled complex Wishart law.
#'
#' @inheritParams dwishart_sar
#' @param alpha Roughness parameter (\eqn{\alpha < -1}).
#'
#' @return `dg0p()` returns (log-)densities and `rg0p()` the simulated
#'   sample.
#'
#' @references
#' Freitas, C. C., Frery, A. C. and Correia, A. H. (2005). The polarimetric
#' G distribution for SAR data analysis. *Environmetrics*, 16(1), 13--31.
#'
#' Nascimento, A. D. C., Ferreira, J. A. and Frery, A. C. (2023).
#' Unsupervised segmentation of PolSAR data with complex Wishart and
#' G0m distributions and Shannon entropy. *Statistics and Computing*, 33,
#' 129. \doi{10.1007/s11222-023-10298-0}
#'
#' @seealso [G0Polarimetric] for the R6 interface.
#'
#' @examples
#' Sigma <- polsar_example_sigma("urban")
#' set.seed(1)
#' Z <- rg0p(3, alpha = -3, Sigma = Sigma, L = 4)
#' dg0p(Z, alpha = -3, Sigma = Sigma, L = 4, log = TRUE)
#'
#' @name dg0p
#' @aliases dg0p rg0p
NULL

#' @rdname dg0p
#' @export
dg0p <- function(x, alpha, Sigma, L, log = FALSE) {
  .check_hpd(Sigma)
  m <- nrow(Sigma)
  .check_scalar(alpha, "alpha", upper = -1)
  .check_scalar(L, "L", lower = m - 1)
  V <- .as_vec(x)
  ldz <- log(.hdet_vec(V, m))
  tr <- .htrace_vec(.hinv(Sigma), V)
  out <- m * L * log(L) + lgamma(m * L - alpha) - L * .hlogdet(Sigma) -
    .lmvgamma(L, m) - lgamma(-alpha) - alpha * log(-alpha - 1) +
    (L - m) * ldz + (alpha - m * L) * log(L * tr - alpha - 1)
  if (log) out else exp(out)
}

#' @rdname dg0p
#' @export
rg0p <- function(n, alpha, Sigma, L, output = c("array", "vec")) {
  output <- match.arg(output)
  .check_scalar(alpha, "alpha", upper = -1)
  Y <- rwishart_sar(n, Sigma, L)
  X <- (-alpha - 1) / stats::rgamma(n, shape = -alpha)
  out <- sweep(Y, 3, X, `*`)
  if (output == "vec") as_polsar_vec(out) else out
}


#' Compound complex Wishart distributions
#'
#' Density and random generation for the compound truncated Poisson complex
#' Wishart (CTPCW) and compound geometric complex Wishart (CGCW) laws
#' (Ferreira, Nascimento and Frery, 2022). They describe the PolSAR return
#' as a random sum \eqn{S = \sum_{i=1}^N Z_i} of independent (unscaled)
#' complex Wishart matrices \eqn{Z_i \sim \mathcal{W}^C_m(\Sigma, L)} with
#' \eqn{E(Z_i) = L\Sigma}, where the number of scatterers \eqn{N} follows a
#' zero-truncated Poisson law with mean parameter \eqn{\lambda} (CTPCW) or a
#' geometric law on \eqn{\{1, 2, \ldots\}} with success probability \eqn{p}
#' (CGCW). Their densities are
#' \deqn{f(S) = \sum_{k=1}^{\infty} P(N = k)
#'   \frac{|S|^{kL - m}}{|\Sigma|^{kL}\Gamma_m(kL)}
#'   \exp\{-\mathrm{tr}(\Sigma^{-1}S)\},}
#' which are able to describe multimodal intensities. The series is
#' evaluated in log-scale and truncated when the terms become negligible.
#'
#' @inheritParams dwishart_sar
#' @param lambda Mean parameter of the zero-truncated Poisson law
#'   (\eqn{\lambda > 0}).
#' @param p Success probability of the geometric law (\eqn{0 < p < 1}).
#' @param L Number of looks of each Wishart component (\eqn{L > m - 1}).
#' @param kmax Maximum number of terms of the series.
#'
#' @return `dctpcw()`/`dcgcw()` return (log-)densities; `rctpcw()`/`rcgcw()`
#'   return the simulated sample, with the simulated numbers of components in
#'   the attribute `"N"`.
#'
#' @references
#' Ferreira, J. A., Nascimento, A. D. C. and Frery, A. C. (2022). PolSAR
#' models with multimodal intensities. *Remote Sensing*, 14(20), 5083.
#' \doi{10.3390/rs14205083}
#'
#' @seealso [CompoundWishart] for the R6 interface and EM estimation.
#'
#' @examples
#' Sigma <- polsar_example_sigma("forest")
#' set.seed(1)
#' S <- rctpcw(3, lambda = 2, Sigma = Sigma, L = 4)
#' attr(S, "N")
#' dctpcw(S, lambda = 2, Sigma = Sigma, L = 4, log = TRUE)
#' S2 <- rcgcw(3, p = 0.6, Sigma = Sigma, L = 4)
#' dcgcw(S2, p = 0.6, Sigma = Sigma, L = 4, log = TRUE)
#'
#' @name compound_wishart
#' @aliases dctpcw rctpcw dcgcw rcgcw
NULL

# log P(N = k) for the counting laws
.log_pmf_count <- function(k, counting, par) {
  if (counting == "poisson") {
    k * log(par) - lgamma(k + 1) - log(expm1(par))
  } else {
    log(par) + (k - 1) * log1p(-par)
  }
}

# log of the terms  log P(N=k) + kL (log|S| - log|Sigma|) - lmvgamma(kL, m)
# for all observations; returns a list with the matrix of terms (n x K)
# and the common part of the log density.
.compound_terms <- function(V, counting, par, Sigma, L, kmax = 500L) {
  m <- nrow(Sigma)
  ldz <- log(.hdet_vec(V, m))
  lds <- .hlogdet(Sigma)
  tr <- .htrace_vec(.hinv(Sigma), V)
  common <- -m * ldz - tr
  ratio <- ldz - lds
  terms <- NULL
  k <- 0L
  repeat {
    k <- k + 1L
    tk <- .log_pmf_count(k, counting, par) + k * L * ratio - .lmvgamma(k * L, m)
    terms <- cbind(terms, tk)
    # the terms are unimodal in k: stop once the current term is negligible
    # with respect to the largest one for every observation
    if (k >= 2L) {
      mx <- apply(terms, 1, max)
      if (all(tk < mx - 40, na.rm = TRUE)) break
    }
    if (k >= kmax) break
  }
  list(terms = terms, common = common)
}

.dcompound <- function(x, counting, par, Sigma, L, log, kmax) {
  .check_hpd(Sigma)
  m <- nrow(Sigma)
  .check_scalar(L, "L", lower = m - 1)
  V <- .as_vec(x)
  tt <- .compound_terms(V, counting, par, Sigma + 0i, L, kmax)
  out <- tt$common + apply(tt$terms, 1, .logsumexp)
  if (log) out else exp(out)
}

.rcompound <- function(n, counting, par, Sigma, L, output) {
  .check_hpd(Sigma)
  m <- nrow(Sigma)
  .check_scalar(L, "L", lower = m - 1)
  N <- if (counting == "poisson") {
    stats::qpois(stats::runif(n, stats::dpois(0, par), 1), par)
  } else {
    stats::rgeom(n, par) + 1L
  }
  N <- pmax(N, 1)
  out <- .rcwishart_unscaled(n, Sigma + 0i, N * L)
  if (output == "vec") out <- as_polsar_vec(out)
  attr(out, "N") <- N
  out
}

#' @rdname compound_wishart
#' @export
dctpcw <- function(x, lambda, Sigma, L, log = FALSE, kmax = 500L) {
  .check_scalar(lambda, "lambda", lower = 0)
  .dcompound(x, "poisson", lambda, Sigma, L, log, kmax)
}

#' @rdname compound_wishart
#' @export
rctpcw <- function(n, lambda, Sigma, L, output = c("array", "vec")) {
  .check_scalar(lambda, "lambda", lower = 0)
  .rcompound(n, "poisson", lambda, Sigma, L, match.arg(output))
}

#' @rdname compound_wishart
#' @export
dcgcw <- function(x, p, Sigma, L, log = FALSE, kmax = 500L) {
  .check_scalar(p, "p", lower = 0, upper = 1)
  .dcompound(x, "geometric", p, Sigma, L, log, kmax)
}

#' @rdname compound_wishart
#' @export
rcgcw <- function(n, p, Sigma, L, output = c("array", "vec")) {
  .check_scalar(p, "p", lower = 0, upper = 1)
  .rcompound(n, "geometric", p, Sigma, L, match.arg(output))
}


#' Example PolSAR covariance matrices
#'
#' Typical \eqn{3 \times 3} covariance matrices (scaled to unit magnitude)
#' for urban, forest and pasture areas, taken from the simulation studies of
#' the authors' dissertation and thesis. They are handy to simulate PolSAR
#' data.
#'
#' @param type One of `"urban"`, `"forest"` or `"pasture"`.
#'
#' @return A complex Hermitian positive definite \eqn{3 \times 3} matrix.
#'
#' @examples
#' polsar_example_sigma("forest")
#'
#' @export
polsar_example_sigma <- function(type = c("urban", "forest", "pasture")) {
  type <- match.arg(type)
  S <- switch(
    type,
    urban = matrix(c(962892, 19171 - 3579i, -154638 + 191388i,
                     19171 + 3579i, 56707, -5798 + 16812i,
                     -154638 - 191388i, -5798 - 16812i, 472251),
                   ncol = 3, byrow = TRUE),
    forest = matrix(c(360932, 11050 + 3759i, 63896 + 1581i,
                      11050 - 3759i, 98960, 6593 + 6868i,
                      63896 - 1581i, 6593 - 6868i, 208843),
                    ncol = 3, byrow = TRUE),
    pasture = matrix(c(32556, 556 + 787i, 24046 - 27287i,
                       556 - 787i, 1647, -146 - 482i,
                       24046 + 27287i, -146 + 482i, 61028),
                     ncol = 3, byrow = TRUE)
  )
  S / 1e6
}
