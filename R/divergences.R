#' Divergences between Gaussian and complex Wishart models
#'
#' Closed-form divergences between two multivariate normal laws
#' \eqn{N_p(\mu_1, \Sigma_1)} and \eqn{N_p(\mu_2, \Sigma_2)}, and between two
#' scaled complex Wishart laws \eqn{\mathcal{W}_m(\Sigma_1, L)} and
#' \eqn{\mathcal{W}_m(\Sigma_2, L)} with a common number of looks. They are
#' the building blocks of the divergence-based classifiers of Ferreira,
#' Coelho and Nascimento (2021).
#'
#' Available measures (`type`):
#' \describe{
#'   \item{`"kl"`}{Kullback-Leibler divergence \eqn{D_{KL}(\theta_1 \| \theta_2)}.}
#'   \item{`"jeffreys"`}{Jeffreys divergence
#'     \eqn{D_{KL}(\theta_1 \| \theta_2) + D_{KL}(\theta_2 \| \theta_1)}.}
#'   \item{`"bhattacharyya"`}{Bhattacharyya distance
#'     \eqn{-\log \int \sqrt{f_1 f_2}}.}
#'   \item{`"renyi"`}{Rényi divergence of order \eqn{\alpha} (`order`),
#'     in one of four forms (`kind`): 1, \eqn{D^\alpha_R(\theta_1\|\theta_2)};
#'     2, \eqn{D^\alpha_R(\theta_2\|\theta_1)}; 3, the symmetrised mean
#'     \eqn{[D^\alpha_R(\theta_1\|\theta_2) + D^\alpha_R(\theta_2\|\theta_1)]/2};
#'     and 4, \eqn{c \log\{[K(\theta_1\|\theta_2) + K(\theta_2\|\theta_1)]/2\}},
#'     where \eqn{K(\theta_1\|\theta_2) = \int f_1^\alpha f_2^{1-\alpha}}.}
#' }
#' For Gaussian laws the Rényi divergence uses the normalisation of Pardo
#' (2006), \eqn{c = 1/[\alpha(\alpha - 1)]}, which gives
#' \deqn{D^\alpha_R = \frac{1}{2}\Delta^\top M^{-1}\Delta -
#'   \frac{1}{2\alpha(\alpha-1)}\log\frac{|M|}{|\Sigma_1|^{1-\alpha}
#'   |\Sigma_2|^{\alpha}},}
#' with \eqn{\Delta = \mu_1 - \mu_2} and \eqn{M = \alpha\Sigma_2 +
#' (1-\alpha)\Sigma_1}. For Wishart laws the normalisation of Frery,
#' Nascimento and Cintra (2014), \eqn{c = 1/(\alpha - 1)}, is used:
#' \deqn{D^\alpha_R = \frac{L}{\alpha - 1}\log\frac{|(\alpha\Sigma_1^{-1} +
#'   (1-\alpha)\Sigma_2^{-1})^{-1}|}{|\Sigma_1|^\alpha|\Sigma_2|^{1-\alpha}}.}
#' Both converge to the Kullback-Leibler divergence when
#' \eqn{\alpha \to 1}.
#'
#' @param mu1,mu2 Mean vectors.
#' @param Sigma1,Sigma2 Covariance matrices (real symmetric for Gaussian
#'   laws, Hermitian for Wishart laws).
#' @param L Number of looks (common to both Wishart laws).
#' @param type Divergence: `"kl"`, `"jeffreys"`, `"bhattacharyya"` or
#'   `"renyi"`.
#' @param order Order \eqn{\alpha \in (0, 1)} of the Rényi divergence.
#' @param kind Form of the Rényi divergence (1 to 4).
#'
#' @return A single non-negative number.
#'
#' @references
#' Ferreira, J. A., Coelho, H. and Nascimento, A. D. C. (2021). A family of
#' divergence-based classifiers for polarimetric synthetic aperture radar
#' (PolSAR) imagery vector and matrix features. *International Journal of
#' Remote Sensing*, 42(4), 1201--1229. \doi{10.1080/01431161.2020.1826060}
#'
#' Frery, A. C., Nascimento, A. D. C. and Cintra, R. J. (2014). Analytic
#' expressions for stochastic distances between relaxed complex Wishart
#' distributions. *IEEE Transactions on Geoscience and Remote Sensing*,
#' 52(2), 1213--1226.
#'
#' Pardo, L. (2006). *Statistical Inference Based on Divergence Measures*.
#' Chapman & Hall/CRC.
#'
#' @examples
#' S1 <- diag(3); S2 <- diag(c(2, 1, 0.5))
#' divergence_gaussian(c(0, 0, 0), S1, c(1, 0, 0), S2, type = "kl")
#' divergence_gaussian(c(0, 0, 0), S1, c(1, 0, 0), S2, type = "renyi",
#'                     order = 0.5, kind = 3)
#'
#' W1 <- polsar_example_sigma("forest"); W2 <- polsar_example_sigma("urban")
#' divergence_wishart(W1, W2, L = 4, type = "jeffreys")
#' divergence_wishart(W1, W2, L = 4, type = "renyi", order = 0.9)
#'
#' @name divergences
NULL

.check_div_args <- function(type, order, kind) {
  type <- match.arg(type, c("kl", "jeffreys", "bhattacharyya", "renyi"))
  if (type == "renyi") {
    .check_scalar(order, "order", lower = 0, upper = 1)
    if (!kind %in% 1:4) stop("`kind` must be 1, 2, 3 or 4.", call. = FALSE)
  }
  type
}

#' @rdname divergences
#' @export
divergence_gaussian <- function(mu1, Sigma1, mu2, Sigma2,
                                type = c("kl", "jeffreys", "bhattacharyya", "renyi"),
                                order = 0.5, kind = 1) {
  type <- .check_div_args(type[1], order, kind)
  .div_gauss(as.numeric(mu1), as.matrix(Sigma1), as.numeric(mu2),
             as.matrix(Sigma2), type, order, kind)
}

.div_gauss <- function(mu1, S1, mu2, S2, type, order, kind) {
  p <- length(mu1)
  d <- mu1 - mu2
  ld1 <- .hlogdet(S1); ld2 <- .hlogdet(S2)
  kl <- function(m1, A, ldA, m2, B, ldB) {
    Bi <- solve(B)
    dd <- m1 - m2
    0.5 * (sum(Bi * t(A)) + sum(dd * (Bi %*% dd)) - p + ldB - ldA)
  }
  # log K(theta_a || theta_b; alpha) = log int f_a^alpha f_b^(1-alpha)
  logK <- function(dd, A, ldA, B, ldB, a) {
    M <- a * B + (1 - a) * A
    a * (a - 1) / 2 * sum(dd * solve(M, dd)) -
      0.5 * (.hlogdet(M) - (1 - a) * ldA - a * ldB)
  }
  switch(type,
    kl = kl(mu1, S1, ld1, mu2, S2, ld2),
    jeffreys = kl(mu1, S1, ld1, mu2, S2, ld2) + kl(mu2, S2, ld2, mu1, S1, ld1),
    bhattacharyya = {
      M <- (S1 + S2) / 2
      sum(d * solve(M, d)) / 8 + 0.5 * (.hlogdet(M) - 0.5 * (ld1 + ld2))
    },
    renyi = {
      a <- order
      c0 <- 1 / (a * (a - 1))
      k12 <- logK(d, S1, ld1, S2, ld2, a)
      k21 <- logK(-d, S2, ld2, S1, ld1, a)
      switch(kind,
             c0 * k12,
             c0 * k21,
             c0 * (k12 + k21) / 2,
             c0 * (.logsumexp(c(k12, k21)) - log(2)))
    })
}

#' @rdname divergences
#' @export
divergence_wishart <- function(Sigma1, Sigma2, L,
                               type = c("kl", "jeffreys", "bhattacharyya", "renyi"),
                               order = 0.5, kind = 1) {
  type <- .check_div_args(type[1], order, kind)
  .check_hpd(Sigma1, "Sigma1"); .check_hpd(Sigma2, "Sigma2")
  .check_scalar(L, "L", lower = 0)
  .div_wish(Sigma1 + 0i, Sigma2 + 0i, L, type, order, kind)
}

.div_wish <- function(S1, S2, L, type, order, kind,
                      ld1 = .hlogdet(S1), ld2 = .hlogdet(S2),
                      S1i = .hinv(S1), S2i = .hinv(S2)) {
  m <- nrow(S1)
  tr <- function(A, B) Re(sum(A * t(B)))
  # log int f_1^a f_2^(1-a) for scaled complex Wishart laws
  logK <- function(Ai, ldA, Bi, ldB, a) {
    -L * (.hlogdet(a * Ai + (1 - a) * Bi) + a * ldA + (1 - a) * ldB)
  }
  switch(type,
    kl = L * (tr(S2i, S1) - m + ld2 - ld1),
    jeffreys = L * (tr(S2i, S1) + tr(S1i, S2) - 2 * m),
    bhattacharyya = -logK(S1i, ld1, S2i, ld2, 0.5),
    renyi = {
      a <- order
      c0 <- 1 / (a - 1)
      k12 <- logK(S1i, ld1, S2i, ld2, a)
      k21 <- logK(S2i, ld2, S1i, ld1, a)
      switch(kind,
             c0 * k12,
             c0 * k21,
             c0 * (k12 + k21) / 2,
             c0 * (.logsumexp(c(k12, k21)) - log(2)))
    })
}


#' Stochastic distances between G0I distributions
#'
#' Computes stochastic distances between two \eqn{\mathcal{G}_I^0} laws by
#' numerical integration (Nascimento, Cintra and Frery, 2010).
#'
#' @param par_1,par_2 Parameter vectors `c(alpha, gamma, L)`.
#' @param p Order of the Rényi distance, in (0, 1).
#' @param type One of `"KL"`, `"Renyi"`, `"Bhattacharyya"`, `"Hellinger"`,
#'   `"Arithmetic-Geometric"`, `"Triangular"` or `"Harmonic-Mean"`.
#'
#' @return An object of class `distancesg0i` with the parameters, the type
#'   and the `distance_value`.
#'
#' @references
#' Nascimento, A. D. C., Cintra, R. J. and Frery, A. C. (2010). Hypothesis
#' testing in speckled data with stochastic distances. *IEEE Transactions on
#' Geoscience and Remote Sensing*, 48(1), 373--385.
#'
#' @examples
#' distancesg0i(c(-1.5, 2, 4), c(-1.5, 2, 4))
#' res <- distancesg0i(c(-1.5, 2, 8), c(-3, 2, 4), type = "Hellinger")
#' res
#'
#' @export
distancesg0i <- function(par_1, par_2, p = 1 / 2,
                         type = c("KL", "Renyi", "Bhattacharyya", "Hellinger",
                                  "Arithmetic-Geometric", "Triangular", "Harmonic-Mean")) {
  type <- match.arg(type)
  if (length(par_1) != 3 || length(par_2) != 3 || par_1[1] >= 0 || par_2[1] >= 0 ||
      par_1[2] <= 0 || par_2[2] <= 0 || par_1[3] <= 0 || par_2[3] <= 0) {
    stop("Parameters must be c(alpha < 0, gamma > 0, L > 0).", call. = FALSE)
  }
  .check_scalar(p, "p", lower = 0, upper = 1)
  lf1 <- function(x) dg0i(x, par_1[1], par_1[2], par_1[3], log = TRUE)
  lf2 <- function(x) dg0i(x, par_2[1], par_2[2], par_2[3], log = TRUE)
  int <- function(f) stats::integrate(function(x) {
    v <- f(x)
    v[!is.finite(v)] <- 0
    v
  }, 0, Inf, subdivisions = 1000L, rel.tol = 1e-8, stop.on.error = FALSE)$value
  bc <- function() int(function(x) exp((lf1(x) + lf2(x)) / 2))
  tri <- function() int(function(x) {
    a <- exp(lf1(x)); b <- exp(lf2(x))
    ifelse(a + b > 0, (a - b)^2 / (a + b), 0)
  })
  distance <- switch(
    type,
    KL = int(function(x) (exp(lf1(x)) - exp(lf2(x))) * (lf1(x) - lf2(x))) / 2,
    Renyi = (log(int(function(x) exp(p * lf1(x) + (1 - p) * lf2(x)))) +
               log(int(function(x) exp(p * lf2(x) + (1 - p) * lf1(x))))) / (2 * (p - 1)),
    Bhattacharyya = -log(bc()),
    Hellinger = 1 - bc(),
    `Arithmetic-Geometric` = int(function(x) {
      a <- exp(lf1(x)); b <- exp(lf2(x))
      ((a + b) / 2) * (log((a + b) / 2) - (lf1(x) + lf2(x)) / 2)
    }),
    Triangular = tri(),
    `Harmonic-Mean` = -log(1 - tri() / 2)
  )
  structure(list(par_1 = par_1, par_2 = par_2, type = type,
                 distance_value = max(distance, 0)),
            class = "distancesg0i")
}

#' @export
print.distancesg0i <- function(x, ...) {
  cat("\nG0I stochastic distance\n\n")
  cat("Type:     ", x$type, "\n")
  cat("Distance: ", format(x$distance_value, digits = 6), "\n")
  invisible(x)
}

#' Table of G0I stochastic distances between two image regions
#'
#' Fits the \eqn{\mathcal{G}_I^0} law to each intensity channel of two
#' rectangular regions of an image and computes all the stochastic
#' distances of [distancesg0i()] between them.
#'
#' @param image A [SARImage] or an array (`nrow x ncol x layers`).
#' @param marks1,marks2 Region limits as `c(left, right, up, down)`, i.e.
#'   first and last column, first and last row.
#' @param p Order of the Rényi distance.
#' @param L Known number of looks; estimated when `NULL`.
#'
#' @return An object of class `distancesg0itable` with the marks and the
#'   table of distances (`distancestable`, channels in rows).
#'
#' @examples
#' \donttest{
#' res <- distancesg0itable(SanFrancisco150, c(10, 30, 10, 30),
#'                          c(10, 30, 120, 140), L = 4)
#' res
#' }
#'
#' @export
distancesg0itable <- function(image, marks1, marks2, p = 1 / 2, L = NULL) {
  if (length(marks1) != 4 || length(marks2) != 4 || any(marks1 < 1) || any(marks2 < 1)) {
    stop("`marks1` and `marks2` must be c(left, right, up, down) with positive indices.",
         call. = FALSE)
  }
  img <- .as_sarimage(image)
  I <- img$intensities()
  types <- c("KL", "Renyi", "Bhattacharyya", "Hellinger",
             "Arithmetic-Geometric", "Triangular", "Harmonic-Mean")
  region <- function(mk, k) as.vector(I[mk[3]:mk[4], mk[1]:mk[2], k])
  res <- t(vapply(seq_len(img$m), function(k) {
    e1 <- G0I$new()$fit(region(marks1, k), L = L)$params
    e2 <- G0I$new()$fit(region(marks2, k), L = L)$params
    vapply(types, function(tp) {
      distancesg0i(unlist(e1), unlist(e2), p = p, type = tp)$distance_value
    }, numeric(1))
  }, numeric(length(types))))
  dimnames(res) <- list(img$channels, types)
  structure(list(marks1 = marks1, marks2 = marks2, distancestable = res),
            class = "distancesg0itable")
}

#' @export
print.distancesg0itable <- function(x, ...) {
  cat("\nG0I distances table\n\n")
  print(x$distancestable)
  invisible(x)
}

#' Estimation of G0I parameters
#'
#' Maximum likelihood estimation of the \eqn{\mathcal{G}_I^0} parameters.
#' This is a convenience wrapper around [G0I]`$fit()`, kept for
#' compatibility with previous versions of the package.
#'
#' @param Theta Initial values `c(alpha, gamma, L)`. When `L_fixed = TRUE`,
#'   `Theta[3]` is the known number of looks.
#' @param x Vector of observations.
#' @param L_fixed Logical; if `TRUE`, \eqn{L} is held fixed at `Theta[3]`.
#'
#' @return Named vector `c(alpha, gamma, L)` of estimates.
#'
#' @examples
#' set.seed(1)
#' y <- rg0i(1000, alpha = -1.5, gamma = 2, L = 8)
#' g0iestimation(c(-1, 1, 1), y)
#' g0iestimation(c(-1, 1, 8), y, L_fixed = TRUE)
#'
#' @export
g0iestimation <- function(Theta, x, L_fixed = FALSE) {
  start <- c(alpha = Theta[1], gamma = Theta[2], L = Theta[3])
  if (!(start[["alpha"]] < 0)) start[["alpha"]] <- -1.5
  if (!(start[["gamma"]] > 0)) start[["gamma"]] <- 1
  fit <- G0I$new()$fit(x, L = if (L_fixed) Theta[3] else NULL, start = start)
  unlist(fit$params)
}
