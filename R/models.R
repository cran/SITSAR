#' Base R6 class for SAR and PolSAR statistical models
#'
#' `SARModel` is the abstract parent of all model classes of `SITSAR`
#' ([GammaSAR], [G0I], [GGammaSAR], [WishartSAR], [G0Polarimetric] and
#' [CompoundWishart]). Every model can evaluate densities, simulate data,
#' be fitted by maximum likelihood (`$fit()`) and report its Shannon entropy
#' and the asymptotic variance of the entropy estimator.
#'
#' Fitted models work with the `generics` verbs [tidy()][generics::tidy],
#' [glance()][generics::glance], as well as with [stats::logLik()] and
#' [stats::AIC()].
#'
#' @field name Model name.
#' @field params Named list of parameters.
#' @field fixed Names of parameters held fixed during estimation.
#' @field n Sample size used in the last fit (`NULL` if not fitted).
#' @field loglik Maximised log-likelihood (`NULL` if not fitted).
#' @field converged Logical convergence flag of the last fit.
#'
#' @return An [R6::R6Class] generator. `SARModel$new()` returns a new `SARModel`
#'   object; its methods are described in the Methods section below.
#' @export
SARModel <- R6::R6Class(
  "SARModel",
  public = list(
    name = "SARModel",
    params = NULL,
    fixed = character(0),
    n = NULL,
    loglik = NULL,
    converged = NA,

    #' @description Abstract constructor.
    #' @param ... Ignored.
    initialize = function(...) {
      stop("`SARModel` is an abstract class.", call. = FALSE)
    },

    #' @description (Log-)density of observations.
    #' @param x Observations.
    #' @param log Return log-densities?
    density = function(x, log = FALSE) stop("Not implemented.", call. = FALSE),

    #' @description Generate random observations.
    #' @param n Number of observations.
    random = function(n) stop("Not implemented.", call. = FALSE),

    #' @description Fit the model by maximum likelihood.
    #' @param x Observations.
    #' @param ... Model specific arguments.
    fit = function(x, ...) stop("Not implemented.", call. = FALSE),

    #' @description Log-likelihood of a sample.
    #' @param x Observations.
    log_likelihood = function(x) sum(self$density(x, log = TRUE)),

    #' @description Shannon entropy at the current parameters.
    entropy = function() stop("Not implemented.", call. = FALSE),

    #' @description Asymptotic variance of the entropy estimator (delta
    #'   method with the inverse observed Fisher information).
    #' @param n Sample size. Defaults to the size of the fitted sample.
    entropy_variance = function(n = self$n) {
      if (is.null(private$.vcov)) {
        stop("Fit the model first with `$fit()`.", call. = FALSE)
      }
      free <- private$.free_names()
      theta <- unlist(self$params[free])
      g <- .num_grad(function(th) private$.entropy_at(th), theta)
      v <- as.numeric(t(g) %*% private$.vcov %*% g) + private$.extra_var
      v * self$n / n
    },

    #' @description Number of free parameters.
    n_parameters = function() length(private$.free_names()) + private$.extra_df,

    #' @description Parameters as a tibble.
    tidy = function() {
      p <- self$params
      scal <- vapply(p, function(v) is.numeric(v) && length(v) == 1L, logical(1))
      se <- rep(NA_real_, sum(scal))
      names(se) <- names(p)[scal]
      if (!is.null(private$.vcov)) {
        free <- intersect(private$.free_names(), names(se))
        se[free] <- sqrt(pmax(diag(private$.vcov)[free], 0))
      }
      tibble::tibble(
        term = names(p)[scal],
        estimate = unname(unlist(p[scal])),
        std.error = unname(se),
        fixed = names(p)[scal] %in% self$fixed
      )
    },

    #' @description One-row summary of the fit.
    glance = function() {
      ll <- self$loglik %||% NA_real_
      k <- self$n_parameters()
      tibble::tibble(
        model = self$name,
        n = self$n %||% NA_integer_,
        logLik = ll,
        AIC = -2 * ll + 2 * k,
        BIC = -2 * ll + log(self$n %||% NA_real_) * k,
        entropy = tryCatch(self$entropy(), error = function(e) NA_real_),
        converged = self$converged
      )
    },

    #' @description Print the model.
    #' @param ... Ignored.
    print = function(...) {
      cat("<", self$name, "> ", if (is.null(self$n)) "(not fitted)" else
        paste0("fitted to n = ", self$n, " observations"), "\n", sep = "")
      for (nm in names(self$params)) {
        v <- self$params[[nm]]
        if (is.matrix(v)) {
          cat("  ", nm, ":\n", sep = "")
          print(signif(v, 4))
        } else {
          cat("  ", nm, " = ", format(signif(v, 6)),
              if (nm %in% self$fixed) " (fixed)" else "", "\n", sep = "")
        }
      }
      if (!is.null(self$loglik)) cat("  logLik = ", format(signif(self$loglik, 6)), "\n", sep = "")
      ent <- tryCatch(self$entropy(), error = function(e) NULL)
      if (!is.null(ent)) cat("  Shannon entropy = ", format(signif(ent, 6)), "\n", sep = "")
      invisible(self)
    }
  ),
  private = list(
    .vcov = NULL,
    .extra_var = 0,
    .extra_df = 0,
    .free_names = function() {
      nm <- names(self$params)
      scal <- vapply(self$params, function(v) is.numeric(v) && length(v) == 1L, logical(1))
      setdiff(nm[scal], self$fixed)
    },
    .entropy_at = function(theta) {
      old <- self$params
      on.exit(self$params <- old)
      for (nm in names(theta)) self$params[[nm]] <- unname(theta[nm])
      self$entropy()
    },
    .set_vcov = function(negll_natural, theta) {
      H <- tryCatch(stats::optimHess(theta, negll_natural), error = function(e) NULL)
      private$.vcov <- if (is.null(H)) NULL else tryCatch(solve(H), error = function(e) NULL)
      if (!is.null(private$.vcov)) dimnames(private$.vcov) <- list(names(theta), names(theta))
      invisible(NULL)
    }
  )
)


# --------------------------------------------------------------------------
#' Gamma model for SAR intensities (R6)
#'
#' R6 class for the multilook Gamma law \eqn{\Gamma(\mu, L)}; see
#' [dgammasar()]. The maximum likelihood estimator of \eqn{\mu} is the
#' sample mean and \eqn{L} solves \eqn{\log L - \psi^{(0)}(L) =
#' \log\bar z - \overline{\log z}}.
#'
#' @examples
#' set.seed(1)
#' z <- rgammasar(500, mu = 2, L = 4)
#' m <- GammaSAR$new()$fit(z)
#' m
#' m$entropy()
#' m$entropy_variance()
#' generics::tidy(m)
#'
#' @return An [R6::R6Class] generator. `GammaSAR$new()` returns a new `GammaSAR`
#'   object; its methods are described in the Methods section below.
#' @export
GammaSAR <- R6::R6Class(
  "GammaSAR",
  inherit = SARModel,
  public = list(
    #' @description Create a Gamma model.
    #' @param mu Mean.
    #' @param L Number of looks.
    initialize = function(mu = 1, L = 1) {
      self$name <- "Gamma"
      self$params <- list(mu = mu, L = L)
    },
    #' @description (Log-)density.
    #' @param x Numeric vector of intensities.
    #' @param log Return log-densities?
    density = function(x, log = FALSE) dgammasar(x, self$params$mu, self$params$L, log = log),
    #' @description Simulate intensities.
    #' @param n Number of observations.
    random = function(n) rgammasar(n, self$params$mu, self$params$L),
    #' @description Fit by maximum likelihood.
    #' @param x Numeric vector of positive intensities.
    #' @param L Known number of looks. If `NULL` (default) it is estimated.
    fit = function(x, L = NULL) {
      x <- as.numeric(Re(x))
      x <- x[is.finite(x) & x > 0]
      if (length(x) < 2) stop("At least two positive observations are needed.", call. = FALSE)
      mu <- mean(x)
      if (is.null(L)) {
        s <- log(mu) - mean(log(x))
        L <- .vec_bisect_log(function(l) log(l) - digamma(l), s, 1e-4, 1e7)
        self$fixed <- character(0)
      } else {
        self$fixed <- "L"
      }
      self$params <- list(mu = mu, L = L)
      self$n <- length(x)
      self$loglik <- sum(self$density(x, log = TRUE))
      self$converged <- is.finite(self$loglik)
      private$.vcov <- diag(c(mu^2 / (self$n * L),
                              if (is.null(self$fixed) || !"L" %in% self$fixed)
                                1 / (self$n * (trigamma(L) - 1 / L))),
                            nrow = 2 - length(self$fixed))
      dimnames(private$.vcov) <- list(private$.free_names(), private$.free_names())
      invisible(self)
    },
    #' @description Shannon entropy.
    entropy = function() entropy_gamma(self$params$mu, self$params$L),
    #' @description Asymptotic variance of the entropy estimator (closed form).
    #' @param n Sample size.
    entropy_variance = function(n = self$n) {
      if (is.null(n)) stop("Fit the model first or supply `n`.", call. = FALSE)
      .entropy_var_wishart(self$params$L, 1, n, !"L" %in% self$fixed)
    }
  )
)


# --------------------------------------------------------------------------
#' G0I model for SAR intensities (R6)
#'
#' R6 class for the \eqn{\mathcal{G}_I^0(\alpha, \gamma, L)} law (see
#' [dg0i()]) with maximum likelihood estimation (analytic score, quasi-Newton
#' optimisation on a transformed scale with several starting points),
#' closed-form Shannon entropy (Ferreira and Nascimento, 2020) and the
#' asymptotic variance of its estimator.
#'
#' @references
#' Ferreira, J. A. and Nascimento, A. D. C. (2020). Shannon entropy for the
#' G0I model: a new segmentation approach. *IEEE Journal of Selected Topics
#' in Applied Earth Observations and Remote Sensing*, 13, 2547--2553.
#' \doi{10.1109/JSTARS.2020.2997666}
#'
#' @examples
#' set.seed(1)
#' z <- rg0i(1000, alpha = -3, gamma = 2, L = 4)
#' m <- G0I$new()$fit(z, L = 4)
#' m
#' m$entropy_variance()
#' G0I$new()$fit(z)$params # L estimated
#'
#' @return An [R6::R6Class] generator. `G0I$new()` returns a new `G0I`
#'   object; its methods are described in the Methods section below.
#' @export
G0I <- R6::R6Class(
  "G0I",
  inherit = SARModel,
  public = list(
    #' @description Create a G0I model.
    #' @param alpha Roughness (< 0).
    #' @param gamma Scale (> 0).
    #' @param L Number of looks (> 0).
    initialize = function(alpha = -2, gamma = 1, L = 1) {
      self$name <- "G0I"
      self$params <- list(alpha = alpha, gamma = gamma, L = L)
    },
    #' @description (Log-)density.
    #' @param x Numeric vector of intensities.
    #' @param log Return log-densities?
    density = function(x, log = FALSE) {
      p <- self$params
      dg0i(x, p$alpha, p$gamma, p$L, log = log)
    },
    #' @description Simulate intensities.
    #' @param n Number of observations.
    random = function(n) {
      p <- self$params
      rg0i(n, p$alpha, p$gamma, p$L)
    },
    #' @description Fit by maximum likelihood.
    #' @param x Numeric vector of positive intensities.
    #' @param L Known number of looks. If `NULL` (default) it is estimated.
    #' @param start Optional named vector of starting values
    #'   (`alpha`, `gamma` and, if estimated, `L`).
    fit = function(x, L = NULL, start = NULL) {
      x <- as.numeric(Re(x))
      x <- x[is.finite(x) & x > 0]
      if (length(x) < 3) stop("At least three positive observations are needed.", call. = FALSE)
      est <- .fit_g0i(x, L = L, start = start)
      self$params <- list(alpha = est$par[["alpha"]], gamma = est$par[["gamma"]],
                          L = est$par[["L"]])
      self$fixed <- if (is.null(L)) character(0) else "L"
      self$n <- length(x)
      self$loglik <- est$loglik
      self$converged <- est$converged
      free <- private$.free_names()
      fixedL <- self$params$L
      negll <- function(th) {
        a <- th[["alpha"]]; g <- th[["gamma"]]
        l <- if ("L" %in% names(th)) th[["L"]] else fixedL
        if (a >= 0 || g <= 0 || l <= 0) return(1e100)
        -sum(dg0i(x, a, g, l, log = TRUE))
      }
      private$.set_vcov(negll, unlist(self$params[free]))
      invisible(self)
    },
    #' @description Shannon entropy.
    entropy = function() {
      p <- self$params
      entropy_g0i(p$alpha, p$gamma, p$L)
    }
  )
)

# G0I log-likelihood score (natural parameters)
.g0i_score <- function(alpha, gamma, L, x) {
  n <- length(x)
  lg <- log(gamma + L * x)
  c(alpha = -n * log(gamma) - n * digamma(L - alpha) + n * digamma(-alpha) + sum(lg),
    gamma = -n * alpha / gamma + (alpha - L) * sum(1 / (gamma + L * x)),
    L = n * (1 + log(L)) + n * digamma(L - alpha) - n * digamma(L) +
      sum(log(x)) - sum(lg) + (alpha - L) * sum(x / (gamma + L * x)))
}

.fit_g0i <- function(x, L = NULL, start = NULL) {
  estimate_L <- is.null(L)
  m1 <- mean(x)
  if (estimate_L) {
    s <- log(m1) - mean(log(x))
    L0 <- .vec_bisect_log(function(l) log(l) - digamma(l), s, 1e-3, 1e6)
    L0 <- min(max(L0, 0.5), 50)
  } else {
    .check_scalar(L, "L", lower = 0)
    L0 <- L
  }
  # moment based starting value for alpha
  r <- mean(x^2) / m1^2
  den <- r * L0 - L0 - 1
  a_mom <- if (den > 0) -(2 * r * L0 - L0 - 1) / den else -20
  a_mom <- min(max(a_mom, -50), -1.1)
  alphas <- unique(c(start[["alpha"]], a_mom, -1.5, -3, -8, -20))
  negll <- function(th) {
    a <- -exp(th[1]); g <- exp(th[2]); l <- if (estimate_L) exp(th[3]) else L
    -sum(dg0i(x, a, g, l, log = TRUE))
  }
  grad <- function(th) {
    a <- -exp(th[1]); g <- exp(th[2]); l <- if (estimate_L) exp(th[3]) else L
    sc <- .g0i_score(a, g, l, x)
    out <- -c(sc[1] * a, sc[2] * g)
    if (estimate_L) out <- c(out, -sc[3] * l)
    out
  }
  lower <- c(log(1e-3), -Inf); upper <- c(log(1e3), Inf)
  if (estimate_L) { lower <- c(lower, log(0.05)); upper <- c(upper, log(500)) }
  best <- NULL
  for (a0 in alphas) {
    g0 <- if (!is.null(start[["gamma"]]) && a0 == start[["alpha"]]) start[["gamma"]] else m1 * (-a0 - 1)
    g0 <- max(g0, 1e-8 * m1)
    th0 <- c(log(-a0), log(g0))
    if (estimate_L) th0 <- c(th0, log(start[["L"]] %||% L0))
    op <- tryCatch(
      suppressWarnings(stats::optim(th0, negll, grad, method = "L-BFGS-B",
                                    lower = lower, upper = upper,
                                    control = list(maxit = 500))),
      error = function(e) NULL)
    if (!is.null(op) && is.finite(op$value) && (is.null(best) || op$value < best$value)) {
      best <- op
    }
  }
  if (is.null(best)) stop("G0I estimation failed.", call. = FALSE)
  th <- best$par
  par <- c(alpha = -exp(th[1]), gamma = exp(th[2]),
           L = if (estimate_L) exp(th[3]) else L)
  list(par = par, loglik = -best$value, converged = best$convergence == 0)
}


# --------------------------------------------------------------------------
#' Generalized Gamma model for SAR intensities (R6)
#'
#' R6 class for the generalized Gamma law (see [dggamma()]). Estimation is
#' carried out by profiling: for a given power \eqn{\alpha}, the maximum
#' likelihood estimates of \eqn{\phi} and \eqn{\mu} are available in closed
#' form (up to a one-dimensional root), so the profile likelihood is
#' maximised over \eqn{\alpha} only.
#'
#' @examples
#' set.seed(1)
#' z <- rggamma(1000, alpha = 1.5, phi = 3, mu = 2)
#' m <- GGammaSAR$new()$fit(z)
#' m
#'
#' @return An [R6::R6Class] generator. `GGammaSAR$new()` returns a new `GGammaSAR`
#'   object; its methods are described in the Methods section below.
#' @export
GGammaSAR <- R6::R6Class(
  "GGammaSAR",
  inherit = SARModel,
  public = list(
    #' @description Create a generalized Gamma model.
    #' @param alpha Power (> 0).
    #' @param phi Shape (> 0).
    #' @param mu Inverse scale (> 0).
    initialize = function(alpha = 1, phi = 1, mu = 1) {
      self$name <- "Generalized Gamma"
      self$params <- list(alpha = alpha, phi = phi, mu = mu)
    },
    #' @description (Log-)density.
    #' @param x Numeric vector of intensities.
    #' @param log Return log-densities?
    density = function(x, log = FALSE) {
      p <- self$params
      dggamma(x, p$alpha, p$phi, p$mu, log = log)
    },
    #' @description Simulate intensities.
    #' @param n Number of observations.
    random = function(n) {
      p <- self$params
      rggamma(n, p$alpha, p$phi, p$mu)
    },
    #' @description Fit by (profile) maximum likelihood.
    #' @param x Numeric vector of positive intensities.
    #' @param phi Known shape (number of looks). If `NULL` it is estimated.
    fit = function(x, phi = NULL) {
      x <- as.numeric(Re(x))
      x <- x[is.finite(x) & x > 0]
      if (length(x) < 3) stop("At least three positive observations are needed.", call. = FALSE)
      lx <- log(x)
      # rescale for numerical stability
      sc <- exp(mean(lx))
      xs <- x / sc
      lxs <- lx - log(sc)
      prof <- function(la) {
        a <- exp(la)
        xa <- xs^a
        mxa <- mean(xa)
        ph <- if (is.null(phi)) {
          s <- log(mxa) - a * mean(lxs)
          .vec_bisect_log(function(f) log(f) - digamma(f), s, 1e-4, 1e7)
        } else phi
        mu <- (ph / mxa)^(1 / a)
        ll <- sum(dggamma(xs, a, ph, mu, log = TRUE))
        attr(ll, "par") <- c(alpha = a, phi = ph, mu = mu / sc)
        ll
      }
      op <- stats::optimize(function(la) -prof(la), c(log(0.02), log(30)))
      ll <- prof(op$minimum)
      par <- attr(ll, "par")
      self$params <- list(alpha = par[["alpha"]], phi = par[["phi"]], mu = par[["mu"]])
      self$fixed <- if (is.null(phi)) character(0) else "phi"
      self$n <- length(x)
      self$loglik <- sum(self$density(x, log = TRUE))
      self$converged <- is.finite(self$loglik)
      free <- private$.free_names()
      fixed_phi <- self$params$phi
      negll <- function(th) {
        a <- th[["alpha"]]; m <- th[["mu"]]
        f <- if ("phi" %in% names(th)) th[["phi"]] else fixed_phi
        if (a <= 0 || f <= 0 || m <= 0) return(1e100)
        -sum(dggamma(x, a, f, m, log = TRUE))
      }
      private$.set_vcov(negll, unlist(self$params[free]))
      invisible(self)
    },
    #' @description Shannon entropy.
    entropy = function() {
      p <- self$params
      entropy_ggamma(p$alpha, p$phi, p$mu)
    }
  )
)


# --------------------------------------------------------------------------
#' Scaled complex Wishart model for PolSAR data (R6)
#'
#' R6 class for the scaled complex Wishart law \eqn{\mathcal{W}_m(\Sigma, L)}
#' (see [dwishart_sar()]). The maximum likelihood estimator of \eqn{\Sigma}
#' is the sample mean and \eqn{L} solves
#' \eqn{m\log L - \psi^{(0)}_m(L) = \log|\bar Z| - \overline{\log|Z|}}
#' (Anfinsen et al., 2009). The entropy variance is available in closed
#' form (Nascimento, Ferreira and Frery, 2023).
#'
#' @references
#' Anfinsen, S. N., Doulgeris, A. P. and Eltoft, T. (2009). Estimation of
#' the equivalent number of looks in polarimetric synthetic aperture radar
#' imagery. *IEEE Transactions on Geoscience and Remote Sensing*, 47(11),
#' 3795--3809.
#'
#' @examples
#' set.seed(1)
#' Z <- rwishart_sar(300, polsar_example_sigma("forest"), L = 4)
#' m <- WishartSAR$new()$fit(Z)
#' m
#' m$entropy_variance()
#'
#' @return An [R6::R6Class] generator. `WishartSAR$new()` returns a new `WishartSAR`
#'   object; its methods are described in the Methods section below.
#' @export
WishartSAR <- R6::R6Class(
  "WishartSAR",
  inherit = SARModel,
  public = list(
    #' @description Create a Wishart model.
    #' @param Sigma Hermitian positive definite matrix (default identity 3x3).
    #' @param L Number of looks.
    initialize = function(Sigma = diag(3) + 0i, L = 4) {
      self$name <- "Scaled complex Wishart"
      self$params <- list(Sigma = Sigma, L = L)
    },
    #' @description (Log-)density.
    #' @param x PolSAR observations (array `m x m x n` or vec matrix).
    #' @param log Return log-densities?
    density = function(x, log = FALSE) dwishart_sar(x, self$params$Sigma, self$params$L, log = log),
    #' @description Simulate PolSAR observations.
    #' @param n Number of observations.
    #' @param output `"array"` or `"vec"`.
    random = function(n, output = "array") rwishart_sar(n, self$params$Sigma, self$params$L, output),
    #' @description Fit by maximum likelihood.
    #' @param x PolSAR observations (array `m x m x n` or vec matrix).
    #' @param L Known number of looks. If `NULL` (default) it is estimated.
    fit = function(x, L = NULL) {
      V <- .as_vec(x)
      m <- .m_from_p(ncol(V))
      Sigma <- .vec_to_herm(colMeans(V))
      if (is.null(L)) {
        s <- .hlogdet(Sigma) - mean(log(.hdet_vec(V, m)))
        L <- .vec_bisect_log(function(l) m * log(l) - .mvpsigamma(l, m), s,
                             m - 1 + 1e-6, 1e7)
        self$fixed <- character(0)
      } else {
        .check_scalar(L, "L", lower = m - 1)
        self$fixed <- "L"
      }
      self$params <- list(Sigma = Sigma, L = L)
      self$n <- nrow(V)
      self$loglik <- sum(self$density(V, log = TRUE))
      self$converged <- is.finite(self$loglik)
      private$.vcov <- if (!"L" %in% self$fixed) {
        matrix(1 / (self$n * (.mvpsigamma(L, m, 1L) - m / L)), 1, 1,
               dimnames = list("L", "L"))
      } else matrix(numeric(0), 0, 0)
      private$.extra_df <- m^2
      invisible(self)
    },
    #' @description Shannon entropy.
    entropy = function() entropy_wishart(self$params$Sigma, self$params$L),
    #' @description Asymptotic variance of the entropy estimator (closed form).
    #' @param n Sample size.
    entropy_variance = function(n = self$n) {
      if (is.null(n)) stop("Fit the model first or supply `n`.", call. = FALSE)
      .entropy_var_wishart(self$params$L, nrow(self$params$Sigma), n,
                           !"L" %in% self$fixed)
    }
  )
)


# --------------------------------------------------------------------------
#' G0 polarimetric model for PolSAR data (R6)
#'
#' R6 class for the \eqn{\mathcal{G}^0_m(\alpha, \Sigma, L)} law (see
#' [dg0p()]). Following Nascimento, Ferreira and Frery (2023),
#' \eqn{\Sigma} is estimated by the sample mean (since \eqn{E(Z) = \Sigma})
#' and \eqn{(\alpha, L)} by maximising the likelihood (analytic score).
#' The variance of the entropy estimator combines the inverse observed
#' information of \eqn{(\alpha, L)} with a moment estimate of the variance
#' of \eqn{m\log|\hat\Sigma|}.
#'
#' @examples
#' set.seed(1)
#' Z <- rg0p(400, alpha = -4, Sigma = polsar_example_sigma("urban"), L = 4)
#' m <- G0Polarimetric$new()$fit(Z, L = 4)
#' m
#' m$entropy_variance()
#'
#' @return An [R6::R6Class] generator. `G0Polarimetric$new()` returns a new `G0Polarimetric`
#'   object; its methods are described in the Methods section below.
#' @export
G0Polarimetric <- R6::R6Class(
  "G0Polarimetric",
  inherit = SARModel,
  public = list(
    #' @description Create a G0 polarimetric model.
    #' @param alpha Roughness (< -1).
    #' @param Sigma Hermitian positive definite matrix.
    #' @param L Number of looks.
    initialize = function(alpha = -3, Sigma = diag(3) + 0i, L = 4) {
      self$name <- "G0 polarimetric"
      self$params <- list(alpha = alpha, Sigma = Sigma, L = L)
    },
    #' @description (Log-)density.
    #' @param x PolSAR observations.
    #' @param log Return log-densities?
    density = function(x, log = FALSE) {
      p <- self$params
      dg0p(x, p$alpha, p$Sigma, p$L, log = log)
    },
    #' @description Simulate PolSAR observations.
    #' @param n Number of observations.
    #' @param output `"array"` or `"vec"`.
    random = function(n, output = "array") {
      p <- self$params
      rg0p(n, p$alpha, p$Sigma, p$L, output)
    },
    #' @description Fit by maximum likelihood.
    #' @param x PolSAR observations.
    #' @param L Known number of looks. If `NULL` (default) it is estimated.
    fit = function(x, L = NULL) {
      V <- .as_vec(x)
      m <- .m_from_p(ncol(V))
      Sigma <- .vec_to_herm(colMeans(V))
      est <- .fit_g0p(V, Sigma, m, L)
      self$params <- list(alpha = est$par[["alpha"]], Sigma = Sigma, L = est$par[["L"]])
      self$fixed <- if (is.null(L)) character(0) else "L"
      self$n <- nrow(V)
      self$loglik <- est$loglik
      self$converged <- est$converged
      ld <- .hlogdet(Sigma)
      ldz <- log(.hdet_vec(V, m))
      tr <- .htrace_vec(.hinv(Sigma), V)
      fixedL <- self$params$L
      negll <- function(th) {
        a <- th[["alpha"]]
        l <- if ("L" %in% names(th)) th[["L"]] else fixedL
        if (a >= -1 || l <= m - 1) return(1e100)
        -.g0p_loglik(a, l, ld, ldz, tr, m)
      }
      free <- private$.free_names()
      private$.set_vcov(negll, unlist(self$params[free]))
      private$.extra_var <- m^2 * stats::var(tr) / self$n
      private$.extra_df <- m^2
      invisible(self)
    },
    #' @description Shannon entropy.
    entropy = function() {
      p <- self$params
      entropy_g0p(p$alpha, p$Sigma, p$L)
    }
  )
)

.g0p_loglik <- function(alpha, L, ld, ldz, tr, m) {
  n <- length(ldz)
  n * (m * L * log(L) + lgamma(m * L - alpha) - L * ld - .lmvgamma(L, m) -
         lgamma(-alpha) - alpha * log(-alpha - 1)) +
    (L - m) * sum(ldz) + (alpha - m * L) * sum(log(L * tr - alpha - 1))
}

.g0p_score <- function(alpha, L, ld, ldz, tr, m) {
  n <- length(ldz)
  u <- L * tr - alpha - 1
  c(alpha = n * (-digamma(m * L - alpha) + digamma(-alpha) - log(-alpha - 1) +
                   alpha / (-alpha - 1)) + sum(log(u)) - (alpha - m * L) * sum(1 / u),
    L = n * (m * log(L) + m + m * digamma(m * L - alpha) - ld - .mvpsigamma(L, m)) +
      sum(ldz) - m * sum(log(u)) + (alpha - m * L) * sum(tr / u))
}

.fit_g0p <- function(V, Sigma, m, L = NULL) {
  estimate_L <- is.null(L)
  if (!estimate_L) .check_scalar(L, "L", lower = m - 1)
  ld <- .hlogdet(Sigma)
  ldz <- log(.hdet_vec(V, m))
  tr <- .htrace_vec(.hinv(Sigma), V)
  if (estimate_L) {
    s <- ld - mean(ldz)
    L0 <- .vec_bisect_log(function(l) m * log(l) - .mvpsigamma(l, m), s, m - 1 + 1e-6, 1e7)
    L0 <- min(max(L0, m - 1 + 0.5), 100)
  }
  to_nat <- function(th) {
    c(alpha = -1 - exp(th[1]), L = if (estimate_L) m - 1 + exp(th[2]) else L)
  }
  negll <- function(th) {
    p <- to_nat(th)
    -.g0p_loglik(p[["alpha"]], p[["L"]], ld, ldz, tr, m)
  }
  grad <- function(th) {
    p <- to_nat(th)
    sc <- .g0p_score(p[["alpha"]], p[["L"]], ld, ldz, tr, m)
    out <- -sc[1] * (p[["alpha"]] + 1)
    if (estimate_L) out <- c(out, -sc[2] * (p[["L"]] - m + 1))
    out
  }
  lower <- log(1e-3); upper <- log(1e3)
  if (estimate_L) { lower <- c(lower, log(0.05)); upper <- c(upper, log(500)) }
  best <- NULL
  for (a0 in c(-1.5, -3, -8, -30)) {
    th0 <- log(-a0 - 1)
    if (estimate_L) th0 <- c(th0, log(L0 - m + 1))
    op <- tryCatch(
      suppressWarnings(stats::optim(th0, negll, grad, method = "L-BFGS-B",
                                    lower = lower, upper = upper)),
      error = function(e) NULL)
    if (!is.null(op) && is.finite(op$value) && (is.null(best) || op$value < best$value)) best <- op
  }
  if (is.null(best)) stop("G0 polarimetric estimation failed.", call. = FALSE)
  list(par = to_nat(best$par), loglik = -best$value, converged = best$convergence == 0)
}


# --------------------------------------------------------------------------
#' Compound complex Wishart models (R6)
#'
#' R6 class for the compound truncated Poisson (CTPCW) and compound
#' geometric (CGCW) complex Wishart laws of Ferreira, Nascimento and Frery
#' (2022); see [dctpcw()]. Parameters are estimated by the
#' expectation-maximisation (EM) algorithm, treating the number of
#' scatterers \eqn{N} as missing data:
#' \itemize{
#'   \item E-step: \eqn{\hat N_i = E(N \mid S_i)};
#'   \item M-step: \eqn{\hat\Sigma = \bar S / (L \bar N)} and
#'     \eqn{\hat\lambda/(1 - e^{-\hat\lambda}) = \bar N} (CTPCW) or
#'     \eqn{\hat p = 1/\bar N} (CGCW).
#' }
#' The Shannon entropy has no closed form and is approximated by Monte
#' Carlo.
#'
#' @field counting Counting law: `"poisson"` (zero-truncated) or `"geometric"`.
#' @field iterations Number of EM iterations of the last fit.
#'
#' @references
#' Ferreira, J. A., Nascimento, A. D. C. and Frery, A. C. (2022). PolSAR
#' models with multimodal intensities. *Remote Sensing*, 14(20), 5083.
#' \doi{10.3390/rs14205083}
#'
#' @examples
#' set.seed(1)
#' S <- rctpcw(200, lambda = 1.5, Sigma = polsar_example_sigma("forest"), L = 4)
#' m <- CompoundWishart$new(counting = "poisson")$fit(S, L = 4)
#' m$params$lambda
#' m$iterations
#'
#' @return An [R6::R6Class] generator. `CompoundWishart$new()` returns a new `CompoundWishart`
#'   object; its methods are described in the Methods section below.
#' @export
CompoundWishart <- R6::R6Class(
  "CompoundWishart",
  inherit = SARModel,
  public = list(
    counting = "poisson",
    iterations = NA_integer_,
    #' @description Create a compound Wishart model.
    #' @param counting `"poisson"` (CTPCW) or `"geometric"` (CGCW).
    #' @param Sigma Hermitian positive definite matrix.
    #' @param L Number of looks of each component.
    #' @param lambda,p Parameter of the counting law.
    initialize = function(counting = c("poisson", "geometric"),
                          Sigma = diag(3) + 0i, L = 4, lambda = 1, p = 0.5) {
      self$counting <- match.arg(counting)
      self$name <- if (self$counting == "poisson") {
        "Compound truncated Poisson complex Wishart"
      } else "Compound geometric complex Wishart"
      self$params <- if (self$counting == "poisson") {
        list(lambda = lambda, Sigma = Sigma, L = L)
      } else list(p = p, Sigma = Sigma, L = L)
      self$fixed <- "L"
    },
    #' @description (Log-)density.
    #' @param x PolSAR observations.
    #' @param log Return log-densities?
    density = function(x, log = FALSE) {
      p <- self$params
      if (self$counting == "poisson") dctpcw(x, p$lambda, p$Sigma, p$L, log = log)
      else dcgcw(x, p$p, p$Sigma, p$L, log = log)
    },
    #' @description Simulate PolSAR observations.
    #' @param n Number of observations.
    #' @param output `"array"` or `"vec"`.
    random = function(n, output = "array") {
      p <- self$params
      if (self$counting == "poisson") rctpcw(n, p$lambda, p$Sigma, p$L, output)
      else rcgcw(n, p$p, p$Sigma, p$L, output)
    },
    #' @description Fit by the EM algorithm.
    #' @param x PolSAR observations.
    #' @param L Known number of looks of each component.
    #' @param tol Convergence tolerance (relative change of the parameters).
    #' @param max_iter Maximum number of iterations.
    fit = function(x, L, tol = 1e-6, max_iter = 500L) {
      V <- .as_vec(x)
      m <- .m_from_p(ncol(V))
      .check_scalar(L, "L", lower = m - 1)
      Sbar <- colMeans(V)
      Nbar <- 2
      par <- .count_par_from_mean(Nbar, self$counting)
      Sigma <- .vec_to_herm(Sbar / (L * Nbar))
      it <- 0L
      conv <- FALSE
      repeat {
        it <- it + 1L
        tt <- .compound_terms(V, self$counting, par, Sigma, L)
        K <- ncol(tt$terms)
        w <- exp(tt$terms - apply(tt$terms, 1, max))
        EN <- as.vector(w %*% seq_len(K)) / rowSums(w)
        Nbar_new <- mean(EN)
        par_new <- .count_par_from_mean(Nbar_new, self$counting)
        Sigma_new <- .vec_to_herm(Sbar / (L * Nbar_new))
        delta <- abs(par_new - par) / max(abs(par), 1e-8) +
          max(Mod(Sigma_new - Sigma)) / max(Mod(Sigma))
        par <- par_new
        Sigma <- Sigma_new
        if (delta < tol) { conv <- TRUE; break }
        if (it >= max_iter) break
      }
      self$params <- if (self$counting == "poisson") {
        list(lambda = par, Sigma = Sigma, L = L)
      } else list(p = par, Sigma = Sigma, L = L)
      self$iterations <- it
      self$converged <- conv
      self$n <- nrow(V)
      self$loglik <- sum(self$density(V, log = TRUE))
      private$.extra_df <- m^2
      private$.vcov <- NULL
      invisible(self)
    },
    #' @description Monte Carlo approximation of the Shannon entropy.
    #' @param nsim Number of simulated observations.
    entropy = function(nsim = 2000L) {
      -mean(self$density(self$random(nsim, output = "vec"), log = TRUE))
    },
    #' @description Not available for compound models.
    #' @param n Ignored.
    entropy_variance = function(n = self$n) {
      stop("No closed-form entropy variance for compound Wishart models.", call. = FALSE)
    }
  )
)

.count_par_from_mean <- function(Nbar, counting) {
  if (counting == "geometric") return(min(max(1 / Nbar, 1e-8), 1 - 1e-8))
  if (Nbar <= 1 + 1e-8) return(1e-6)
  f <- function(l) l / (-expm1(-l)) - Nbar
  stats::uniroot(f, c(1e-8, max(2 * Nbar, 10)), extendInt = "upX", tol = 1e-10)$root
}


# S3 bridges ---------------------------------------------------------------

#' Tidy and glance methods for SITSAR models
#'
#' @param x A fitted [SARModel] object.
#' @param object A fitted [SARModel] object.
#' @param ... Not used.
#'
#' @return `tidy()` returns a tibble with one row per scalar parameter;
#'   `glance()` a one-row tibble; `logLik()` an object of class `logLik`.
#'
#' @examples
#' set.seed(1)
#' m <- G0I$new()$fit(rg0i(300, -2, 1, 4), L = 4)
#' generics::tidy(m)
#' generics::glance(m)
#' logLik(m)
#'
#' @name sarmodel-methods
NULL

#' @rdname sarmodel-methods
#' @importFrom generics tidy
#' @export
tidy.SARModel <- function(x, ...) x$tidy()

#' @rdname sarmodel-methods
#' @importFrom generics glance
#' @export
glance.SARModel <- function(x, ...) x$glance()

#' @rdname sarmodel-methods
#' @export
logLik.SARModel <- function(object, ...) {
  if (is.null(object$loglik)) stop("The model has not been fitted.", call. = FALSE)
  structure(object$loglik, df = object$n_parameters(), nobs = object$n,
            class = "logLik")
}

#' @export
generics::tidy

#' @export
generics::glance
