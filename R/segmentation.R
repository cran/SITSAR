#' Entropy-based segmentation by (stochastic) EM (R6)
#'
#' Unsupervised segmentation of SAR/PolSAR images from maps of estimated
#' Shannon entropies (Ferreira and Nascimento, 2020; Nascimento, Ferreira and
#' Frery, 2023). By the asymptotic normality of the entropy estimators, the
#' entropy vector of a pixel of the \eqn{g}-th region is approximately
#' Gaussian, so the image is modelled as a finite mixture of Gaussian laws
#' (with diagonal covariance matrices) fitted by the expectation-maximisation
#' algorithm (`method = "em"`) or by its stochastic version (`method =
#' "sem"`), in which the E-step is followed by a random allocation of every
#' pixel to a class drawn from its posterior probabilities. Pixels are then
#' labelled by the maximum a posteriori rule.
#'
#' @field k Number of classes.
#' @field method `"em"` or `"sem"`.
#' @field weights Mixture weights.
#' @field means Matrix (`k x features`) of class means.
#' @field variances Matrix (`k x features`) of class variances.
#' @field posterior Matrix of posterior probabilities.
#' @field labels Integer vector of labels (classes ordered by the mean of
#'   the first feature).
#' @field loglik Trace of the log-likelihood.
#' @field iterations Number of iterations.
#' @field converged Logical convergence flag.
#' @field dim Spatial dimensions when fitted to a map (`NULL` otherwise).
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
#' H <- entropy_map(SanFrancisco150, model = "wishart", window = 7, L = 4)
#' set.seed(1)
#' seg <- EntropySegmenter$new(k = 3)$fit(H)
#' seg
#' plot_map(seg$label_map(), main = "Entropy-based segmentation")
#' kappa_statistics(SanFrancisco150_gt, seg$label_map())
#'
#' @return An [R6::R6Class] generator. `EntropySegmenter$new()` returns a new `EntropySegmenter`
#'   object; its methods are described in the Methods section below.
#' @export
EntropySegmenter <- R6::R6Class(
  "EntropySegmenter",
  public = list(
    k = NULL,
    method = NULL,
    weights = NULL,
    means = NULL,
    variances = NULL,
    posterior = NULL,
    labels = NULL,
    loglik = NULL,
    iterations = NULL,
    converged = NULL,
    dim = NULL,

    #' @description Create a segmenter.
    #' @param k Number of classes.
    #' @param method `"em"` or `"sem"`.
    #' @param max_iter Maximum number of iterations.
    #' @param tol Tolerance for the relative change of the log-likelihood.
    initialize = function(k = 2, method = c("em", "sem"), max_iter = 200L, tol = 1e-8) {
      if (!is.numeric(k) || length(k) != 1 || k < 1) stop("`k` must be a positive integer.", call. = FALSE)
      self$k <- as.integer(k)
      self$method <- match.arg(method)
      private$.max_iter <- as.integer(max_iter)
      private$.tol <- tol
    },

    #' @description Fit the mixture.
    #' @param x Entropy features: a matrix (one map, or pixels in rows and
    #'   features in columns when `is_map = FALSE`), an array
    #'   `nrow x ncol x features` (several maps), or a numeric vector.
    #' @param init Initial labels (integer vector) or `NULL` for k-means.
    #' @param is_map Treat a matrix `x` as a single map? Default `TRUE`.
    fit = function(x, init = NULL, is_map = TRUE) {
      X <- private$.as_features(x, is_map)
      ok <- stats::complete.cases(X)
      Xo <- X[ok, , drop = FALSE]
      n <- nrow(Xo); p <- ncol(Xo); k <- self$k
      if (n < k) stop("Not enough complete observations.", call. = FALSE)
      if (is.null(init)) {
        Z <- scale(Xo)
        Z[!is.finite(Z)] <- 0
        init <- stats::kmeans(Z, centers = k, nstart = 5, iter.max = 50)$cluster
      } else {
        init <- as.integer(as.factor(init))[ok]
      }
      st <- private$.mstep(Xo, init)
      ll_trace <- numeric(0)
      conv <- FALSE
      for (it in seq_len(private$.max_iter)) {
        e <- private$.estep(Xo, st)
        ll_trace <- c(ll_trace, e$loglik)
        if (self$method == "em") {
          st <- private$.mstep_soft(Xo, e$post)
        } else {
          cls <- apply(e$post, 1, function(pr) sample.int(k, 1, prob = pr))
          st <- private$.mstep(Xo, cls)
        }
        if (it > 1) {
          rel <- abs(ll_trace[it] - ll_trace[it - 1]) / max(1, abs(ll_trace[it]))
          if (rel < private$.tol) { conv <- TRUE; break }
        }
      }
      e <- private$.estep(Xo, st)
      ord <- order(st$means[, 1])
      st$means <- st$means[ord, , drop = FALSE]
      st$vars <- st$vars[ord, , drop = FALSE]
      st$w <- st$w[ord]
      post <- matrix(NA_real_, nrow(X), k)
      post[ok, ] <- e$post[, ord, drop = FALSE]
      self$weights <- st$w
      self$means <- st$means
      self$variances <- st$vars
      colnames(self$means) <- colnames(self$variances) <- colnames(X)
      self$posterior <- post
      self$labels <- apply(post, 1, function(pr) if (anyNA(pr)) NA_integer_ else which.max(pr))
      self$loglik <- ll_trace
      self$iterations <- length(ll_trace)
      self$converged <- conv
      invisible(self)
    },

    #' @description Labels of new observations (maximum a posteriori).
    #' @param x Entropy features in the same format used in `$fit()`.
    #' @param is_map Treat a matrix `x` as a single map?
    predict = function(x, is_map = TRUE) {
      X <- private$.as_features(x, is_map, keep_dim = FALSE)
      ok <- stats::complete.cases(X)
      out <- rep(NA_integer_, nrow(X))
      st <- list(w = self$weights, means = self$means, vars = self$variances)
      e <- private$.estep(X[ok, , drop = FALSE], st)
      out[ok] <- apply(e$post, 1, which.max)
      if (is_map && (is.matrix(x) || is.array(x))) out <- matrix(out, dim(x)[1], dim(x)[2])
      out
    },

    #' @description Labels as a matrix with the spatial dimensions of the
    #'   fitted map.
    label_map = function() {
      if (is.null(self$dim)) stop("The segmenter was not fitted to a map.", call. = FALSE)
      matrix(self$labels, self$dim[1], self$dim[2])
    },

    #' @description Print a summary.
    #' @param ... Ignored.
    print = function(...) {
      cat("<EntropySegmenter> k =", self$k, "| method:", toupper(self$method), "\n")
      if (!is.null(self$labels)) {
        cat("  iterations:", self$iterations, "| converged:", self$converged, "\n")
        cat("  weights:", paste(format(round(self$weights, 3)), collapse = " "), "\n")
        cat("  class sizes:", paste(tabulate(self$labels, self$k), collapse = " "), "\n")
      }
      invisible(self)
    }
  ),
  private = list(
    .max_iter = 200L,
    .tol = 1e-8,
    .as_features = function(x, is_map, keep_dim = TRUE) {
      if (is.array(x) && length(dim(x)) == 3L) {
        d <- dim(x)
        if (keep_dim) self$dim <- d[1:2]
        X <- matrix(x, d[1] * d[2], d[3])
        colnames(X) <- dimnames(x)[[3]] %||% paste0("H", seq_len(d[3]))
      } else if (is.matrix(x) && is_map) {
        if (keep_dim) self$dim <- dim(x)
        X <- matrix(as.vector(x), ncol = 1, dimnames = list(NULL, "H"))
      } else if (is.matrix(x) || is.data.frame(x)) {
        if (keep_dim) self$dim <- NULL
        X <- as.matrix(x)
        if (is.null(colnames(X))) colnames(X) <- paste0("H", seq_len(ncol(X)))
      } else {
        if (keep_dim) self$dim <- NULL
        X <- matrix(as.numeric(x), ncol = 1, dimnames = list(NULL, "H"))
      }
      X[!is.finite(X)] <- NA
      X
    },
    .mstep = function(X, cls) {
      k <- self$k
      cls <- as.integer(cls)
      # keep every class populated
      for (g in seq_len(k)) if (!any(cls == g)) cls[sample.int(length(cls), 1)] <- g
      w <- tabulate(cls, k) / length(cls)
      means <- t(vapply(seq_len(k), function(g) colMeans(X[cls == g, , drop = FALSE]),
                        numeric(ncol(X))))
      vars <- t(vapply(seq_len(k), function(g) {
        v <- apply(X[cls == g, , drop = FALSE], 2, stats::var)
        v[!is.finite(v)] <- 0
        v
      }, numeric(ncol(X))))
      if (ncol(X) == 1L) { means <- t(means); vars <- t(vars) }
      vars <- pmax(vars, 1e-8 * max(1, mean(apply(X, 2, stats::var))))
      list(w = w, means = means, vars = vars)
    },
    .mstep_soft = function(X, post) {
      w <- colMeans(post)
      ng <- colSums(post)
      means <- t(post) %*% X / ng
      vars <- t(post) %*% (X^2) / ng - means^2
      vars <- pmax(vars, 1e-8 * max(1, mean(apply(X, 2, stats::var))))
      list(w = w, means = means, vars = vars)
    },
    .estep = function(X, st) {
      k <- self$k
      lp <- vapply(seq_len(k), function(g) {
        log(st$w[g]) + rowSums(stats::dnorm(X, rep(st$means[g, ], each = nrow(X)),
                                            rep(sqrt(st$vars[g, ]), each = nrow(X)),
                                            log = TRUE) |> matrix(nrow(X)))
      }, numeric(nrow(X)))
      if (!is.matrix(lp)) lp <- matrix(lp, nrow(X))
      mx <- apply(lp, 1, max)
      post <- exp(lp - mx)
      s <- rowSums(post)
      list(post = post / s, loglik = sum(mx + log(s)))
    }
  )
)

#' Entropy-based segmentation of a SAR/PolSAR image
#'
#' Convenience wrapper that computes the entropy map(s) with
#' [entropy_map()] and segments them with [EntropySegmenter].
#'
#' @inheritParams entropy_map
#' @param k Number of classes.
#' @param method `"em"` or `"sem"`.
#' @param ... Passed to `EntropySegmenter$new()`.
#'
#' @return A fitted [EntropySegmenter]; the labels are available with
#'   `$label_map()`.
#'
#' @examples
#' set.seed(1)
#' seg <- segment_sar(SanFrancisco150, k = 3, model = "gamma", window = 7, L = 4)
#' table(seg$label_map(), SanFrancisco150_gt)
#'
#' @export
segment_sar <- function(image, k = 2, model = c("wishart", "g0p", "gamma", "g0i", "ggamma"),
                        window = 7, L = NULL, channels = NULL,
                        method = c("em", "sem"), ...) {
  model <- match.arg(model)
  H <- entropy_map(image, model = model, window = window, L = L, channels = channels)
  EntropySegmenter$new(k = k, method = match.arg(method), ...)$fit(H)
}
