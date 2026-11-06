#' Divergence-based classifier for SAR/PolSAR features (R6)
#'
#' Implements the family of divergence-based classifiers of Ferreira, Coelho
#' and Nascimento (2021). Each pixel is represented by the model fitted to
#' the observations of a moving window (see [sar_features()]) and each class
#' by the model fitted to its training pixels. A pixel is assigned to the
#' class that minimises the chosen divergence between the two fitted models:
#' \deqn{\hat\imath = \arg\min_{i \in \{1, \ldots, g\}} d(\hat\theta_h,
#'   \hat\theta_i).}
#'
#' Two families of features are available:
#' \describe{
#'   \item{`"gaussian"`}{multivariate normal laws for the intensity vector,
#'     using the columns `mean_*` and `cov_*` of [sar_features()]. The class
#'     covariance combines the average local covariance and the covariance
#'     of the local means (law of total variance).}
#'   \item{`"wishart"`}{scaled complex Wishart laws for the full PolSAR
#'     covariance matrix, using the columns `c11`, `c12_re`, `c12_im`, ...
#'     of [sar_features()]. Since all laws share the number of looks, it
#'     only rescales the divergences and does not affect the classification.}
#' }
#'
#' The classifier is also available as the 'parsnip' model
#' [divergence_classifier()], for use in 'tidymodels' workflows.
#'
#' @field family `"gaussian"` or `"wishart"`.
#' @field divergence Divergence used (see [divergence_gaussian()]).
#' @field order Order of the Rényi divergence.
#' @field kind Form (1-4) of the Rényi divergence.
#' @field classes Class labels.
#' @field prototypes List with the fitted model of each class.
#' @field channels Channels (Gaussian) or matrix order (Wishart) used.
#'
#' @references
#' Ferreira, J. A., Coelho, H. and Nascimento, A. D. C. (2021). A family of
#' divergence-based classifiers for polarimetric synthetic aperture radar
#' (PolSAR) imagery vector and matrix features. *International Journal of
#' Remote Sensing*, 42(4), 1201--1229. \doi{10.1080/01431161.2020.1826060}
#'
#' @examples
#' img <- SARImage$new(SanFrancisco150)$crop(1:80, 1:80)
#' feats <- sar_features(img, window = 5)
#' feats$class <- factor(as.vector(SanFrancisco150_gt[1:80, 1:80]),
#'                       labels = c("ocean", "forest", "urban"))
#' set.seed(1)
#' train <- feats[sample(nrow(feats), 300), ]
#'
#' clf <- DivergenceClassifier$new(family = "wishart", divergence = "renyi",
#'                                 order = 0.5, kind = 3)
#' clf$fit(train, train$class)
#' pred <- clf$predict(feats)
#' mean(pred == feats$class)
#' head(clf$predict(feats, type = "distance"))
#'
#' @return An [R6::R6Class] generator. `DivergenceClassifier$new()` returns a new `DivergenceClassifier`
#'   object; its methods are described in the Methods section below.
#' @export
DivergenceClassifier <- R6::R6Class(
  "DivergenceClassifier",
  public = list(
    family = NULL,
    divergence = NULL,
    order = NULL,
    kind = NULL,
    classes = NULL,
    prototypes = NULL,
    channels = NULL,

    #' @description Create a classifier.
    #' @param family `"gaussian"` or `"wishart"`.
    #' @param divergence `"kl"`, `"jeffreys"`, `"bhattacharyya"` or `"renyi"`.
    #' @param order Order of the Rényi divergence, in (0, 1).
    #' @param kind Form of the Rényi divergence (1 to 4).
    initialize = function(family = c("gaussian", "wishart"),
                          divergence = c("renyi", "kl", "jeffreys", "bhattacharyya"),
                          order = 0.5, kind = 1) {
      self$family <- match.arg(family)
      self$divergence <- .check_div_args(match.arg(divergence), order, kind)
      self$order <- order
      self$kind <- as.integer(kind)
    },

    #' @description Fit the class prototypes.
    #' @param x Data frame of features (see [sar_features()]); other columns
    #'   are ignored.
    #' @param y Factor (or vector) of class labels.
    fit = function(x, y) {
      y <- as.factor(y)
      if (length(y) != nrow(x)) stop("`x` and `y` have different lengths.", call. = FALSE)
      if (anyNA(y)) stop("`y` has missing values.", call. = FALSE)
      self$classes <- levels(droplevels(y))
      y <- droplevels(y)
      info <- private$.parse(x)
      self$channels <- info$channels
      self$prototypes <- lapply(self$classes, function(cl) {
        idx <- which(y == cl)
        if (self$family == "gaussian") {
          M <- info$means[idx, , drop = FALSE]
          mu <- colMeans(M)
          within <- Reduce(`+`, lapply(idx, function(r) info$cov(r))) / length(idx)
          between <- if (length(idx) > 1) stats::cov(M) else 0 * within
          S <- .regularize(within + between)
          list(mu = mu, Sigma = S, logdet = .hlogdet(S))
        } else {
          S <- .regularize(.vec_to_herm(colMeans(info$vec[idx, , drop = FALSE])))
          list(Sigma = S, logdet = .hlogdet(S), inv = .hinv(S))
        }
      })
      names(self$prototypes) <- self$classes
      invisible(self)
    },

    #' @description Divergence between every row of `x` and every class.
    #' @param x Data frame of features.
    distances = function(x) {
      if (is.null(self$prototypes)) stop("Fit the classifier first.", call. = FALSE)
      info <- private$.parse(x)
      n <- nrow(x)
      D <- matrix(NA_real_, n, length(self$classes),
                  dimnames = list(NULL, self$classes))
      for (r in seq_len(n)) {
        if (self$family == "gaussian") {
          mu <- info$means[r, ]
          S <- .regularize(info$cov(r))
          if (anyNA(mu) || anyNA(S)) next
          for (k in seq_along(self$classes)) {
            pk <- self$prototypes[[k]]
            D[r, k] <- tryCatch(
              .div_gauss(mu, S, pk$mu, pk$Sigma, self$divergence, self$order, self$kind),
              error = function(e) NA_real_)
          }
        } else {
          v <- info$vec[r, ]
          if (anyNA(v)) next
          S <- .regularize(.vec_to_herm(v))
          ld <- .hlogdet(S)
          Si <- tryCatch(.hinv(S), error = function(e) NULL)
          if (is.null(Si) || !is.finite(ld)) next
          for (k in seq_along(self$classes)) {
            pk <- self$prototypes[[k]]
            D[r, k] <- .div_wish(S, pk$Sigma, 1, self$divergence, self$order, self$kind,
                                 ld1 = ld, ld2 = pk$logdet, S1i = Si, S2i = pk$inv)
          }
        }
      }
      D
    },

    #' @description Predict classes, pseudo-probabilities or distances.
    #' @param x Data frame of features.
    #' @param type `"class"`, `"prob"` (softmax of minus the divergences) or
    #'   `"distance"`.
    predict = function(x, type = c("class", "prob", "distance")) {
      type <- match.arg(type)
      D <- self$distances(x)
      if (type == "distance") return(tibble::as_tibble(as.data.frame(D)))
      if (type == "prob") {
        P <- t(apply(D, 1, function(d) {
          if (all(is.na(d))) return(rep(NA_real_, length(d)))
          d[is.na(d)] <- Inf
          w <- exp(-(d - min(d)))
          w / sum(w)
        }))
        if (length(self$classes) == 1L) P <- t(P)
        colnames(P) <- self$classes
        return(tibble::as_tibble(as.data.frame(P)))
      }
      cls <- apply(D, 1, function(d) if (all(is.na(d))) NA_integer_ else which.min(d))
      factor(self$classes[cls], levels = self$classes)
    },

    #' @description Print the classifier.
    #' @param ... Ignored.
    print = function(...) {
      cat("<DivergenceClassifier> family:", self$family, "| divergence:", self$divergence)
      if (self$divergence == "renyi") cat(" (order ", self$order, ", kind ", self$kind, ")", sep = "")
      cat("\n")
      if (!is.null(self$classes)) {
        cat("  classes:", paste(self$classes, collapse = ", "), "\n")
      } else cat("  (not fitted)\n")
      invisible(self)
    }
  ),
  private = list(
    .parse = function(x) {
      x <- as.data.frame(x)
      nms <- names(x)
      if (self$family == "gaussian") {
        mcols <- grep("^mean_", nms, value = TRUE)
        if (length(mcols) == 0) {
          stop("Gaussian features need `mean_*` and `cov_*` columns (see `sar_features()`).",
               call. = FALSE)
        }
        ch <- sub("^mean_", "", mcols)
        if (!is.null(self$channels) && !identical(ch, self$channels)) {
          stop("The feature columns differ from those used in the fit.", call. = FALSE)
        }
        p <- length(ch)
        cov_cols <- matrix(NA_character_, p, p)
        for (a in seq_len(p)) for (b in a:p) {
          nm <- paste("cov", ch[a], ch[b], sep = "_")
          if (!nm %in% nms) stop("Missing column `", nm, "`.", call. = FALSE)
          cov_cols[a, b] <- cov_cols[b, a] <- nm
        }
        means <- as.matrix(x[, mcols, drop = FALSE])
        C <- as.matrix(x[, as.vector(cov_cols), drop = FALSE])
        list(channels = ch, means = means,
             cov = function(r) matrix(C[r, ], p, p))
      } else {
        dcols <- grep("^c([0-9])\\1$", nms, value = TRUE)
        if (length(dcols) == 0) {
          stop("Wishart features need `c11`, `c12_re`, `c12_im`, ... columns (see `sar_features()`).",
               call. = FALSE)
        }
        m <- length(dcols)
        if (!is.null(self$channels) && !identical(m, self$channels)) {
          stop("The feature columns differ from those used in the fit.", call. = FALSE)
        }
        pr <- .upper_pairs(m)
        V <- matrix(0 + 0i, nrow(x), m + nrow(pr))
        for (k in seq_len(m)) V[, k] <- x[[paste0("c", k, k)]]
        for (k in seq_len(nrow(pr))) {
          nm <- paste0("c", pr[k, 1], pr[k, 2])
          if (!all(paste0(nm, c("_re", "_im")) %in% nms)) {
            stop("Missing columns `", nm, "_re`/`", nm, "_im`.", call. = FALSE)
          }
          V[, m + k] <- complex(real = x[[paste0(nm, "_re")]],
                                imaginary = x[[paste0(nm, "_im")]])
        }
        list(channels = m, vec = V)
      }
    }
  )
)

#' @export
predict.DivergenceClassifier <- function(object, new_data, type = "class", ...) {
  object$predict(new_data, type = type)
}
