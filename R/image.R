#' SAR and PolSAR images (R6)
#'
#' `SARImage` stores a SAR image as an array of dimension
#' `nrow x ncol x layers`. Two types are supported:
#' \describe{
#'   \item{`"polsar"`}{a complex array whose `p = m(m+1)/2` layers are the
#'     distinct entries of the \eqn{m \times m} covariance matrix of each
#'     pixel, in the order \eqn{C_{11}, \ldots, C_{mm}, C_{12}, C_{13},
#'     \ldots} (for \eqn{m = 3}: `HH, HV, VV, HH.HV*, HH.VV*, HV.VV*`), as
#'     in [SanFrancisco150];}
#'   \item{`"intensity"`}{a real (or real-valued complex) array whose layers
#'     are intensity channels.}
#' }
#'
#' The class offers cropping, Pauli colour composition, plotting, local
#' (moving-window) features for classification, entropy maps and conversion
#' to tibbles.
#'
#' @field data The image array.
#' @field type `"polsar"` or `"intensity"`.
#' @field channels Names of the intensity channels.
#' @field looks Nominal number of looks (optional, `NULL` if unknown).
#' @field nrow Number of rows (read-only).
#' @field ncol Number of columns (read-only).
#' @field m Number of channels; matrix order for PolSAR images (read-only).
#'
#' @examples
#' img <- SARImage$new(SanFrancisco150, looks = 4)
#' img
#' small <- img$crop(1:40, 1:40)
#' dim(small$intensities())
#' small$pixel(1, 1)
#' head(small$features(window = 3))
#' \donttest{
#' img$plot()
#' }
#'
#' @return An [R6::R6Class] generator. `SARImage$new()` returns a new `SARImage`
#'   object; its methods are described in the Methods section below.
#' @export
SARImage <- R6::R6Class(
  "SARImage",
  public = list(
    data = NULL,
    type = NULL,
    channels = NULL,
    looks = NULL,

    #' @description Create a SAR image.
    #' @param data Array `nrow x ncol x layers` (a matrix is treated as a
    #'   single intensity channel).
    #' @param type `"polsar"` or `"intensity"`. Guessed when `NULL`: complex
    #'   arrays with `m(m+1)/2` layers (`m >= 2`) and non-zero imaginary
    #'   parts are PolSAR images.
    #' @param looks Nominal number of looks.
    #' @param channels Names of the intensity channels.
    initialize = function(data, type = NULL, looks = NULL, channels = NULL) {
      if (is.matrix(data)) data <- array(data, dim = c(dim(data), 1L))
      if (!is.array(data) || length(dim(data)) != 3L) {
        stop("`data` must be a matrix or a 3-dimensional array.", call. = FALSE)
      }
      p <- dim(data)[3]
      if (is.null(type)) {
        m <- (sqrt(8 * p + 1) - 1) / 2
        type <- if (is.complex(data) && p > 1 && abs(m - round(m)) < 1e-8 &&
                    any(Im(data) != 0)) "polsar" else "intensity"
      }
      type <- match.arg(type, c("polsar", "intensity"))
      if (type == "polsar") {
        m <- .m_from_p(p)
        if (is.null(channels)) {
          channels <- if (m == 3L) c("HH", "HV", "VV") else paste0("C", seq_len(m))
        }
        if (!is.complex(data)) data <- data + 0i
      } else {
        if (is.complex(data)) data <- Re(data)
        m <- p
        if (is.null(channels)) {
          channels <- if (p == 3L) c("HH", "HV", "VV") else paste0("I", seq_len(p))
        }
      }
      if (length(channels) != m) stop("`channels` has the wrong length.", call. = FALSE)
      if (!is.null(looks)) .check_scalar(looks, "looks", lower = 0)
      self$data <- data
      self$type <- type
      self$channels <- channels
      self$looks <- looks
    },

    #' @description Intensity of one channel.
    #' @param channel Channel index or name.
    intensity = function(channel = 1) {
      k <- private$.channel_index(channel)
      Re(self$data[, , k])
    },

    #' @description All intensity channels as a real array.
    intensities = function() {
      Re(self$data[, , seq_len(self$m), drop = FALSE])
    },

    #' @description Crop the image.
    #' @param rows,cols Row and column indices to keep.
    crop = function(rows, cols) {
      SARImage$new(self$data[rows, cols, , drop = FALSE], type = self$type,
                   looks = self$looks, channels = self$channels)
    },

    #' @description Observation at a pixel: the Hermitian covariance matrix
    #'   (PolSAR) or the vector of intensities.
    #' @param row,col Pixel position.
    pixel = function(row, col) {
      v <- self$data[row, col, ]
      if (self$type == "polsar") .vec_to_herm(v) else v
    },

    #' @description Observations inside a square window centred at a pixel
    #'   (truncated at the borders). Returns an `n x layers` matrix.
    #' @param row,col Pixel position.
    #' @param window Odd window side.
    window = function(row, col, window = 3) {
      .window_values(self$data, row, col, .check_window(window))
    },

    #' @description Pauli RGB composition scaled to \[0, 1\] (PolSAR images
    #'   with three channels); for other images the first three intensities
    #'   are used.
    #' @param stretch Upper quantile used for linear stretching.
    #' @param gamma Gamma correction exponent.
    pauli = function(stretch = 0.98, gamma = 0.5) {
      if (self$type == "polsar" && self$m == 3L) {
        C11 <- Re(self$data[, , 1]); C22 <- Re(self$data[, , 2])
        C33 <- Re(self$data[, , 3]); C13 <- Re(self$data[, , 5])
        rgb <- list((C11 + C33 - 2 * C13) / 2, 2 * C22, (C11 + C33 + 2 * C13) / 2)
      } else {
        I <- self$intensities()
        idx <- rep_len(seq_len(dim(I)[3]), 3)
        rgb <- lapply(idx, function(k) I[, , k])
      }
      sc <- function(v) {
        v[!is.finite(v) | v < 0] <- 0
        q <- stats::quantile(v, stretch, names = FALSE)
        if (q <= 0) q <- max(v, 1)
        pmin(v / q, 1)^gamma
      }
      array(unlist(lapply(rgb, sc)), dim = c(self$nrow, self$ncol, 3L))
    },

    #' @description Plot the image (Pauli composition or one channel in dB).
    #' @param type `"pauli"` or `"intensity"`.
    #' @param channel Channel for `type = "intensity"`.
    #' @param main Title.
    #' @param ... Passed to [graphics::rasterImage()].
    plot = function(type = c("pauli", "intensity"), channel = 1, main = NULL, ...) {
      type <- match.arg(type)
      if (type == "pauli") {
        ras <- grDevices::as.raster(self$pauli())
        main <- main %||% "Pauli RGB"
      } else {
        db <- 10 * log10(pmax(self$intensity(channel), .Machine$double.eps))
        q <- stats::quantile(db, c(0.02, 0.98), names = FALSE)
        v <- pmin(pmax((db - q[1]) / diff(q), 0), 1)
        ras <- grDevices::as.raster(matrix(grDevices::grey(v), nrow(v)))
        main <- main %||% paste(self$channels[private$.channel_index(channel)], "(dB)")
      }
      .plot_raster(ras, main, ...)
      invisible(self)
    },

    #' @description Local features for classification; see [sar_features()].
    #' @param window Odd window side.
    #' @param type Which features: `"all"`, `"gaussian"` or `"wishart"`.
    features = function(window = 3, type = c("all", "gaussian", "wishart")) {
      sar_features(self, window = window, type = type)
    },

    #' @description Entropy map; see [entropy_map()].
    #' @param ... Passed to [entropy_map()].
    entropy_map = function(...) entropy_map(self, ...),

    #' @description Convert to a tibble with one row per pixel.
    as_tibble = function() {
      out <- tibble::tibble(row = rep(seq_len(self$nrow), times = self$ncol),
                            col = rep(seq_len(self$ncol), each = self$nrow))
      I <- self$intensities()
      for (k in seq_len(self$m)) out[[self$channels[k]]] <- as.vector(I[, , k])
      if (self$type == "polsar" && self$m > 1) {
        pr <- .upper_pairs(self$m)
        for (k in seq_len(nrow(pr))) {
          nm <- paste(self$channels[pr[k, 1]], self$channels[pr[k, 2]], sep = "_")
          z <- as.vector(self$data[, , self$m + k])
          out[[paste0(nm, "_re")]] <- Re(z)
          out[[paste0(nm, "_im")]] <- Im(z)
        }
      }
      out
    },

    #' @description Print a summary.
    #' @param ... Ignored.
    print = function(...) {
      cat("<SARImage> ", self$nrow, " x ", self$ncol, " ",
          if (self$type == "polsar") paste0("PolSAR (", self$m, " x ", self$m, " covariance)")
          else paste0("intensity image (", self$m, " channels)"), "\n", sep = "")
      cat("  channels: ", paste(self$channels, collapse = ", "), "\n", sep = "")
      if (!is.null(self$looks)) cat("  looks: ", self$looks, "\n", sep = "")
      I <- self$intensities()
      for (k in seq_len(self$m)) {
        cat(sprintf("  mean %-3s = %.4g\n", self$channels[k], mean(I[, , k])))
      }
      invisible(self)
    }
  ),
  active = list(
    nrow = function() dim(self$data)[1],
    ncol = function() dim(self$data)[2],
    m = function() {
      if (self$type == "polsar") .m_from_p(dim(self$data)[3]) else dim(self$data)[3]
    }
  ),
  private = list(
    .channel_index = function(channel) {
      if (is.character(channel)) {
        k <- match(channel, self$channels)
        if (is.na(k)) stop("Unknown channel: ", channel, call. = FALSE)
        return(k)
      }
      if (channel < 1 || channel > self$m) stop("Invalid channel index.", call. = FALSE)
      as.integer(channel)
    }
  )
)

.as_sarimage <- function(image, looks = NULL) {
  if (inherits(image, "SARImage")) return(image)
  SARImage$new(image, looks = looks)
}

.plot_raster <- function(ras, main = NULL, legend = NULL, ...) {
  op <- graphics::par(mar = c(0.5, 0.5, if (is.null(main)) 0.5 else 2, 0.5))
  on.exit(graphics::par(op))
  nr <- nrow(ras); nc <- ncol(ras)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, nc), ylim = c(0, nr), asp = 1, xaxs = "i", yaxs = "i")
  graphics::rasterImage(ras, 0, 0, nc, nr, interpolate = FALSE, ...)
  if (!is.null(main)) graphics::title(main)
  if (!is.null(legend)) {
    graphics::legend("bottomright", legend = legend$labels, fill = legend$fill,
                     bg = "#FFFFFFD9", box.col = "#666666", cex = 0.8, inset = 0.02)
  }
}

#' Plot a map
#'
#' Displays a matrix as an image with the first row on top: continuous maps
#' (e.g. entropy maps) use a sequential palette and integer/factor maps
#' (e.g. segmentation or classification results) use a qualitative one.
#'
#' @param x A numeric, integer, logical or factor matrix.
#' @param palette Optional vector of colours.
#' @param main Title.
#' @param legend Add a legend for discrete maps?
#'
#' @return Invisibly, `x`.
#'
#' @examples
#' plot_map(SanFrancisco150_gt, main = "Ground truth", legend = TRUE)
#' plot_map(matrix(rnorm(100), 10), main = "Noise")
#'
#' @export
plot_map <- function(x, palette = NULL, main = NULL, legend = FALSE) {
  if (!is.matrix(x)) stop("`x` must be a matrix.", call. = FALSE)
  discrete <- is.factor(x) || is.logical(x) || is.integer(x) || is.character(x)
  if (discrete) {
    f <- factor(as.vector(x))
    lev <- levels(f)
    palette <- palette %||% c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                              "#66A61E", "#E6AB02", "#A6761D", "#666666")
    palette <- rep_len(palette, length(lev))
    cols <- palette[as.integer(f)]
    cols[is.na(cols)] <- "#FFFFFF"
  } else {
    v <- as.vector(x)
    rng <- range(v[is.finite(v)])
    palette <- palette %||% grDevices::hcl.colors(256, "viridis")
    idx <- 1 + floor((v - rng[1]) / max(diff(rng), .Machine$double.eps) * (length(palette) - 1))
    cols <- palette[idx]
    cols[is.na(cols)] <- "#FFFFFF"
  }
  ras <- grDevices::as.raster(matrix(cols, nrow(x)))
  leg <- if (discrete && legend) list(labels = lev, fill = palette) else NULL
  .plot_raster(ras, main, legend = leg)
  invisible(x)
}

#' Local features of SAR images for pixel classification
#'
#' Computes moving-window statistics around every pixel, in the format used
#' by [DivergenceClassifier] and the [divergence_classifier()] 'parsnip'
#' model. Window statistics are computed with integral images, so the cost
#' does not depend on the window size. Windows are truncated at the image
#' borders.
#'
#' Two groups of columns are produced:
#' \describe{
#'   \item{Gaussian features}{`mean_<ch>`: local mean of each intensity
#'     channel; `cov_<ch1>_<ch2>`: local (unbiased) covariance between the
#'     intensities. Together they define the multivariate normal fitted to
#'     the intensity vector in the window.}
#'   \item{Wishart features (PolSAR only)}{`c<i><i>` and
#'     `c<i><j>_re`/`c<i><j>_im`: entries of the local mean covariance
#'     matrix, i.e. the maximum likelihood estimate of \eqn{\Sigma} for the
#'     scaled complex Wishart law in the window.}
#' }
#'
#' @param image A [SARImage] or an array accepted by `SARImage$new()`.
#' @param window Odd window side.
#' @param type Which features: `"all"`, `"gaussian"` or `"wishart"`.
#'
#' @return A tibble with the columns `row`, `col` and the features.
#'
#' @examples
#' img <- SARImage$new(SanFrancisco150)$crop(1:20, 1:20)
#' f <- sar_features(img, window = 3)
#' names(f)
#'
#' @export
sar_features <- function(image, window = 3, type = c("all", "gaussian", "wishart")) {
  type <- match.arg(type)
  img <- .as_sarimage(image)
  window <- .check_window(window)
  nr <- img$nrow; nc <- img$ncol; m <- img$m
  out <- tibble::tibble(row = rep(seq_len(nr), times = nc),
                        col = rep(seq_len(nc), each = nr))
  I <- img$intensities()
  ch <- img$channels
  if (type %in% c("all", "gaussian")) {
    cnt <- .box_count(nr, nc, window)
    corr <- ifelse(cnt > 1, cnt / (cnt - 1), 1)
    mu <- lapply(seq_len(m), function(k) .box_mean(I[, , k], window))
    for (k in seq_len(m)) out[[paste0("mean_", ch[k])]] <- as.vector(mu[[k]])
    for (a in seq_len(m)) for (b in a:m) {
      cv <- (.box_mean(I[, , a] * I[, , b], window) - mu[[a]] * mu[[b]]) * corr
      out[[paste("cov", ch[a], ch[b], sep = "_")]] <- as.vector(cv)
    }
  }
  if (type %in% c("all", "wishart")) {
    if (img$type != "polsar") {
      if (type == "wishart") stop("Wishart features need a PolSAR image.", call. = FALSE)
    } else {
      for (k in seq_len(m)) {
        out[[paste0("c", k, k)]] <- as.vector(.box_mean(Re(img$data[, , k]), window))
      }
      pr <- .upper_pairs(m)
      for (k in seq_len(nrow(pr))) {
        z <- .box_mean(img$data[, , m + k], window)
        nm <- paste0("c", pr[k, 1], pr[k, 2])
        out[[paste0(nm, "_re")]] <- as.vector(Re(z))
        out[[paste0(nm, "_im")]] <- as.vector(Im(z))
      }
    }
  }
  out
}

#' @importFrom tibble as_tibble
#' @export
as_tibble.SARImage <- function(x, ...) x$as_tibble()

#' @export
tibble::as_tibble
