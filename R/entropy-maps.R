.intensity_models <- c("gamma", "g0i", "ggamma")
.matrix_models <- c("wishart", "g0p")

.new_model <- function(model) {
  switch(model,
         gamma = GammaSAR$new(),
         g0i = G0I$new(),
         ggamma = GGammaSAR$new(),
         wishart = WishartSAR$new(),
         g0p = G0Polarimetric$new(),
         stop("Unknown model: ", model, call. = FALSE))
}

# fit a model; L is passed with the name expected by each class
.fit_model <- function(model, x, L = NULL) {
  obj <- .new_model(model)
  if (model == "ggamma") obj$fit(x, phi = L) else obj$fit(x, L = L)
}

#' Shannon entropy maps
#'
#' Estimates, at every pixel, the parameters of a SAR/PolSAR model with the
#' observations inside a square moving window and returns the map of
#' Shannon entropies (Ferreira and Nascimento, 2020; Nascimento, Ferreira
#' and Frery, 2023). Optionally, the map of asymptotic variances of the
#' entropy estimators is returned as well.
#'
#' For the Gamma and Wishart models the computations are fully vectorised
#' (local sufficient statistics are obtained from integral images), so the
#' maps are computed in a fraction of a second. The G0I, generalized Gamma
#' and G0 polarimetric models require a numerical maximisation per pixel and
#' are therefore slower (use small images or crops to explore them).
#'
#' @param image A [SARImage] or an array accepted by `SARImage$new()`.
#' @param model One of `"gamma"`, `"g0i"`, `"ggamma"` (intensity models,
#'   applied to each channel) or `"wishart"`, `"g0p"` (matrix models, for
#'   PolSAR images).
#' @param window Odd side of the moving window.
#' @param L Known number of looks. If `NULL` (default) it is estimated in
#'   every window.
#' @param channels Channels used by intensity models (default: all).
#' @param variance Also return the variance map?
#'
#' @return A matrix (one channel or matrix model) or an array
#'   `nrow x ncol x channels` of entropies. If `variance = TRUE`, a list with
#'   elements `entropy`, `variance` and `n` (window sizes).
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
#' img <- SARImage$new(SanFrancisco150)
#' H <- entropy_map(img, model = "wishart", window = 7, L = 4)
#' dim(H)
#' plot_map(H, main = "Wishart entropy")
#'
#' Hg <- entropy_map(img, model = "gamma", window = 7, channels = "HV")
#' \donttest{
#' small <- img$crop(1:15, 1:15)
#' Hi <- entropy_map(small, model = "g0i", window = 5, L = 4, channels = 1)
#' }
#'
#' @export
entropy_map <- function(image, model = c("gamma", "g0i", "ggamma", "wishart", "g0p"),
                        window = 7, L = NULL, channels = NULL, variance = FALSE) {
  model <- match.arg(model)
  img <- .as_sarimage(image)
  window <- .check_window(window)
  if (!is.null(L)) .check_scalar(L, "L", lower = 0)
  nr <- img$nrow; nc <- img$ncol
  cnt <- .box_count(nr, nc, window)
  if (model %in% .matrix_models) {
    if (img$type != "polsar") stop("Matrix models need a PolSAR image.", call. = FALSE)
    res <- if (model == "wishart") {
      .map_wishart(img, window, L, cnt)
    } else {
      .map_loop(nr, nc, function(i, j) img$window(i, j, window), model, L, variance)
    }
    return(if (variance) list(entropy = res$entropy, variance = res$variance, n = cnt)
           else res$entropy)
  }
  channels <- channels %||% seq_len(img$m)
  idx <- vapply(channels, function(ch) {
    if (is.character(ch)) match(ch, img$channels) else as.integer(ch)
  }, integer(1))
  if (anyNA(idx) || any(idx < 1 | idx > img$m)) stop("Invalid `channels`.", call. = FALSE)
  I <- img$intensities()
  ent <- array(NA_real_, dim = c(nr, nc, length(idx)),
               dimnames = list(NULL, NULL, img$channels[idx]))
  vr <- ent
  for (k in seq_along(idx)) {
    X <- I[, , idx[k]]
    res <- if (model == "gamma") {
      .map_gamma(X, window, L, cnt)
    } else {
      .map_loop(nr, nc, function(i, j) .window_values(X, i, j, window), model, L, variance)
    }
    ent[, , k] <- res$entropy
    if (variance) vr[, , k] <- res$variance
  }
  if (length(idx) == 1L) { ent <- ent[, , 1]; vr <- vr[, , 1] }
  if (variance) list(entropy = ent, variance = vr, n = cnt) else ent
}

.map_gamma <- function(X, window, L, cnt) {
  X[!(X > 0)] <- NA
  mu <- .box_mean(X, window)
  if (is.null(L)) {
    lm <- .box_mean(log(X), window)
    s <- pmax(log(mu) - lm, 1e-12)
    Lh <- .vec_bisect_log(function(l) log(l) - digamma(l), as.vector(s), 1e-4, 1e7)
    Lh <- matrix(Lh, nrow(X))
  } else {
    Lh <- matrix(L, nrow(X), ncol(X))
  }
  H <- entropy_gamma(pmax(mu, .Machine$double.xmin), Lh)
  H[!is.finite(mu)] <- NA
  V <- matrix(.entropy_var_wishart(as.vector(Lh), 1, as.vector(cnt), is.null(L)), nrow(X))
  list(entropy = matrix(H, nrow(X)), variance = V)
}

.map_wishart <- function(img, window, L, cnt) {
  m <- img$m
  nr <- img$nrow; nc <- img$ncol
  p <- dim(img$data)[3]
  Cbar <- matrix(0 + 0i, nr * nc, p)
  for (k in seq_len(p)) Cbar[, k] <- as.vector(.box_mean(img$data[, , k], window))
  ld <- suppressWarnings(log(.hdet_vec(Cbar, m)))
  if (is.null(L)) {
    pix <- matrix(img$data, nr * nc, p)
    lz <- suppressWarnings(log(.hdet_vec(pix, m)))
    lzm <- as.vector(.box_mean(matrix(lz, nr), window))
    s <- pmax(ld - lzm, 1e-12)
    Lh <- .vec_bisect_log(function(l) m * log(l) - .mvpsigamma(l, m), s,
                          m - 1 + 1e-6, 1e7)
  } else {
    Lh <- rep(L, nr * nc)
  }
  H <- rep(NA_real_, nr * nc)
  ok <- is.finite(ld) & is.finite(Lh)
  H[ok] <- .entropy_wishart_ld(ld[ok], Lh[ok], m)
  V <- rep(NA_real_, nr * nc)
  V[ok] <- .entropy_var_wishart(Lh[ok], m, as.vector(cnt)[ok], is.null(L))
  list(entropy = matrix(H, nr), variance = matrix(V, nr))
}

.map_loop <- function(nr, nc, get_values, model, L, variance) {
  H <- matrix(NA_real_, nr, nc)
  V <- matrix(NA_real_, nr, nc)
  for (j in seq_len(nc)) for (i in seq_len(nr)) {
    x <- get_values(i, j)
    r <- tryCatch({
      obj <- .fit_model(model, x, L)
      c(obj$entropy(), if (variance) obj$entropy_variance() else NA_real_)
    }, error = function(e) c(NA_real_, NA_real_))
    H[i, j] <- r[1]
    V[i, j] <- r[2]
  }
  list(entropy = H, variance = V)
}


#' Entropy-based test for two SAR samples
#'
#' Tests whether two samples come from the same distribution by comparing
#' their estimated Shannon entropies. Under the null hypothesis of equal
#' entropies, and by the asymptotic normality of the maximum likelihood
#' estimator of the entropy,
#' \deqn{S = \frac{(\hat H_1 - \hat H_2)^2}{\widehat{\mathrm{Var}}(\hat H_1) +
#'   \widehat{\mathrm{Var}}(\hat H_2)}}
#' follows asymptotically a \eqn{\chi^2_1} distribution. This statistic is
#' the basis of the entropy-based change detectors of the thesis of
#' J. A. Ferreira (2021).
#'
#' @param x,y Samples: numeric vectors of intensities (intensity models) or
#'   PolSAR observations (matrix models).
#' @param model One of `"gamma"`, `"g0i"`, `"ggamma"`, `"wishart"`, `"g0p"`.
#' @param L Known number of looks; estimated when `NULL`.
#'
#' @return An object of class `htest`.
#'
#' @references
#' Ferreira, J. A. (2021). *Entropias e somas compostas para matrizes
#' estocásticas aplicadas a sistemas PolSAR*. PhD thesis, Universidade
#' Federal de Pernambuco.
#'
#' @examples
#' set.seed(1)
#' x <- rg0i(200, alpha = -3, gamma = 2, L = 4)
#' y <- rg0i(200, alpha = -3, gamma = 2, L = 4)
#' z <- rg0i(200, alpha = -1.5, gamma = 1, L = 4)
#' entropy_test(x, y, model = "g0i", L = 4)
#' entropy_test(x, z, model = "g0i", L = 4)$p.value
#'
#' @export
entropy_test <- function(x, y, model = c("gamma", "g0i", "ggamma", "wishart", "g0p"),
                         L = NULL) {
  model <- match.arg(model)
  dname <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  mx <- .fit_model(model, x, L)
  my <- .fit_model(model, y, L)
  hx <- mx$entropy(); hy <- my$entropy()
  vx <- mx$entropy_variance(); vy <- my$entropy_variance()
  stat <- (hx - hy)^2 / (vx + vy)
  structure(list(
    statistic = c(S = stat),
    parameter = c(df = 1),
    p.value = stats::pchisq(stat, df = 1, lower.tail = FALSE),
    estimate = c(`entropy of x` = hx, `entropy of y` = hy),
    null.value = c(`difference in entropies` = 0),
    alternative = "two.sided",
    method = paste0("Shannon entropy test (", mx$name, " model)"),
    data.name = dname
  ), class = "htest")
}

#' Entropy-based change detection between two images
#'
#' Applies the entropy test of [entropy_test()] pixel by pixel, using the
#' observations of a moving window in two co-registered images, and returns
#' the map of p-values and the binary change map.
#'
#' @param image1,image2 Co-registered images ([SARImage] objects or arrays)
#'   of the same size.
#' @inheritParams entropy_map
#' @param channel Channel used by intensity models.
#' @param level Significance level used to declare a change.
#'
#' @return A list of class `sitsar_change` with the matrices `statistic`,
#'   `p_value` and `change` (logical), and the `level`.
#'
#' @examples
#' set.seed(1)
#' a <- matrix(rgammasar(40 * 40, mu = 1, L = 4), 40)
#' b <- a
#' b[11:30, 11:30] <- rgammasar(400, mu = 4, L = 4)
#' cd <- change_map(a, b, model = "gamma", window = 5, L = 4)
#' cd
#' plot_map(cd$change, main = "Detected changes")
#'
#' @export
change_map <- function(image1, image2, model = c("gamma", "g0i", "ggamma", "wishart", "g0p"),
                       window = 5, L = NULL, channel = 1, level = 0.01) {
  model <- match.arg(model)
  img1 <- .as_sarimage(image1); img2 <- .as_sarimage(image2)
  if (img1$nrow != img2$nrow || img1$ncol != img2$ncol) {
    stop("Images must have the same dimensions.", call. = FALSE)
  }
  .check_scalar(level, "level", lower = 0, upper = 1)
  ch <- if (model %in% .matrix_models) NULL else channel
  e1 <- entropy_map(img1, model, window, L, channels = ch, variance = TRUE)
  e2 <- entropy_map(img2, model, window, L, channels = ch, variance = TRUE)
  S <- (e1$entropy - e2$entropy)^2 / (e1$variance + e2$variance)
  p <- stats::pchisq(S, df = 1, lower.tail = FALSE)
  structure(list(statistic = S, p_value = p, change = p < level, level = level,
                 model = model),
            class = "sitsar_change")
}

#' @export
print.sitsar_change <- function(x, ...) {
  cat("<Entropy-based change map> model:", x$model, "| level:", x$level, "\n")
  cat("  size:", nrow(x$change), "x", ncol(x$change), "\n")
  cat("  changed pixels:", sum(x$change, na.rm = TRUE),
      sprintf("(%.2f%%)\n", 100 * mean(x$change, na.rm = TRUE)))
  invisible(x)
}
