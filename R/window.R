# Moving-window tools -------------------------------------------------------

.check_window <- function(window) {
  if (!is.numeric(window) || length(window) != 1L || window < 1 ||
      window %% 2 != 1) {
    stop("`window` must be an odd positive integer (e.g. 3, 5, 7).", call. = FALSE)
  }
  as.integer(window)
}

# Index ranges of the (border-truncated) square window around each row/col
.window_ranges <- function(n, half) {
  i <- seq_len(n)
  list(lo = pmax(1L, i - half), hi = pmin(n, i + half))
}

# Local mean of a (real or complex) matrix over a square window, with the
# window truncated at the image borders. Uses an integral image, so the cost
# does not depend on the window size.
.box_mean <- function(M, window) {
  half <- (window - 1L) %/% 2L
  nr <- nrow(M); nc <- ncol(M)
  S <- matrix(if (is.complex(M)) 0 + 0i else 0, nr + 1L, nc + 1L)
  cs <- matrix(apply(M, 2, cumsum), nr)
  S[-1, -1] <- matrix(t(apply(cs, 1, cumsum)), nr)
  r <- .window_ranges(nr, half)
  c <- .window_ranges(nc, half)
  r1 <- rep(r$lo, times = nc); r2 <- rep(r$hi, times = nc)
  c1 <- rep(c$lo, each = nr);  c2 <- rep(c$hi, each = nr)
  tot <- S[cbind(r2 + 1L, c2 + 1L)] - S[cbind(r1, c2 + 1L)] -
    S[cbind(r2 + 1L, c1)] + S[cbind(r1, c1)]
  cnt <- (r2 - r1 + 1L) * (c2 - c1 + 1L)
  matrix(tot / cnt, nr, nc)
}

# Number of pixels in each (truncated) window
.box_count <- function(nr, nc, window) {
  half <- (window - 1L) %/% 2L
  r <- .window_ranges(nr, half)
  c <- .window_ranges(nc, half)
  outer(r$hi - r$lo + 1L, c$hi - c$lo + 1L)
}

# Values of the layers of `arr` (nr x nc x p) inside the window at (i, j)
.window_values <- function(arr, i, j, window) {
  half <- (window - 1L) %/% 2L
  d <- dim(arr)
  rows <- max(1L, i - half):min(d[1], i + half)
  cols <- max(1L, j - half):min(d[2], j + half)
  if (length(d) == 2L) return(as.vector(arr[rows, cols]))
  matrix(arr[rows, cols, , drop = FALSE], ncol = d[3])
}
