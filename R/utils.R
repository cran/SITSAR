# Internal helpers ---------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.check_scalar <- function(x, name, lower = -Inf, upper = Inf,
                          lower_open = TRUE, upper_open = TRUE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) {
    ok <- if (lower_open) x > lower else x >= lower
    ok <- ok && (if (upper_open) x < upper else x <= upper)
  }
  if (!ok) {
    rng <- paste0(if (lower_open) "(" else "[", lower, ", ", upper,
                  if (upper_open) ")" else "]")
    stop(sprintf("`%s` must be a single finite number in %s.", name, rng),
         call. = FALSE)
  }
  invisible(x)
}

.check_positive_vector <- function(x, name) {
  if (!is.numeric(x) || anyNA(x) || any(x <= 0)) {
    stop(sprintf("`%s` must contain positive numbers only.", name), call. = FALSE)
  }
  invisible(x)
}

# Log of the multivariate (complex) gamma function
# Gamma_m(a) = pi^{m(m-1)/2} prod_{i=0}^{m-1} Gamma(a - i)
.lmvgamma <- function(a, m) {
  m * (m - 1) / 2 * log(pi) +
    vapply(a, function(ai) sum(lgamma(ai - 0:(m - 1))), numeric(1))
}

# Multivariate polygamma: sum_{i=0}^{m-1} psi^{(deriv)}(a - i)
.mvpsigamma <- function(a, m, deriv = 0L) {
  vapply(a, function(ai) sum(psigamma(ai - 0:(m - 1), deriv = deriv)),
         numeric(1))
}

.logsumexp <- function(x) {
  mx <- max(x)
  if (!is.finite(mx)) return(mx)
  mx + log(sum(exp(x - mx)))
}

# Hermitian matrices --------------------------------------------------------

# number of channels m from the number of distinct entries p = m(m+1)/2
.m_from_p <- function(p) {
  m <- (sqrt(8 * p + 1) - 1) / 2
  if (abs(m - round(m)) > 1e-8) {
    stop("The number of layers is not of the form m(m+1)/2.", call. = FALSE)
  }
  as.integer(round(m))
}

# indices of the upper triangle in the order (1,2), (1,3), ..., (2,3), ...
.upper_pairs <- function(m) {
  if (m < 2) return(matrix(integer(0), ncol = 2))
  pr <- NULL
  for (i in 1:(m - 1)) for (j in (i + 1):m) pr <- rbind(pr, c(i, j))
  pr
}

#' Convert between PolSAR storage formats
#'
#' PolSAR observations are Hermitian positive definite matrices. `SITSAR`
#' accepts them either as a complex array of dimension `m x m x n` or as an
#' `n x p` complex matrix in which each row stores the \eqn{p = m(m+1)/2}
#' distinct entries in the order \eqn{C_{11}, \ldots, C_{mm}} (the
#' intensities) followed by the upper triangle \eqn{C_{12}, C_{13}, \ldots,
#' C_{23}, \ldots}. For \eqn{m = 3} this is the usual layer order
#' `HH, HV, VV, HH.HV*, HH.VV*, HV.VV*` of the [SanFrancisco150] image.
#'
#' @param x A complex array (`m x m x n`), a single `m x m` matrix, or an
#'   `n x p` complex matrix (or vector of length `p`).
#'
#' @return `as_polsar_array()` returns a complex array `m x m x n`;
#'   `as_polsar_vec()` returns an `n x p` complex matrix.
#'
#' @examples
#' S <- diag(3) + 0i
#' S[1, 2] <- 0.2 + 0.1i
#' S[2, 1] <- Conj(S[1, 2])
#' v <- as_polsar_vec(S)
#' v
#' as_polsar_array(v)[, , 1]
#'
#' @name polsar-format
#' @export
as_polsar_array <- function(x) {
  if (is.array(x) && length(dim(x)) == 3L) {
    if (dim(x)[1] != dim(x)[2]) stop("Each slice must be a square matrix.", call. = FALSE)
    return(x)
  }
  if (.is_single_hermitian(x)) return(array(x, dim = c(dim(x), 1L)))
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  p <- ncol(x)
  m <- .m_from_p(p)
  n <- nrow(x)
  pr <- .upper_pairs(m)
  out <- array(0 + 0i, dim = c(m, m, n))
  for (i in seq_len(m)) out[i, i, ] <- Re(x[, i])
  if (m > 1) {
    for (k in seq_len(nrow(pr))) {
      out[pr[k, 1], pr[k, 2], ] <- x[, m + k]
      out[pr[k, 2], pr[k, 1], ] <- Conj(x[, m + k])
    }
  }
  out
}

#' @rdname polsar-format
#' @export
as_polsar_vec <- function(x) {
  if (.is_single_hermitian(x)) x <- array(x, dim = c(dim(x), 1L))
  if (!(is.array(x) && length(dim(x)) == 3L)) {
    if (is.null(dim(x))) return(matrix(x, nrow = 1L))
    return(x)
  }
  m <- dim(x)[1]
  n <- dim(x)[3]
  pr <- .upper_pairs(m)
  out <- matrix(0 + 0i, nrow = n, ncol = m + nrow(pr))
  for (i in seq_len(m)) out[, i] <- x[i, i, ]
  if (m > 1) for (k in seq_len(nrow(pr))) out[, m + k] <- x[pr[k, 1], pr[k, 2], ]
  out
}

# vec (length p) -> Hermitian m x m matrix
.vec_to_herm <- function(v) {
  as_polsar_array(matrix(v, nrow = 1L))[, , 1]
}

# A single m x m Hermitian complex matrix (one PolSAR observation)
.is_single_hermitian <- function(x) {
  is.matrix(x) && is.complex(x) && nrow(x) == ncol(x) && nrow(x) > 1 &&
    max(Mod(x - Conj(t(x)))) <= 1e-10 * max(1, max(Mod(x)))
}

# Make sure input is a vec matrix (n x p)
.as_vec <- function(x) {
  if ((is.array(x) && length(dim(x)) == 3L) || .is_single_hermitian(x)) {
    return(as_polsar_vec(x))
  }
  if (is.null(dim(x))) x <- matrix(x, nrow = 1L)
  x
}

# Real determinants of Hermitian matrices stored as rows of a vec matrix
.hdet_vec <- function(V, m = .m_from_p(ncol(V))) {
  if (m == 1L) return(Re(V[, 1]))
  if (m == 2L) {
    return(Re(V[, 1]) * Re(V[, 2]) - Mod(V[, 3])^2)
  }
  if (m == 3L) {
    a <- Re(V[, 1]); b <- Re(V[, 2]); c <- Re(V[, 3])
    z12 <- V[, 4]; z13 <- V[, 5]; z23 <- V[, 6]
    return(a * b * c - a * Mod(z23)^2 - b * Mod(z13)^2 - c * Mod(z12)^2 +
             2 * Re(z12 * z23 * Conj(z13)))
  }
  A <- as_polsar_array(V)
  vapply(seq_len(dim(A)[3]), function(i) .hdet(A[, , i]), numeric(1))
}

# tr(A Z_i) for a Hermitian matrix A and the rows Z_i of a vec matrix
.htrace_vec <- function(A, V) {
  m <- nrow(A)
  out <- as.vector(Re(V[, seq_len(m), drop = FALSE]) %*% Re(diag(A)))
  if (m > 1) {
    pr <- .upper_pairs(m)
    for (k in seq_len(nrow(pr))) {
      out <- out + 2 * Re(Conj(A[pr[k, 1], pr[k, 2]]) * V[, m + k])
    }
  }
  out
}

.hdet <- function(M) {
  if (!is.complex(M)) return(det(M))
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  prod(Re(ev))
}

.hlogdet <- function(M) {
  if (!is.complex(M)) {
    return(as.numeric(determinant(M, logarithm = TRUE)$modulus))
  }
  ev <- Re(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
  if (any(ev <= 0)) return(NaN)
  sum(log(ev))
}

.hinv <- function(M) {
  out <- solve(M)
  if (is.complex(out)) out <- (out + Conj(t(out))) / 2
  out
}

# Hermitian square root
.hsqrt <- function(M) {
  e <- eigen(M, symmetric = TRUE)
  vals <- pmax(Re(e$values), 0)
  e$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% Conj(t(e$vectors))
}

.check_hpd <- function(Sigma, name = "Sigma") {
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop(sprintf("`%s` must be a square matrix.", name), call. = FALSE)
  }
  if (max(Mod(Sigma - Conj(t(Sigma)))) > 1e-8 * max(1, max(Mod(Sigma)))) {
    stop(sprintf("`%s` must be Hermitian (symmetric).", name), call. = FALSE)
  }
  ev <- Re(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values)
  if (any(ev <= 0)) {
    stop(sprintf("`%s` must be positive definite.", name), call. = FALSE)
  }
  invisible(Sigma)
}

# Small ridge to keep matrices invertible
.regularize <- function(S, eps = 1e-8) {
  d <- mean(Re(diag(S)))
  if (!is.finite(d) || d <= 0) d <- 1
  S + diag(eps * d, nrow(S))
}

# Vectorised bisection for a decreasing function f(L) - s = 0 in log scale.
# f must accept a vector L and return a vector.
.vec_bisect_log <- function(f, s, lower, upper, iter = 80L) {
  n <- length(s)
  lo <- rep(log(lower), n)
  hi <- rep(log(upper), n)
  for (k in seq_len(iter)) {
    mid <- (lo + hi) / 2
    val <- f(exp(mid)) - s
    go_up <- val > 0
    go_up[is.na(go_up)] <- FALSE
    lo[go_up] <- mid[go_up]
    hi[!go_up] <- mid[!go_up]
  }
  out <- exp((lo + hi) / 2)
  out[!is.finite(s)] <- NA_real_
  out
}

# Numerical gradient (central differences)
.num_grad <- function(f, x, ...) {
  h <- 1e-5 * pmax(abs(x), 1e-2)
  vapply(seq_along(x), function(i) {
    xp <- x; xm <- x
    xp[i] <- xp[i] + h[i]
    xm[i] <- xm[i] - h[i]
    (f(xp, ...) - f(xm, ...)) / (2 * h[i])
  }, numeric(1))
}
