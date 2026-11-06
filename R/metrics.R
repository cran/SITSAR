#' Kappa coefficient and overall accuracy
#'
#' Computes Cohen's kappa coefficient of agreement, its asymptotic variance
#' (Congalton and Green, 2009) and the overall accuracy from a confusion
#' matrix or from two vectors (or matrices) of labels. Unlike
#' [yardstick::kap()], the variance of the estimator is also returned, which
#' allows comparing classifiers with Z tests.
#'
#' Labels of unsupervised methods are arbitrary: use `match = TRUE` to
#' relabel the estimated classes so that they best agree with the truth
#' (an exhaustive search over permutations for up to 8 classes).
#'
#' @param truth Either a square confusion matrix (estimates in rows, truth
#'   in columns) or a vector/matrix of reference labels.
#' @param estimate Vector/matrix of estimated labels (when `truth` is not a
#'   confusion matrix).
#' @param match Relabel the estimated classes to maximise agreement?
#'
#' @return A tibble with `kappa`, `variance`, `std_error`,
#'   `overall_accuracy` and the number of observations `n`.
#'
#' @references
#' Congalton, R. G. and Green, K. (2009). *Assessing the Accuracy of
#' Remotely Sensed Data: Principles and Practices*. CRC Press.
#'
#' @examples
#' truth <- c(1, 1, 2, 2, 3, 3, 3)
#' est <- c(1, 2, 2, 2, 3, 3, 1)
#' kappa_statistics(truth, est)
#' kappa_statistics(table(est, truth))
#'
#' @export
kappa_statistics <- function(truth, estimate = NULL, match = FALSE) {
  if (is.null(estimate)) {
    tab <- as.matrix(truth)
    if (nrow(tab) != ncol(tab)) stop("The confusion matrix must be square.", call. = FALSE)
  } else {
    tv <- as.vector(truth); ev <- as.vector(estimate)
    if (length(tv) != length(ev)) stop("`truth` and `estimate` differ in length.", call. = FALSE)
    ok <- !is.na(tv) & !is.na(ev)
    lev <- sort(unique(c(as.character(tv[ok]), as.character(ev[ok]))))
    tab <- table(factor(ev[ok], levels = lev), factor(tv[ok], levels = lev))
    tab <- unclass(tab)
    if (match) tab <- .match_labels(tab)
  }
  tab <- tab / 1
  N <- sum(tab)
  P <- tab / N
  r <- rowSums(P); cc <- colSums(P)
  t1 <- sum(diag(P))
  t2 <- sum(r * cc)
  t3 <- sum(diag(P) * (r + cc))
  # theta_4 = sum_ij p_ij (p_{j+} + p_{+i})^2
  t4 <- sum(P * outer(cc, r, `+`)^2)
  kap <- (t1 - t2) / (1 - t2)
  v <- (t1 * (1 - t1) / (1 - t2)^2 +
          2 * (1 - t1) * (2 * t1 * t2 - t3) / (1 - t2)^3 +
          (1 - t1)^2 * (t4 - 4 * t2^2) / (1 - t2)^4) / N
  tibble::tibble(kappa = kap, variance = v, std_error = sqrt(v),
                 overall_accuracy = t1, n = N)
}

.match_labels <- function(tab) {
  k <- nrow(tab)
  if (k > 8) return(tab)
  perms <- .permutations(k)
  score <- apply(perms, 1, function(pm) sum(tab[cbind(pm, seq_len(k))]))
  best <- perms[which.max(score), ]
  out <- tab[best, , drop = FALSE]
  rownames(out) <- colnames(tab)
  out
}

.permutations <- function(k) {
  if (k == 1) return(matrix(1L, 1, 1))
  sub <- .permutations(k - 1)
  do.call(rbind, lapply(seq_len(k), function(i) {
    cbind(i, matrix(setdiff(seq_len(k), i)[sub], nrow(sub)))
  }))
}

#' Accuracy measures for change detection
#'
#' Computes the usual accuracy measures of binary change maps: detection
#' rate (DR, percentage of true changes that were detected), false alarm
#' rate (FA, percentage of unchanged pixels flagged as changed), and the
#' percentages of false positives (FP) and false negatives (FN) with
#' respect to all pixels.
#'
#' @param truth Logical (or 0/1) vector/matrix of true changes.
#' @param estimate Logical (or 0/1) vector/matrix of detected changes.
#'
#' @return A tibble with `DR`, `FA`, `FP`, `FN` (in percent) and the
#'   overall accuracy `OA` (in percent).
#'
#' @examples
#' truth <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
#' est <- c(TRUE, FALSE, TRUE, FALSE, FALSE)
#' change_detection_metrics(truth, est)
#'
#' @export
change_detection_metrics <- function(truth, estimate) {
  tv <- as.logical(as.vector(truth)); ev <- as.logical(as.vector(estimate))
  if (length(tv) != length(ev)) stop("`truth` and `estimate` differ in length.", call. = FALSE)
  ok <- !is.na(tv) & !is.na(ev)
  tv <- tv[ok]; ev <- ev[ok]
  tp <- sum(tv & ev); tn <- sum(!tv & !ev)
  fp <- sum(!tv & ev); fn <- sum(tv & !ev)
  N <- length(tv)
  tibble::tibble(DR = 100 * tp / max(tp + fn, 1),
                 FA = 100 * fp / max(fp + tn, 1),
                 FP = 100 * fp / N,
                 FN = 100 * fn / N,
                 OA = 100 * (tp + tn) / N)
}
