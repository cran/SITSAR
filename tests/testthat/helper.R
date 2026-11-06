dplyr_free_bind <- function(truth, estimate) {
  data.frame(truth = truth, estimate = factor(estimate, levels = levels(truth)))
}
