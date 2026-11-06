test_that("kappa statistics", {
  tab <- matrix(c(30, 5, 2, 4, 40, 6, 1, 3, 25), 3)
  ks <- kappa_statistics(tab)
  po <- sum(diag(tab)) / sum(tab)
  pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  expect_equal(ks$kappa, (po - pe) / (1 - pe))
  expect_equal(ks$overall_accuracy, po)
  expect_gt(ks$variance, 0)
  truth <- c(1, 1, 2, 2, 3, 3)
  expect_equal(kappa_statistics(truth, truth)$kappa, 1)
  swapped <- c(2, 2, 3, 3, 1, 1)
  expect_equal(kappa_statistics(truth, swapped, match = TRUE)$kappa, 1)
  expect_lt(kappa_statistics(truth, swapped)$kappa, 0)
  skip_if_not_installed("yardstick")
  set.seed(70)
  a <- factor(sample(1:3, 200, TRUE)); b <- a
  b[1:60] <- factor(sample(1:3, 60, TRUE), levels = levels(a))
  expect_equal(kappa_statistics(a, b)$kappa,
               yardstick::kap_vec(a, b), tolerance = 1e-10)
})

test_that("change detection metrics", {
  m <- change_detection_metrics(c(TRUE, TRUE, FALSE, FALSE), c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(m$DR, 50)
  expect_equal(m$FA, 0)
  expect_equal(m$FN, 25)
  expect_equal(m$OA, 75)
})
