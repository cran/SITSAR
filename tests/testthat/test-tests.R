test_that("entropy test has the right size and power", {
  set.seed(30)
  p_null <- replicate(100, entropy_test(rgammasar(200, 1, 4), rgammasar(200, 1, 4),
                                        model = "gamma", L = 4)$p.value)
  expect_lt(abs(mean(p_null < 0.05) - 0.05), 0.06)
  tt <- entropy_test(rg0i(300, -3, 2, 4), rg0i(300, -1.5, 1, 4), model = "g0i", L = 4)
  expect_s3_class(tt, "htest")
  expect_lt(tt$p.value, 0.01)
  S <- polsar_example_sigma("forest")
  tw <- entropy_test(rwishart_sar(200, S, 4), rwishart_sar(200, 2 * S, 4), "wishart", L = 4)
  expect_lt(tw$p.value, 1e-6)
})

test_that("change maps detect a simulated change", {
  set.seed(31)
  a <- matrix(rgammasar(40 * 40, 1, 4), 40)
  b <- matrix(rgammasar(40 * 40, 1, 4), 40)
  b[11:30, 11:30] <- rgammasar(400, 5, 4)
  cd <- change_map(a, b, model = "gamma", window = 5, L = 4)
  expect_s3_class(cd, "sitsar_change")
  truth <- matrix(FALSE, 40, 40); truth[11:30, 11:30] <- TRUE
  met <- change_detection_metrics(truth, cd$change)
  expect_gt(met$DR, 90)
  expect_lt(met$FA, 15)
  expect_output(print(cd), "changed pixels")
})
