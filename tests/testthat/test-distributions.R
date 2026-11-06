test_that("G0I density, cdf, quantile and generator are consistent", {
  expect_equal(integrate(function(x) dg0i(x, -3, 2, 4), 0, Inf)$value, 1, tolerance = 1e-6)
  q <- c(0.2, 1, 3)
  num <- vapply(q, function(u) integrate(function(x) dg0i(x, -3, 2, 4), 0, u)$value, 1)
  expect_equal(pg0i(q, -3, 2, 4), num, tolerance = 1e-6)
  p <- c(0.1, 0.5, 0.9)
  expect_equal(pg0i(qg0i(p, -2, 1, 3), -2, 1, 3), p)
  expect_equal(pg0i(1, -2, 1, 3, lower.tail = FALSE), 1 - pg0i(1, -2, 1, 3))
  expect_equal(dg0i(c(0.5, 2), -2, 1, 3, log = TRUE), log(dg0i(c(0.5, 2), -2, 1, 3)))
  set.seed(1)
  z <- rg0i(1e5, -4, 3, 4)
  expect_equal(mean(z), 3 / (4 - 1), tolerance = 0.02)
  expect_equal(dg0i(c(-1, 0), -2, 1, 1), c(0, 0))
  expect_error(dg0i(1, 2, 1, 1), "negative")
  expect_error(rg0i(1, -2, -1, 1), "positive")
})

test_that("Gamma and generalized Gamma laws", {
  expect_equal(dgammasar(1.5, mu = 2, L = 3), dgamma(1.5, 3, 1.5))
  expect_equal(pgammasar(qgammasar(0.3, 2, 3), 2, 3), 0.3)
  expect_equal(integrate(function(x) dggamma(x, 1.7, 2, 0.8), 0, Inf)$value, 1, tolerance = 1e-6)
  expect_equal(dggamma(1.2, 1, 3, 2), dgamma(1.2, 3, 2))
  expect_equal(pggamma(qggamma(0.7, 1.3, 2, 1), 1.3, 2, 1), 0.7)
  set.seed(2)
  expect_equal(mean(rgammasar(1e5, 2, 4)), 2, tolerance = 0.01)
})

test_that("Wishart and G0 polarimetric laws", {
  S <- polsar_example_sigma("forest")
  set.seed(3)
  W <- rwishart_sar(5000, S, L = 4, output = "vec")
  expect_equal(dim(W), c(5000L, 6L))
  expect_lt(max(Mod(colMeans(W) - as_polsar_vec(S))), 0.01)
  # all simulated matrices are positive definite
  expect_true(all(SITSAR:::.hdet_vec(W) > 0))
  # m = 1 reduces to the Gamma / G0I laws
  expect_equal(dwishart_sar(matrix(0.7 + 0i), matrix(2 + 0i), 3), dgammasar(0.7, 2, 3))
  expect_equal(dg0p(matrix(0.7 + 0i), -3, matrix(1 + 0i), 4), dg0i(0.7, -3, 2, 4))
  # G0P converges to the Wishart law
  Z <- rwishart_sar(3, S, 4)
  expect_equal(dg0p(Z, -1e7, S, 4, log = TRUE), dwishart_sar(Z, S, 4, log = TRUE),
               tolerance = 1e-4)
  G <- rg0p(4000, -5, S, 4, output = "vec")
  expect_lt(max(Mod(colMeans(G) - as_polsar_vec(S))), 0.02)
  expect_error(dwishart_sar(Z, diag(3) - 2, 4), "positive definite")
})

test_that("compound Wishart laws", {
  S <- polsar_example_sigma("pasture")
  set.seed(4)
  x <- rctpcw(3000, lambda = 2, Sigma = S, L = 4, output = "vec")
  expect_equal(mean(attr(x, "N")), 2 / (1 - exp(-2)), tolerance = 0.05)
  y <- rcgcw(3000, p = 0.5, Sigma = S, L = 4)
  expect_equal(mean(attr(y, "N")), 2, tolerance = 0.05)
  # with lambda -> 0 / p -> 1 the CTPCW / CGCW reduce to an unscaled Wishart
  Z <- rwishart_sar(3, S, 4) * 4
  ref <- dwishart_sar(Z / 4, S, 4, log = TRUE) - 9 * log(4)
  expect_equal(dctpcw(Z, 1e-8, S, 4, log = TRUE), ref, tolerance = 1e-6)
  expect_equal(dcgcw(Z, 1 - 1e-10, S, 4, log = TRUE), ref, tolerance = 1e-6)
})

test_that("PolSAR storage conversions round trip", {
  v <- SanFrancisco150[1, 1, ]
  A <- as_polsar_array(matrix(v, 1))
  expect_equal(dim(A), c(3L, 3L, 1L))
  expect_equal(A[, , 1], Conj(t(A[, , 1])))
  expect_equal(as.vector(as_polsar_vec(A)), v)
  expect_equal(as_polsar_vec(A[, , 1]), matrix(v, 1))
})
