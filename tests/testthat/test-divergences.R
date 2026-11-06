test_that("Gaussian divergences", {
  m1 <- c(0, 1, 2); m2 <- c(1, 1, 0)
  S1 <- diag(c(1, 2, 3)); S2 <- matrix(c(2, 0.3, 0, 0.3, 1, 0.2, 0, 0.2, 1.5), 3)
  kl <- divergence_gaussian(m1, S1, m2, S2, "kl")
  ref <- 0.5 * (sum(diag(solve(S2) %*% S1)) + t(m2 - m1) %*% solve(S2) %*% (m2 - m1) -
                  3 + log(det(S2) / det(S1)))
  expect_equal(kl, as.numeric(ref))
  expect_equal(divergence_gaussian(m1, S1, m1, S1, "kl"), 0)
  expect_equal(divergence_gaussian(m1, S1, m2, S2, "jeffreys"),
               kl + divergence_gaussian(m2, S2, m1, S1, "kl"))
  # Renyi -> KL when the order -> 1
  expect_equal(divergence_gaussian(m1, S1, m2, S2, "renyi", order = 1 - 1e-6), kl,
               tolerance = 1e-4)
  # Bhattacharyya = Renyi(1/2) / 4 in Pardo's normalisation
  expect_equal(divergence_gaussian(m1, S1, m2, S2, "bhattacharyya"),
               divergence_gaussian(m1, S1, m2, S2, "renyi", order = 0.5) / 4)
  r1 <- divergence_gaussian(m1, S1, m2, S2, "renyi", 0.3, kind = 1)
  r2 <- divergence_gaussian(m1, S1, m2, S2, "renyi", 0.3, kind = 2)
  expect_equal(divergence_gaussian(m1, S1, m2, S2, "renyi", 0.3, kind = 3), (r1 + r2) / 2)
  expect_gt(divergence_gaussian(m1, S1, m2, S2, "renyi", 0.3, kind = 4), 0)
  expect_error(divergence_gaussian(m1, S1, m2, S2, "renyi", order = 1.5))
})

test_that("Wishart divergences", {
  S1 <- polsar_example_sigma("forest"); S2 <- polsar_example_sigma("urban")
  kl <- divergence_wishart(S1, S2, 4, "kl")
  expect_gt(kl, 0)
  expect_equal(divergence_wishart(S1, S1, 4, "jeffreys"), 0, tolerance = 1e-10)
  expect_equal(divergence_wishart(S1, S2, 4, "jeffreys"),
               kl + divergence_wishart(S2, S1, 4, "kl"))
  expect_equal(divergence_wishart(S1, S2, 4, "renyi", order = 1 - 1e-6), kl,
               tolerance = 1e-4)
  # Monte Carlo check of the KL divergence
  set.seed(50)
  Z <- rwishart_sar(20000, S1, 4)
  mc <- mean(dwishart_sar(Z, S1, 4, log = TRUE) - dwishart_sar(Z, S2, 4, log = TRUE))
  expect_equal(kl, mc, tolerance = 0.02)
  # Bhattacharyya distance: -log of the Bhattacharyya coefficient
  bc <- mean(exp(0.5 * (dwishart_sar(Z, S2, 4, log = TRUE) - dwishart_sar(Z, S1, 4, log = TRUE))))
  expect_equal(divergence_wishart(S1, S2, 4, "bhattacharyya"), -log(bc), tolerance = 0.05)
})

test_that("G0I stochastic distances (compatibility functions)", {
  d0 <- distancesg0i(c(-1.5, 2, 4), c(-1.5, 2, 4))
  expect_s3_class(d0, "distancesg0i")
  expect_equal(d0$distance_value, 0, tolerance = 1e-8)
  for (tp in c("KL", "Renyi", "Bhattacharyya", "Hellinger", "Arithmetic-Geometric",
               "Triangular", "Harmonic-Mean")) {
    expect_gt(distancesg0i(c(-1.5, 2, 8), c(-4, 2, 4), type = tp)$distance_value, 0)
  }
  expect_output(print(d0), "Distance")
  expect_error(distancesg0i(c(1, 2, 4), c(-1.5, 2, 4)))
  tab <- distancesg0itable(SanFrancisco150, c(10, 25, 10, 25), c(10, 25, 125, 140), L = 4)
  expect_equal(dim(tab$distancestable), c(3L, 7L))
  expect_true(all(tab$distancestable >= 0))
  expect_output(print(tab), "distances")
})
