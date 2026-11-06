## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4,
                      fig.align = "center")

## ----setup--------------------------------------------------------------------
library(SITSAR)

## -----------------------------------------------------------------------------
set.seed(1)
z <- rg0i(5, alpha = -3, gamma = 2, L = 4)
dg0i(z, alpha = -3, gamma = 2, L = 4)
pg0i(1, alpha = -3, gamma = 2, L = 4)
qg0i(0.5, alpha = -3, gamma = 2, L = 4)

## -----------------------------------------------------------------------------
Sigma <- polsar_example_sigma("forest")
Z <- rwishart_sar(3, Sigma, L = 4)
dwishart_sar(Z, Sigma, L = 4, log = TRUE)
dg0p(Z, alpha = -3, Sigma = Sigma, L = 4, log = TRUE)
dctpcw(Z, lambda = 1.5, Sigma = Sigma, L = 4, log = TRUE)

## -----------------------------------------------------------------------------
set.seed(2)
x <- rg0i(800, alpha = -3, gamma = 2, L = 4)
m <- G0I$new()$fit(x)          # alpha, gamma and L estimated
m
m$entropy_variance()
generics::tidy(m)
generics::glance(m)

## -----------------------------------------------------------------------------
entropy_g0i(alpha = -3, gamma = 2, L = 4)
entropy_wishart(Sigma, L = 4)
entropy_g0p(alpha = -3, Sigma = Sigma, L = 4)

## -----------------------------------------------------------------------------
set.seed(3)
S <- rctpcw(300, lambda = 2, Sigma = Sigma, L = 4)
cw <- CompoundWishart$new("poisson")$fit(S, L = 4)
cw$params$lambda
cw$iterations

## ----fig.show = "hold", out.width = "48%", fig.width = 4, fig.height = 4------
img <- SARImage$new(SanFrancisco150, looks = 4)
img
img$plot()
plot_map(SanFrancisco150_gt, main = "Ground truth",
         palette = c("#1f78b4", "#33a02c", "#e31a1c"))

## ----fig.show = "hold", out.width = "48%", fig.width = 4, fig.height = 4------
Hw <- entropy_map(img, model = "wishart", window = 7, L = 4)
plot_map(Hw, main = "Wishart entropy (L = 4)")
Hhv <- entropy_map(img, model = "gamma", window = 7, channels = "HV")
plot_map(Hhv, main = "Gamma entropy, HV (L estimated)")

## ----fig.width = 4, fig.height = 4--------------------------------------------
set.seed(4)
seg <- EntropySegmenter$new(k = 3, method = "em")$fit(Hw)
seg
plot_map(seg$label_map(), main = "Segmentation (Wishart entropy)",
         palette = c("#1f78b4", "#33a02c", "#e31a1c"))
kappa_statistics(SanFrancisco150_gt, seg$label_map(), match = TRUE)

## -----------------------------------------------------------------------------
set.seed(5)
entropy_test(rg0i(300, -3, 2, 4), rg0i(300, -1.5, 1, 4), model = "g0i", L = 4)

## ----fig.width = 4, fig.height = 4--------------------------------------------
set.seed(6)
a <- matrix(rgammasar(60 * 60, mu = 1, L = 4), 60)
b <- matrix(rgammasar(60 * 60, mu = 1, L = 4), 60)
b[21:40, 21:45] <- rgammasar(500, mu = 3, L = 4)
cd <- change_map(a, b, model = "gamma", window = 5, L = 4, level = 0.01)
cd
plot_map(cd$change, main = "Detected changes")
truth <- matrix(FALSE, 60, 60); truth[21:40, 21:45] <- TRUE
change_detection_metrics(truth, cd$change)

## -----------------------------------------------------------------------------
divergence_wishart(polsar_example_sigma("forest"), polsar_example_sigma("urban"),
                   L = 4, type = "renyi", order = 0.5, kind = 3)

