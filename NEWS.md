# SITSAR 1.0.0

* Complete rewrite of the package, gathering the methods of the master's
  dissertation and PhD thesis of Jodavid A. Ferreira.
* Intensity distributions: Gamma, generalized Gamma and G0I (exact
  distribution and quantile functions through the F distribution).
* Matrix distributions: scaled complex Wishart, G0 polarimetric, compound
  truncated Poisson and compound geometric complex Wishart.
* R6 classes `GammaSAR`, `G0I`, `GGammaSAR`, `WishartSAR`, `G0Polarimetric` and
  `CompoundWishart` with maximum likelihood / EM estimation, Shannon entropy and
  asymptotic variance of the entropy estimator, plus `tidy()`, `glance()` and
  `logLik()` methods.
* `SARImage` R6 class, `sar_features()`, fast vectorised `entropy_map()` and
  `plot_map()`.
* Entropy-based segmentation (`EntropySegmenter`, `segment_sar()`), hypothesis
  tests (`entropy_test()`) and change detection (`change_map()`).
* Divergences between Gaussian and complex Wishart laws, the
  `DivergenceClassifier` R6 class and the `divergence_classifier()` 'parsnip'
  model for 'tidymodels' workflows.
* Accuracy measures: `kappa_statistics()` and `change_detection_metrics()`.
* New dataset `SanFrancisco150_gt` (ground truth) and new octagonal logo.
* Functions of previous versions (`dg0i()`, `pg0i()`, `qg0i()`, `rg0i()`,
  `g0iestimation()`, `distancesg0i()` and `distancesg0itable()`) are kept.
