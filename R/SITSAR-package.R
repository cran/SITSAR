#' SITSAR: Information Theory, Geometry and Statistics for SAR Data
#'
#' Statistical models, information-theoretic measures and image analysis
#' tools for Synthetic Aperture Radar (SAR) and polarimetric SAR (PolSAR)
#' data, gathering the methods developed in the master's dissertation
#' (divergence-based classifiers) and the PhD thesis (Shannon entropies,
#' entropy-based segmentation and change detection, and compound Wishart
#' models) of Jodavid A. Ferreira, supervised by Abraão D. C. Nascimento.
#'
#' @section Main tools:
#' \describe{
#'   \item{Distributions}{[dg0i], [dgammasar], [dggamma], [dwishart_sar],
#'     [dg0p], [compound_wishart].}
#'   \item{R6 models}{[GammaSAR], [G0I], [GGammaSAR], [WishartSAR],
#'     [G0Polarimetric], [CompoundWishart] (estimation, entropy and
#'     asymptotic variance).}
#'   \item{Images}{[SARImage], [sar_features()], [entropy_map()],
#'     [plot_map()].}
#'   \item{Segmentation and change detection}{[EntropySegmenter],
#'     [segment_sar()], [entropy_test()], [change_map()].}
#'   \item{Classification}{[divergences], [DivergenceClassifier] and the
#'     'parsnip' model [divergence_classifier()].}
#'   \item{Accuracy}{[kappa_statistics()], [change_detection_metrics()].}
#' }
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom R6 R6Class
#' @importFrom utils globalVariables
## usethis namespace: end
NULL

globalVariables(c("object", "new_data"))
