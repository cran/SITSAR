#' San Francisco PolSAR image (150 x 150 pixels)
#'
#' A \eqn{150 \times 150} crop of the fully polarimetric L-band AIRSAR image
#' of the San Francisco Bay area (nominal 4 looks), used in the
#' dissertation and thesis of J. A. Ferreira. The scene contains ocean,
#' vegetation (forest/park) and urban areas; see [SanFrancisco150_gt].
#'
#' @format A complex array of dimension `150 x 150 x 6`. For each pixel the
#'   six layers are the distinct entries of the \eqn{3 \times 3} Hermitian
#'   covariance matrix: \eqn{|S_{HH}|^2}, \eqn{|S_{HV}|^2}, \eqn{|S_{VV}|^2}
#'   (real intensities), \eqn{S_{HH}S_{HV}^*}, \eqn{S_{HH}S_{VV}^*} and
#'   \eqn{S_{HV}S_{VV}^*} (complex products).
#'
#' @source NASA/JPL AIRSAR mission.
#'
#' @examples
#' dim(SanFrancisco150)
#' img <- SARImage$new(SanFrancisco150, looks = 4)
#' img$plot()
"SanFrancisco150"

#' Ground truth of the San Francisco image
#'
#' Reference map of the [SanFrancisco150] image with three classes.
#'
#' @format An integer matrix `150 x 150` with values `1` (ocean), `2`
#'   (vegetation) and `3` (urban area).
#'
#' @examples
#' table(SanFrancisco150_gt)
#' plot_map(SanFrancisco150_gt)
"SanFrancisco150_gt"
