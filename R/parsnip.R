#' Divergence-based classifier as a 'parsnip' model
#'
#' `divergence_classifier()` defines the divergence-based classifier of
#' Ferreira, Coelho and Nascimento (2021) as a 'parsnip' model
#' specification, so it can be used with 'workflows', 'recipes', 'rsample',
#' 'tune' and 'yardstick'. The engine selects the family of features:
#' `"gaussian"` (multivariate normal intensity features, the default) or
#' `"wishart"` (complex Wishart PolSAR matrix features). The predictors must
#' be in the format of [sar_features()]; extra columns such as `row`/`col`
#' are ignored by the engine (it is good practice to give them an `"id"`
#' role in a recipe).
#'
#' @param mode A single character string for the type of model. The only
#'   possible value is `"classification"`.
#' @param divergence Divergence: `"renyi"` (default), `"kl"`, `"jeffreys"`
#'   or `"bhattacharyya"`.
#' @param order Order of the Rényi divergence, in (0, 1). Tunable with
#'   [renyi_order()].
#' @param kind Form of the Rényi divergence (1 to 4); see [divergences].
#' @param engine `"gaussian"` or `"wishart"`.
#'
#' @return A model specification object with classes
#'   `divergence_classifier` and `model_spec`.
#'
#' @seealso [DivergenceClassifier], [sar_features()]
#'
#' @examples
#' spec <- divergence_classifier(divergence = "renyi", order = 0.5) |>
#'   parsnip::set_engine("wishart")
#' spec
#'
#' feats <- sar_features(SARImage$new(SanFrancisco150)$crop(1:60, 1:60), window = 5)
#' feats$class <- factor(as.vector(SanFrancisco150_gt[1:60, 1:60]))
#' fit <- parsnip::fit(spec, class ~ ., data = feats)
#' head(predict(fit, feats))
#' head(predict(fit, feats, type = "prob"))
#'
#' @export
divergence_classifier <- function(mode = "classification", divergence = NULL,
                                  order = NULL, kind = NULL, engine = "gaussian") {
  if (mode != "classification") {
    stop("`mode` should be 'classification'.", call. = FALSE)
  }
  args <- list(divergence = rlang::enquo(divergence),
               order = rlang::enquo(order),
               kind = rlang::enquo(kind))
  parsnip::new_model_spec(
    "divergence_classifier",
    args = args, eng_args = NULL, mode = mode,
    user_specified_mode = !missing(mode),
    method = NULL, engine = engine,
    user_specified_engine = !missing(engine)
  )
}

#' @param object A `divergence_classifier` model specification.
#' @param parameters A 1-row tibble or named list with the main parameters.
#' @param fresh Replace the arguments instead of updating them?
#' @param ... Not used.
#' @rdname divergence_classifier
#' @importFrom stats update
#' @export
update.divergence_classifier <- function(object, parameters = NULL, divergence = NULL,
                                         order = NULL, kind = NULL, fresh = FALSE, ...) {
  args <- list(divergence = rlang::enquo(divergence),
               order = rlang::enquo(order),
               kind = rlang::enquo(kind))
  parsnip::update_spec(object = object, parameters = parameters, args_enquo_list = args,
                       fresh = fresh, cls = "divergence_classifier", ...)
}

#' Engine fit function for `divergence_classifier()`
#'
#' Internal fit function used by the 'parsnip' engines; it is exported only
#' so that 'parsnip' can call it.
#'
#' @param x Data frame of predictors.
#' @param y Factor of classes.
#' @param family `"gaussian"` or `"wishart"`.
#' @param divergence,order,kind See [divergence_classifier()].
#'
#' @return A fitted [DivergenceClassifier].
#' @keywords internal
#' @export
sitsar_divergence_fit <- function(x, y, family = "gaussian", divergence = "renyi",
                                  order = 0.5, kind = 1) {
  DivergenceClassifier$new(family = family, divergence = divergence,
                           order = order, kind = kind)$fit(x, y)
}

#' Tuning parameters for the divergence classifier
#'
#' 'dials' parameter objects for [divergence_classifier()]:
#' `renyi_order()` for the order of the Rényi divergence and
#' `divergence_type()` for the divergence measure. Both require the 'dials'
#' package.
#'
#' @param range Range of the Rényi order.
#' @param trans A transformation object (or `NULL`).
#' @param values Candidate divergences.
#'
#' @return A 'dials' parameter object.
#'
#' @examples
#' if (requireNamespace("dials", quietly = TRUE)) {
#'   renyi_order()
#'   divergence_type()
#' }
#'
#' @export
renyi_order <- function(range = c(0.05, 0.95), trans = NULL) {
  rlang::check_installed("dials")
  dials::new_quant_param(type = "double", range = range, inclusive = c(TRUE, TRUE),
                         trans = trans, label = c(order = "Renyi Order"),
                         finalize = NULL)
}

#' @rdname renyi_order
#' @export
divergence_type <- function(values = c("renyi", "kl", "jeffreys", "bhattacharyya")) {
  rlang::check_installed("dials")
  dials::new_qual_param(type = "character", values = values,
                        label = c(divergence = "Divergence"), finalize = NULL)
}

#' @rdname renyi_order
#' @export
renyi_kind <- function(range = c(1L, 4L)) {
  rlang::check_installed("dials")
  dials::new_quant_param(type = "integer", range = range, inclusive = c(TRUE, TRUE),
                         trans = NULL, label = c(kind = "Renyi Form"), finalize = NULL)
}

.make_divergence_classifier <- function() {
  current <- parsnip::get_model_env()
  if ("divergence_classifier" %in% current$models) return(invisible(FALSE))
  parsnip::set_new_model("divergence_classifier")
  parsnip::set_model_mode("divergence_classifier", "classification")
  for (eng in c("gaussian", "wishart")) {
    parsnip::set_model_engine("divergence_classifier", "classification", eng)
    parsnip::set_dependency("divergence_classifier", eng, "SITSAR",
                            mode = "classification")
    parsnip::set_model_arg(
      model = "divergence_classifier", eng = eng, parsnip = "divergence",
      original = "divergence", func = list(pkg = "SITSAR", fun = "divergence_type"),
      has_submodel = FALSE)
    parsnip::set_model_arg(
      model = "divergence_classifier", eng = eng, parsnip = "order",
      original = "order", func = list(pkg = "SITSAR", fun = "renyi_order"),
      has_submodel = FALSE)
    parsnip::set_model_arg(
      model = "divergence_classifier", eng = eng, parsnip = "kind",
      original = "kind", func = list(pkg = "SITSAR", fun = "renyi_kind"),
      has_submodel = FALSE)
    parsnip::set_fit(
      model = "divergence_classifier", eng = eng, mode = "classification",
      value = list(interface = "data.frame", protect = c("x", "y"),
                   func = c(pkg = "SITSAR", fun = "sitsar_divergence_fit"),
                   defaults = list(family = eng)))
    parsnip::set_encoding(
      model = "divergence_classifier", eng = eng, mode = "classification",
      options = list(predictor_indicators = "none", compute_intercept = FALSE,
                     remove_intercept = FALSE, allow_sparse_x = FALSE))
    parsnip::set_pred(
      model = "divergence_classifier", eng = eng, mode = "classification",
      type = "class",
      value = list(pre = NULL, post = NULL, func = c(fun = "predict"),
                   args = list(object = rlang::expr(object$fit),
                               new_data = rlang::expr(new_data),
                               type = "class")))
    parsnip::set_pred(
      model = "divergence_classifier", eng = eng, mode = "classification",
      type = "prob",
      value = list(pre = NULL, post = NULL, func = c(fun = "predict"),
                   args = list(object = rlang::expr(object$fit),
                               new_data = rlang::expr(new_data),
                               type = "prob")))
  }
  invisible(TRUE)
}

#' @export
print.divergence_classifier <- function(x, ...) {
  parsnip::print_model_spec(x, desc = "Divergence Classifier", ...)
  invisible(x)
}
