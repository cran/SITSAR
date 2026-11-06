make_data <- function() {
  img <- SARImage$new(SanFrancisco150)
  f <- sar_features(img, window = 5)
  f$class <- factor(as.vector(SanFrancisco150_gt), labels = c("ocean", "forest", "urban"))
  set.seed(60)
  idx <- unlist(lapply(split(seq_len(nrow(f)), f$class), sample, 100))
  list(train = f[idx, ], test = f[sample(nrow(f), 600), ])
}

test_that("R6 divergence classifier (Wishart and Gaussian)", {
  d <- make_data()
  for (fam in c("wishart", "gaussian")) {
    for (dv in c("kl", "jeffreys", "bhattacharyya", "renyi")) {
      clf <- DivergenceClassifier$new(family = fam, divergence = dv, order = 0.5, kind = 3)
      clf$fit(d$train, d$train$class)
      pred <- clf$predict(d$test)
      expect_s3_class(pred, "factor")
      expect_equal(levels(pred), c("ocean", "forest", "urban"))
      expect_gt(mean(pred == d$test$class), 0.7)
    }
  }
  pr <- clf$predict(d$test, type = "prob")
  expect_equal(names(pr), c("ocean", "forest", "urban"))
  expect_equal(rowSums(pr), rep(1, nrow(pr)))
  ds <- clf$predict(d$test, type = "distance")
  expect_equal(dim(ds), c(600L, 3L))
  expect_output(print(clf), "DivergenceClassifier")
  expect_error(DivergenceClassifier$new("wishart")$fit(d$train[, 1:5], d$train$class),
               "Wishart features")
})

test_that("parsnip model specification", {
  skip_if_not_installed("parsnip")
  d <- make_data()
  spec <- divergence_classifier(divergence = "renyi", order = 0.5, kind = 1) |>
    parsnip::set_engine("wishart")
  expect_s3_class(spec, "divergence_classifier")
  expect_output(print(spec), "Divergence Classifier")
  fit <- parsnip::fit(spec, class ~ ., data = d$train)
  p <- predict(fit, d$test)
  expect_named(p, ".pred_class")
  expect_gt(mean(p$.pred_class == d$test$class), 0.7)
  pp <- predict(fit, d$test, type = "prob")
  expect_named(pp, c(".pred_ocean", ".pred_forest", ".pred_urban"))
  spec2 <- update(spec, order = 0.8)
  expect_equal(rlang::eval_tidy(spec2$args$order), 0.8)
  g <- parsnip::fit(divergence_classifier(divergence = "kl"), class ~ ., data = d$train)
  expect_gt(mean(predict(g, d$test)$.pred_class == d$test$class), 0.7)
})

test_that("tidymodels workflow and tuning", {
  skip_if_not_installed("workflows")
  skip_if_not_installed("recipes")
  skip_if_not_installed("rsample")
  skip_if_not_installed("tune")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("dials")
  d <- make_data()
  rec <- recipes::recipe(class ~ ., data = d$train) |>
    recipes::update_role(row, col, new_role = "id")
  spec <- divergence_classifier(divergence = "renyi", order = tune::tune()) |>
    parsnip::set_engine("wishart")
  wf <- workflows::workflow(rec, spec)
  set.seed(61)
  folds <- rsample::vfold_cv(d$train, v = 3, strata = class)
  grid <- data.frame(order = c(0.2, 0.8))
  res <- tune::tune_grid(wf, resamples = folds, grid = grid,
                         metrics = yardstick::metric_set(yardstick::accuracy))
  best <- tune::select_best(res, metric = "accuracy")
  final <- tune::finalize_workflow(wf, best) |> parsnip::fit(d$train)
  acc <- yardstick::accuracy(
    dplyr_free_bind(d$test$class, predict(final, d$test)$.pred_class),
    truth, estimate)
  expect_gt(acc$.estimate, 0.7)
  expect_s3_class(renyi_order(), "quant_param")
  expect_s3_class(divergence_type(), "qual_param")
})
