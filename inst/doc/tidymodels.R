## ----include = FALSE----------------------------------------------------------
pkgs <- c("parsnip", "recipes", "workflows", "rsample", "tune", "yardstick", "dials")
has_pkgs <- all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = has_pkgs,
                      fig.width = 4, fig.height = 4, fig.align = "center")

## ----setup, message = FALSE---------------------------------------------------
library(SITSAR)
library(parsnip)
library(recipes)
library(workflows)
library(rsample)
library(tune)
library(yardstick)

## -----------------------------------------------------------------------------
img <- SARImage$new(SanFrancisco150, looks = 4)
feats <- sar_features(img, window = 5)
feats$class <- factor(as.vector(SanFrancisco150_gt),
                      labels = c("ocean", "vegetation", "urban"))
dim(feats)

set.seed(123)
split <- initial_split(feats, prop = 0.02, strata = class)
train <- training(split)
test <- testing(split)
nrow(train)

## -----------------------------------------------------------------------------
spec <- divergence_classifier(divergence = "renyi", order = 0.5, kind = 3) |>
  set_engine("wishart")
spec

rec <- recipe(class ~ ., data = train) |>
  update_role(row, col, new_role = "id")

wf <- workflow(rec, spec)
fitted <- fit(wf, data = train)
pred <- augment(fitted, test)
accuracy(pred, truth = class, estimate = .pred_class)
conf_mat(pred, truth = class, estimate = .pred_class)

## -----------------------------------------------------------------------------
spec_tune <- divergence_classifier(divergence = "renyi", order = tune()) |>
  set_engine("wishart")
set.seed(1)
folds <- vfold_cv(train, v = 3, strata = class)
res <- tune_grid(workflow(rec, spec_tune), resamples = folds,
                 grid = data.frame(order = c(0.1, 0.5, 0.9)),
                 metrics = metric_set(accuracy, kap))
collect_metrics(res)
best <- select_best(res, metric = "accuracy")
best

## -----------------------------------------------------------------------------
final <- finalize_workflow(workflow(rec, spec_tune), best) |> fit(train)
map <- predict(final, feats)$.pred_class
plot_map(matrix(as.integer(map), 150), main = "Rényi-Wishart classifier",
         palette = c("#1f78b4", "#33a02c", "#e31a1c"))
kappa_statistics(SanFrancisco150_gt, matrix(as.integer(map), 150))

## -----------------------------------------------------------------------------
grid <- expand.grid(engine = c("gaussian", "wishart"),
                    divergence = c("kl", "jeffreys", "bhattacharyya", "renyi"),
                    stringsAsFactors = FALSE)
grid$accuracy <- mapply(function(eng, dv) {
  f <- workflow(rec, divergence_classifier(divergence = dv) |> set_engine(eng)) |>
    fit(train)
  accuracy_vec(test$class, predict(f, test)$.pred_class)
}, grid$engine, grid$divergence)
grid

