#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rsample)
  library(ranger)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript predict_future_one_target.R <target> <split_id> <mtry> <min_n> <trees> <train_csv> <future_csv> <out_csv> [seed_base] [n_workers]\n\n",
      "Example:\n",
      "  Rscript predict_future_one_target.R mean_loadD_relax 1534 1 1 979 pop_geo_niche_predictors_from_TSS.csv future_env_ssp245_2061_2080_mean.csv out.csv 1234 8\n"
    )
  )
}

target     <- args[[1]]
split_id   <- as.integer(args[[2]])
mtry       <- as.integer(args[[3]])
min_n      <- as.integer(args[[4]])
trees      <- as.integer(args[[5]])
train_csv  <- args[[6]]
future_csv <- args[[7]]
out_csv    <- args[[8]]
seed_base  <- ifelse(length(args) >= 9,  as.integer(args[[9]]), 1234)
n_workers  <- ifelse(length(args) >= 10, as.integer(args[[10]]), 8)

predictors <- c(
  "bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9",
  "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"
)

safe_rsq <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  truth <- truth[ok]
  pred  <- pred[ok]
  if (length(truth) < 3) return(NA_real_)
  if (sd(truth) == 0 || sd(pred) == 0) return(NA_real_)
  cor(truth, pred)^2
}

safe_rmse <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  truth <- truth[ok]
  pred  <- pred[ok]
  if (length(truth) == 0) return(NA_real_)
  sqrt(mean((truth - pred)^2))
}

safe_mae <- function(truth, pred) {
  ok <- is.finite(truth) & is.finite(pred)
  truth <- truth[ok]
  pred  <- pred[ok]
  if (length(truth) == 0) return(NA_real_)
  mean(abs(truth - pred))
}

project_to_svd_space <- function(df, predictors, svd_obj = NULL, center = NULL, scale_ = NULL) {
  X <- as.matrix(df[, predictors, drop = FALSE])

  if (is.null(svd_obj)) {
    Xs <- scale(X)
    center <- attr(Xs, "scaled:center")
    scale_ <- attr(Xs, "scaled:scale")
    svd_obj <- svd(Xs)
  } else {
    Xs <- scale(X, center = center, scale = scale_)
  }

  d_inv <- ifelse(svd_obj$d > 1e-8, 1 / svd_obj$d, 0)

  Z <- Xs %*% svd_obj$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(svd_obj$v)
  Z <- as.data.frame(Z)
  colnames(Z) <- paste0("SV", seq_len(ncol(Z)))

  list(
    data   = Z,
    svd    = svd_obj,
    center = center,
    scale  = scale_
  )
}

cat("Reading training data:", train_csv, "\n")
train_df <- read.csv(train_csv, header = TRUE, check.names = FALSE)

required_cols <- c("pop", predictors, target)
missing_cols <- setdiff(required_cols, colnames(train_df))
if (length(missing_cols) > 0) {
  stop("Missing columns in training csv: ", paste(missing_cols, collapse = ", "))
}

train_df[predictors] <- lapply(train_df[predictors], as.numeric)
train_df[[target]]   <- as.numeric(train_df[[target]])

model_df <- train_df %>%
  dplyr::select(pop, all_of(predictors), all_of(target)) %>%
  tidyr::drop_na()

if (nrow(model_df) < 10) {
  stop("Too few complete rows in training data for target: ", target)
}

cat("Rebuilding split:", split_id, "\n")
set.seed(seed_base + split_id)
split_obj  <- initial_split(model_df, prop = 0.8)
train_data <- training(split_obj)
test_data  <- testing(split_obj)

cat("Train n =", nrow(train_data), "; Test n =", nrow(test_data), "\n")

# SVD fitted only on training data
svd_train <- project_to_svd_space(
  df = train_data,
  predictors = predictors
)

train_z <- svd_train$data
train_z$y <- train_data[[target]]

test_z <- project_to_svd_space(
  df = test_data,
  predictors = predictors,
  svd_obj = svd_train$svd,
  center = svd_train$center,
  scale_ = svd_train$scale
)$data

# fit final model with fixed params
cat("Fitting ranger model...\n")
fit <- ranger(
  formula = y ~ .,
  data = train_z,
  mtry = mtry,
  min.node.size = min_n,
  num.trees = trees,
  importance = "permutation",
  num.threads = n_workers,
  seed = seed_base + split_id
)

# sanity check on test split
test_pred <- predict(fit, data = test_z)$predictions
test_truth <- test_data[[target]]

test_rsq  <- safe_rsq(test_truth, test_pred)
test_rmse <- safe_rmse(test_truth, test_pred)
test_mae  <- safe_mae(test_truth, test_pred)

cat("Sanity check on recreated test split:\n")
cat("  target    =", target, "\n")
cat("  split_id  =", split_id, "\n")
cat("  test_rsq  =", round(test_rsq, 6), "\n")
cat("  test_rmse =", round(test_rmse, 6), "\n")
cat("  test_mae  =", round(test_mae, 6), "\n")

# current prediction for all populations
all_z <- project_to_svd_space(
  df = model_df,
  predictors = predictors,
  svd_obj = svd_train$svd,
  center = svd_train$center,
  scale_ = svd_train$scale
)$data

current_pred <- predict(fit, data = all_z)$predictions

current_out <- data.frame(
  pop = model_df$pop,
  observed = model_df[[target]],
  predicted = current_pred,
  target = target,
  split_id = split_id,
  mtry = mtry,
  min_n = min_n,
  trees = trees,
  stringsAsFactors = FALSE
)

current_check_csv <- sub("\\.csv$", ".current_check.csv", out_csv)
if (identical(current_check_csv, out_csv)) {
  current_check_csv <- paste0(out_csv, ".current_check.csv")
}
write.csv(current_out, current_check_csv, row.names = FALSE)

# future prediction
cat("Reading future data:", future_csv, "\n")
future_df <- read.csv(future_csv, header = TRUE, check.names = FALSE)

missing_future_cols <- setdiff(predictors, colnames(future_df))
if (length(missing_future_cols) > 0) {
  stop("Missing predictor columns in future csv: ", paste(missing_future_cols, collapse = ", "))
}

future_df[predictors] <- lapply(future_df[predictors], as.numeric)

if (!("pop" %in% colnames(future_df))) {
  future_df$pop <- paste0("future_row_", seq_len(nrow(future_df)))
}

future_z <- project_to_svd_space(
  df = future_df,
  predictors = predictors,
  svd_obj = svd_train$svd,
  center = svd_train$center,
  scale_ = svd_train$scale
)$data

future_pred <- predict(fit, data = future_z)$predictions

out_df <- future_df
out_df$predicted_target <- target
out_df$predicted_value  <- future_pred
out_df$split_id         <- split_id
out_df$mtry             <- mtry
out_df$min_n            <- min_n
out_df$trees            <- trees
out_df$recreated_test_rsq  <- test_rsq
out_df$recreated_test_rmse <- test_rmse
out_df$recreated_test_mae  <- test_mae

write.csv(out_df, out_csv, row.names = FALSE)

cat("Saved future predictions to:", out_csv, "\n")
cat("Saved current check to:", current_check_csv, "\n")
cat("Done.\n")
