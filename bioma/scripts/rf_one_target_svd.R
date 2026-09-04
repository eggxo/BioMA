#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidymodels)
  library(ranger)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(tibble)
  library(future)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript rf_one_target_svd.R <target> <workdir> [n_splits] [cv_v] [cv_repeats] [n_workers] [grid_size] [seed]\n\n",
      "Example:\n",
      "  Rscript rf_one_target_svd.R mean_loadM_strict /path/to/load_results 1000 5 10 16\n"
    )
  )
}

target      <- args[[1]]
workdir     <- args[[2]]
n_splits    <- ifelse(length(args) >= 3, as.integer(args[[3]]), 1000)
cv_v        <- ifelse(length(args) >= 4, as.integer(args[[4]]), 5)
cv_repeats  <- ifelse(length(args) >= 5, as.integer(args[[5]]), 10)
n_workers   <- ifelse(length(args) >= 6, as.integer(args[[6]]), max(1, min(32, future::availableCores() - 1)))
grid_size   <- ifelse(length(args) >= 7, as.integer(args[[7]]), 50L)
seed_base   <- ifelse(length(args) >= 8, as.integer(args[[8]]), 1234L)
if (!is.finite(grid_size) || grid_size < 1) stop("grid_size must be a positive integer")
if (!is.finite(seed_base)) stop("seed must be an integer")

setwd(workdir)

# =====================
# 1. 读入数据
# =====================
data <- read.csv(
  "pop_geo_niche_predictors_from_TSS.csv",
  header = TRUE,
  check.names = FALSE
)

predictors <- c(
  "bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio9",
  "bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19"
)

allowed_targets <- c(
  "mean_loadM_strict",
  "mean_loadD_strict",
  "mean_loadM_relax",
  "mean_loadD_relax",
  "Deleterious_synonymous_homozygous_variant_ratio",
  "Deleterious_synonymous_heterozygous_variant_ratio",
  "LoF_synonymous_homozygous_variant_ratio",
  "LoF_synonymous_heterozygous_variant_ratio",
  "Tolerance_synonymous_homozygous_variant_ratio",
  "Tolerance_synonymous_heterozygous_variant_ratio"
)

if (!(target %in% allowed_targets)) {
  stop(paste0(
    "Unknown target: ", target, "\nAllowed targets:\n  ",
    paste(allowed_targets, collapse = "\n  ")
  ))
}

required_cols <- c("pop", predictors, target)
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

data[predictors] <- lapply(data[predictors], as.numeric)
data[[target]]   <- as.numeric(data[[target]])

# =====================
# 2. 颜色映射
# =====================
target_colors <- c(
  "mean_loadM_strict" = "#49C0C9",
  "mean_loadD_strict" = "#6C9DF2",
  "mean_loadM_relax"  = "#7E57C2",
  "mean_loadD_relax"  = "#E07AD9",
  "Deleterious_synonymous_homozygous_variant_ratio"   = "#B89F1D",
  "Deleterious_synonymous_heterozygous_variant_ratio" = "#E98473",
  "LoF_synonymous_homozygous_variant_ratio"           = "#43B84F",
  "LoF_synonymous_heterozygous_variant_ratio"         = "#2E8B57",
  "Tolerance_synonymous_homozygous_variant_ratio"     = "#A68A7A",
  "Tolerance_synonymous_heterozygous_variant_ratio"   = "#BD7BB1"
)

this_color <- target_colors[[target]]
if (is.null(this_color)) this_color <- "#4C78A8"

# =====================
# 3. 输出目录
# =====================
safe_target <- gsub("[^A-Za-z0-9]+", "_", target)
outdir <- file.path(workdir, paste0("RF_", safe_target))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots_tuning_each_split"),       showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots_test_eval_each_split"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots_current_eval_each_split"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "plots_summary"),                 showWarnings = FALSE, recursive = TRUE)

# =====================
# 4. 并行
# =====================
plan(multisession, workers = n_workers)

# =====================
# 5. 工具函数
# =====================
save_plot_both <- function(plot_obj, file_stub, width, height, dpi = 600) {
  ggsave(
    paste0(file_stub, ".png"),
    plot_obj,
    width = width, height = height,
    dpi = dpi, bg = "white"
  )
  ggsave(
    paste0(file_stub, ".pdf"),
    plot_obj,
    width = width, height = height,
    bg = "white"
  )
}

project_to_svd_space <- function(df, predictors, target = NULL,
                                 svd_obj = NULL, center = NULL, scale_ = NULL) {

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

  if (!is.null(target)) {
    Z[[target]] <- df[[target]]
  }

  list(
    data   = Z,
    svd    = svd_obj,
    center = center,
    scale  = scale_
  )
}

safe_lm_stats <- function(df) {
  df <- df %>% dplyr::filter(is.finite(truth), is.finite(predicted))

  if (nrow(df) < 3) {
    return(data.frame(r2 = NA_real_, p = NA_real_, x = NA_real_, y = NA_real_))
  }

  if (sd(df$truth) == 0 || sd(df$predicted) == 0) {
    return(data.frame(r2 = NA_real_, p = NA_real_, x = NA_real_, y = NA_real_))
  }

  fit <- tryCatch(lm(predicted ~ truth, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.frame(r2 = NA_real_, p = NA_real_, x = NA_real_, y = NA_real_))
  }

  p_val <- tryCatch(summary(fit)$coefficients[2, 4], error = function(e) NA_real_)

  x_rng <- range(df$truth, na.rm = TRUE)
  y_rng <- range(df$predicted, na.rm = TRUE)

  x_span <- diff(x_rng)
  y_span <- diff(y_rng)

  if (x_span == 0) x_span <- 1
  if (y_span == 0) y_span <- 1

  data.frame(
    r2 = summary(fit)$r.squared,
    p  = p_val,
    x  = x_rng[1] + 0.02 * x_span,
    y  = y_rng[2] - 0.08 * y_span
  )
}

make_stat_label <- function(stat_df, split_id, n_points) {
  r2_txt <- ifelse(is.na(stat_df$r2), "NA", formatC(stat_df$r2, format = "f", digits = 3))
  p_txt  <- ifelse(is.na(stat_df$p),  "NA", formatC(stat_df$p,  format = "e", digits = 2))

  paste0(
    "split = ", split_id,
    "\nR² = ", r2_txt,
    "\np = ", p_txt,
    "\nn = ", n_points
  )
}

safe_rsq_vec <- function(truth, estimate) {
  truth <- as.numeric(truth)
  estimate <- as.numeric(estimate)

  if (length(unique(truth[is.finite(truth)])) < 2) return(NA_real_)
  if (length(unique(estimate[is.finite(estimate)])) < 2) return(NA_real_)

  suppressWarnings(yardstick::rsq_vec(truth = truth, estimate = estimate))
}

plot_tuning_facet <- function(rf_results_rmse, best_params_plot, target_name, this_color, n_folds) {

  plot_data_long <- rf_results_rmse %>%
    dplyr::select(mean, mtry, min_n, trees) %>%
    tidyr::pivot_longer(
      cols = c(trees, mtry, min_n),
      names_to = "param",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      param = factor(param, levels = c("mtry", "min_n", "trees")),
      value = as.numeric(value)
    )

  best_params_long <- best_params_plot %>%
    dplyr::select(mean, mtry, min_n, trees) %>%
    tidyr::pivot_longer(
      cols = c(trees, mtry, min_n),
      names_to = "param",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      param = factor(param, levels = c("mtry", "min_n", "trees")),
      value = as.numeric(value)
    )

  ggplot(plot_data_long, aes(x = value, y = mean)) +
    geom_point(color = this_color, alpha = 0.95, size = 3) +
    geom_smooth(
      method = "loess", se = TRUE,
      color = this_color, fill = this_color,
      linewidth = 1.1, alpha = 0.15, span = 0.9
    ) +
    geom_point(
      data = best_params_long,
      aes(x = value, y = mean),
      inherit.aes = FALSE,
      shape = 18, size = 3.8, color = "red"
    ) +
    facet_wrap(~ param, scales = "free_x", nrow = 1) +
    labs(
      title = paste0(target_name, " | tuning"),
      x = "Parameter value",
      y = paste0("Mean RMSE across ", n_folds, "-fold CV")
    ) +
    theme_classic() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      axis.ticks.length = unit(0.14, "cm"),
      panel.spacing.x = unit(0.8, "lines"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin = margin(6, 8, 6, 6)
    )
}

plot_eval_single <- function(plot_df, stat_df, target_name, this_color, y_lab = NULL) {

  if (is.null(y_lab)) y_lab <- paste0("Predicted ", target_name)

  ggplot(plot_df, aes(x = truth, y = predicted)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55", linewidth = 0.8) +
    geom_point(color = this_color, size = 3, alpha = 0.95) +
    geom_smooth(
      method = "lm", se = TRUE,
      color = this_color, linewidth = 1.1,
      alpha = 0.15, fill = this_color
    ) +
    geom_text(
      data = stat_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 3.5, color = "black"
    ) +
    labs(
      title = target_name,
      x = paste0("Observed ", target_name),
      y = y_lab
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

run_one_split <- function(model_df, predictors, target, split_id, seed_base = 1234) {

  set.seed(seed_base + split_id)
  split_obj  <- initial_split(model_df, prop = 0.8)
  train_data <- training(split_obj)
  test_data  <- testing(split_obj)

  svd_train <- project_to_svd_space(
    df         = train_data,
    predictors = predictors,
    target     = target
  )
  train_z <- svd_train$data

  svd_test <- project_to_svd_space(
    df         = test_data,
    predictors = predictors,
    target     = target,
    svd_obj    = svd_train$svd,
    center     = svd_train$center,
    scale_     = svd_train$scale
  )
  test_z <- svd_test$data

  svd_all <- project_to_svd_space(
    df         = model_df,
    predictors = predictors,
    target     = NULL,
    svd_obj    = svd_train$svd,
    center     = svd_train$center,
    scale_     = svd_train$scale
  )
  all_z <- svd_all$data

  rf_mod <- rand_forest(
    mtry  = tune(),
    min_n = tune(),
    trees = tune()
  ) %>%
    set_engine("ranger", importance = "permutation", num.threads = 1) %>%
    set_mode("regression")

  rf_wf <- workflow() %>%
    add_formula(as.formula(paste(target, "~ ."))) %>%
    add_model(rf_mod)

  max_mtry <- ncol(train_z) - 1
  if (max_mtry < 1) stop("No predictors remain in SVD space.")

  max_min_n <- max(2L, min(10L, nrow(train_z) - 1L))

  set.seed(seed_base + 10000 + split_id)
  rf_grid <- grid_latin_hypercube(
    mtry(range = c(1L, min(6L, max_mtry))),
    min_n(range = c(1L, max_min_n)),
    trees(range = c(600L, 1200L)),
    size = grid_size
  )

  set.seed(seed_base + 20000 + split_id)
  cv_splits <- vfold_cv(train_z, v = cv_v, repeats = cv_repeats)

  ctrl <- control_grid(
    verbose = FALSE,
    save_pred = FALSE,
    save_workflow = FALSE,
    parallel_over = "resamples",
    allow_par = TRUE
  )

  rf_tune <- tune_grid(
    rf_wf,
    resamples = cv_splits,
    grid = rf_grid,
    metrics = metric_set(rmse, rsq, mae),
    control = ctrl
  )

  best_rf <- select_best(rf_tune, metric = "rmse")

  cv_rsq_best <- collect_metrics(rf_tune) %>%
    dplyr::filter(.metric == "rsq") %>%
    dplyr::inner_join(best_rf, by = c("mtry", "min_n", "trees")) %>%
    dplyr::slice(1)

  final_rf  <- finalize_workflow(rf_wf, best_rf)
  final_fit <- fit(final_rf, data = train_z)

  fit_obj <- extract_fit_parsnip(final_fit)$fit

  oob_mse  <- fit_obj$prediction.error
  oob_rmse <- if (!is.null(oob_mse) && is.finite(oob_mse)) sqrt(oob_mse) else NA_real_
  oob_rsq  <- fit_obj$r.squared
  if (is.null(oob_rsq)) oob_rsq <- NA_real_

  test_pred <- predict(final_fit, new_data = test_z) %>%
    bind_cols(test_data %>% dplyr::select(pop, all_of(target))) %>%
    dplyr::rename(obs = all_of(target)) %>%
    dplyr::rename(predicted = .pred) %>%
    dplyr::mutate(split_id = split_id, target = target) %>%
    dplyr::select(pop, split_id, target, obs, predicted)

  current_pred <- predict(final_fit, new_data = all_z) %>%
    bind_cols(model_df %>% dplyr::select(pop, all_of(target))) %>%
    dplyr::rename(obs = all_of(target)) %>%
    dplyr::rename(predicted = .pred) %>%
    dplyr::mutate(split_id = split_id, target = target) %>%
    dplyr::select(pop, split_id, target, obs, predicted)

  out_metrics <- tibble(
    split_id = split_id,
    target   = target,
    n_train  = nrow(train_z),
    n_test   = nrow(test_z),
    mtry     = best_rf$mtry,
    min_n    = best_rf$min_n,
    trees    = best_rf$trees,
    oob_mse  = as.numeric(oob_mse),
    oob_rmse = as.numeric(oob_rmse),
    oob_rsq  = as.numeric(oob_rsq),
    oob_var_explained_pct = as.numeric(oob_rsq) * 100,
    test_rsq   = safe_rsq_vec(truth = test_pred$obs, estimate = test_pred$predicted),
    test_rmse  = yardstick::rmse_vec(truth = test_pred$obs, estimate = test_pred$predicted),
    test_mae   = yardstick::mae_vec(truth = test_pred$obs, estimate = test_pred$predicted),
    cv_rsq_mean = cv_rsq_best$mean,
    cv_rsq_se   = cv_rsq_best$std_err
  )

  out_imp <- NULL
  if (!is.null(fit_obj$variable.importance)) {
    imp_sv   <- fit_obj$variable.importance
    sv_names <- setdiff(colnames(train_z), target)
    imp_sv   <- imp_sv[sv_names]
    imp_sv[is.na(imp_sv)] <- 0

    lambda <- svd_train$svd$v %*%
      diag(svd_train$svd$d, nrow = length(svd_train$svd$d)) %*%
      t(svd_train$svd$v)

    imp_orig <- as.numeric((lambda^2) %*% matrix(imp_sv, ncol = 1))
    denom <- sum(abs(imp_orig))
    if (denom == 0) denom <- 1

    out_imp <- tibble(
      split_id = split_id,
      target   = target,
      variable = predictors,
      importance_raw = imp_orig,
      importance_abs_rel = abs(imp_orig) / denom
    ) %>%
      arrange(desc(importance_abs_rel))
  }

  list(
    metrics       = out_metrics,
    preds_test    = test_pred,
    preds_current = current_pred,
    imp           = out_imp,
    final_fit     = final_fit,
    best_rf       = best_rf,
    rf_tune       = rf_tune,
    train_pop     = train_data$pop,
    test_pop      = test_data$pop
  )
}

# =====================
# 6. 建模
# =====================
cat("Running target:", target, "\n")
cat("Output dir:", outdir, "\n")
cat("n_splits:", n_splits, " cv_v:", cv_v, " cv_repeats:", cv_repeats,
    " workers:", n_workers, " grid_size:", grid_size, " seed:", seed_base, "\n")

model_df <- data %>%
  dplyr::select(pop, all_of(predictors), all_of(target)) %>%
  tidyr::drop_na()

if (nrow(model_df) < 10) {
  stop(paste("Too few rows for target:", target))
}

all_metrics_list        <- list()
all_preds_test_list     <- list()
all_preds_current_list  <- list()
all_imp_list            <- list()
all_model_list          <- list()
all_splitinfo_list      <- list()

for (split_id in seq_len(n_splits)) {
  cat("  Split:", split_id, "\n")

  res_i <- run_one_split(
    model_df   = model_df,
    predictors = predictors,
    target     = target,
    split_id   = split_id,
    seed_base  = seed_base
  )

  obj_name <- paste(target, split_id, sep = "_")

  all_metrics_list[[obj_name]]       <- res_i$metrics
  all_preds_test_list[[obj_name]]    <- res_i$preds_test
  all_preds_current_list[[obj_name]] <- res_i$preds_current
  all_imp_list[[obj_name]]           <- res_i$imp

  all_model_list[[obj_name]] <- list(
    target    = target,
    split_id  = split_id,
    best_rf   = res_i$best_rf,
    final_fit = res_i$final_fit
  )

  all_splitinfo_list[[obj_name]] <- tibble(
    target    = target,
    split_id  = split_id,
    train_pop = paste(res_i$train_pop, collapse = ";"),
    test_pop  = paste(res_i$test_pop, collapse = ";")
  )

  rf_results <- collect_metrics(res_i$rf_tune)

  best_params_plot <- rf_results %>%
    dplyr::filter(.metric == "rmse") %>%
    dplyr::semi_join(res_i$best_rf, by = c("mtry", "min_n", "trees")) %>%
    dplyr::slice(1)

  plot_data_tune <- rf_results %>%
    dplyr::filter(.metric == "rmse")

  figure_tune <- plot_tuning_facet(
    rf_results_rmse = plot_data_tune,
    best_params_plot = best_params_plot,
    target_name = target,
    this_color = this_color,
    n_folds = cv_v
  )

  save_plot_both(
    figure_tune,
    file.path(outdir, "plots_tuning_each_split",
              paste0("Figure_Tuning_", safe_target, "_split", sprintf("%03d", split_id), "_facet_SVD")),
    width = 12, height = 4
  )

  plot_df_test_i <- res_i$preds_test %>%
    dplyr::transmute(truth = obs, predicted = predicted)

  stat_test_i <- safe_lm_stats(plot_df_test_i)
  stat_test_i$label <- make_stat_label(stat_test_i, split_id = split_id, n_points = nrow(plot_df_test_i))

  figure_test_i <- plot_eval_single(
    plot_df = plot_df_test_i,
    stat_df = stat_test_i,
    target_name = target,
    this_color = this_color,
    y_lab = paste0("Predicted ", target)
  )

  save_plot_both(
    figure_test_i,
    file.path(outdir, "plots_test_eval_each_split",
              paste0("Figure_TestEval_", safe_target, "_split", sprintf("%03d", split_id), "_SVD")),
    width = 6, height = 5
  )

  plot_df_current_i <- res_i$preds_current %>%
    dplyr::transmute(truth = obs, predicted = predicted)

  stat_current_i <- safe_lm_stats(plot_df_current_i)
  stat_current_i$label <- make_stat_label(stat_current_i, split_id = split_id, n_points = nrow(plot_df_current_i))

  figure_current_i <- plot_eval_single(
    plot_df = plot_df_current_i,
    stat_df = stat_current_i,
    target_name = target,
    this_color = this_color,
    y_lab = paste0("Predicted ", target)
  )

  save_plot_both(
    figure_current_i,
    file.path(outdir, "plots_current_eval_each_split",
              paste0("Figure_CurrentTruth_vs_CurrentPrediction_", safe_target, "_split", sprintf("%03d", split_id), "_SVD")),
    width = 6, height = 5
  )
}

# =====================
# 7. 汇总
# =====================
all_metrics       <- bind_rows(all_metrics_list)
all_preds_test    <- bind_rows(all_preds_test_list)
all_preds_current <- bind_rows(all_preds_current_list)
all_imp           <- bind_rows(all_imp_list)
split_info_df     <- bind_rows(all_splitinfo_list)

summary_metrics <- all_metrics %>%
  summarise(
    target          = target[1],
    n_splits        = n(),
    mean_n_train    = mean(n_train, na.rm = TRUE),
    mean_n_test     = mean(n_test, na.rm = TRUE),
    mean_test_rsq   = mean(test_rsq, na.rm = TRUE),
    sd_test_rsq     = sd(test_rsq, na.rm = TRUE),
    mean_test_rmse  = mean(test_rmse, na.rm = TRUE),
    sd_test_rmse    = sd(test_rmse, na.rm = TRUE),
    mean_test_mae   = mean(test_mae, na.rm = TRUE),
    sd_test_mae     = sd(test_mae, na.rm = TRUE),
    mean_cv_rsq     = mean(cv_rsq_mean, na.rm = TRUE),
    sd_cv_rsq       = sd(cv_rsq_mean, na.rm = TRUE),
    mean_oob_rsq    = mean(oob_rsq, na.rm = TRUE),
    sd_oob_rsq      = sd(oob_rsq, na.rm = TRUE),
    best_test_rsq   = max(test_rsq, na.rm = TRUE)
  )

report_eval_by_split <- all_metrics %>%
  dplyr::select(
    split_id, target, n_train, n_test,
    test_rsq, test_rmse, test_mae,
    cv_rsq_mean, cv_rsq_se,
    mtry, min_n, trees,
    oob_rsq, oob_rmse
  ) %>%
  arrange(desc(test_rsq))

report_eval_by_split_full <- report_eval_by_split %>%
  left_join(split_info_df, by = c("target", "split_id"))

imp_summary <- all_imp %>%
  group_by(target, variable) %>%
  summarise(
    importance_abs_rel_mean = mean(importance_abs_rel, na.rm = TRUE),
    importance_abs_rel_sd   = sd(importance_abs_rel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(importance_abs_rel_mean))

write.csv(all_metrics,               file.path(outdir, paste0("RF_metrics_by_split_", safe_target, ".csv")), row.names = FALSE)
write.csv(summary_metrics,           file.path(outdir, paste0("RF_metrics_summary_", safe_target, ".csv")), row.names = FALSE)
write.csv(all_preds_test,            file.path(outdir, paste0("RF_test_predictions_by_split_", safe_target, ".csv")), row.names = FALSE)
write.csv(all_preds_current,         file.path(outdir, paste0("RF_current_predictions_by_split_", safe_target, ".csv")), row.names = FALSE)
write.csv(all_imp,                   file.path(outdir, paste0("RF_importance_by_split_", safe_target, ".csv")), row.names = FALSE)
write.csv(imp_summary,               file.path(outdir, paste0("RF_importance_summary_", safe_target, ".csv")), row.names = FALSE)
write.csv(report_eval_by_split_full, file.path(outdir, paste0("RF_eval_by_split_full_", safe_target, ".csv")), row.names = FALSE)

# 汇总图1：test R² 分布
p_box <- ggplot(all_metrics, aes(x = target, y = test_rsq, fill = target)) +
  geom_boxplot(alpha = 0.8, width = 0.55, outlier.shape = 16) +
  geom_jitter(width = 0.08, alpha = 0.65, size = 2) +
  scale_fill_manual(values = setNames(this_color, target)) +
  theme_bw() +
  labs(
    title = paste0("Repeated split test R²: ", target),
    x = NULL,
    y = "Withheld test R²"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

save_plot_both(
  p_box,
  file.path(outdir, "plots_summary", paste0("Figure_Test_R2_RepeatedSplits_", safe_target)),
  width = 7, height = 5
)

# 最佳 split：test 图
best_split_for_plot <- report_eval_by_split %>%
  slice_max(order_by = test_rsq, n = 1, with_ties = FALSE)

plot_df_best_test <- all_preds_test %>%
  dplyr::inner_join(
    best_split_for_plot %>% dplyr::select(target, split_id),
    by = c("target", "split_id")
  ) %>%
  dplyr::transmute(truth = obs, predicted = predicted)

stat_df_best_test <- safe_lm_stats(plot_df_best_test)
stat_df_best_test$label <- make_stat_label(
  stat_df_best_test,
  split_id = best_split_for_plot$split_id[1],
  n_points = nrow(plot_df_best_test)
)

figure_eval_test_best <- plot_eval_single(
  plot_df = plot_df_best_test,
  stat_df = stat_df_best_test,
  target_name = target,
  this_color = this_color,
  y_lab = paste0("Predicted ", target)
)

save_plot_both(
  figure_eval_test_best,
  file.path(outdir, "plots_summary", paste0("Figure_BestSplit_TestTruth_vs_TestPrediction_", safe_target, "_SVD")),
  width = 6, height = 5
)

# 最佳 split：current 图
plot_df_best_current <- all_preds_current %>%
  dplyr::inner_join(
    best_split_for_plot %>% dplyr::select(target, split_id),
    by = c("target", "split_id")
  ) %>%
  dplyr::transmute(truth = obs, predicted = predicted)

stat_df_best_current <- safe_lm_stats(plot_df_best_current)
stat_df_best_current$label <- make_stat_label(
  stat_df_best_current,
  split_id = best_split_for_plot$split_id[1],
  n_points = nrow(plot_df_best_current)
)

figure_eval_current_best <- plot_eval_single(
  plot_df = plot_df_best_current,
  stat_df = stat_df_best_current,
  target_name = target,
  this_color = this_color,
  y_lab = paste0("Predicted ", target)
)

save_plot_both(
  figure_eval_current_best,
  file.path(outdir, "plots_summary", paste0("Figure_BestSplit_CurrentTruth_vs_CurrentPrediction_", safe_target, "_SVD")),
  width = 6, height = 5
)

saveRDS(all_model_list,     file.path(outdir, paste0("RF_model_objects_", safe_target, ".rds")))
saveRDS(all_splitinfo_list, file.path(outdir, paste0("RF_split_info_list_", safe_target, ".rds")))

cat("Done:", target, "\n")
cat("Output dir:", outdir, "\n")
