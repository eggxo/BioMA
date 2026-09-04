#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: niche_pipeline.R niche_config.tsv")

suppressPackageStartupMessages({
  library(raster)
  library(sf)
})

cfg_raw <- read.delim(args[[1]], header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
cfg <- setNames(cfg_raw$V2, cfg_raw$V1)
get_cfg <- function(key, default = NULL) if (key %in% names(cfg)) cfg[[key]] else default
csv_cfg <- function(key) {
  value <- get_cfg(key, "")
  if (!nzchar(value)) character() else trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}
num_cfg <- function(key) as.numeric(csv_cfg(key))

work_dir <- get_cfg("work_dir")
occurrence_csv <- get_cfg("occurrence_csv")
current_env_dir <- get_cfg("current_env_dir")
future_root <- get_cfg("future_root")
mask_path <- get_cfg("mask_shp")
out_dir <- get_cfg("output_dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rasters"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(occurrence_csv)) stop("Occurrence CSV does not exist: ", occurrence_csv)
if (!dir.exists(current_env_dir)) stop("Current climate directory does not exist: ", current_env_dir)
if (!file.exists(mask_path)) stop("Mask shapefile does not exist: ", mask_path)

maxent_jar <- get_cfg("maxent_jar", "")
if (nzchar(maxent_jar) && !file.exists(maxent_jar)) stop("maxent_jar does not exist: ", maxent_jar)
dismo_package <- system.file(package = "dismo")
if (!nzchar(dismo_package) || !dir.exists(dismo_package)) stop("The dismo R package is not installed")
if (nzchar(maxent_jar)) {
  # dismo only discovers maxent.jar inside its own package. Use a temporary
  # package copy so a licensed user-supplied jar works without modifying the
  # installed R library.
  shadow_library <- tempfile("bioma-dismo-")
  shadow_package <- file.path(shadow_library, "dismo")
  dir.create(shadow_package, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(dismo_package, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  copied <- file.copy(entries, shadow_package, recursive = TRUE)
  if (!length(copied) || any(!copied)) stop("Could not stage the dismo package for the configured MaxEnt jar")
  dir.create(file.path(shadow_package, "java"), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(maxent_jar, file.path(shadow_package, "java", "maxent.jar"), overwrite = TRUE)) {
    stop("Could not stage configured maxent_jar")
  }
  .libPaths(c(shadow_library, .libPaths()))
}
suppressPackageStartupMessages(library(dismo))
installed_jar <- system.file("java/maxent.jar", package = "dismo")
if (!nzchar(installed_jar) || !file.exists(installed_jar)) {
  stop("dismo MaxEnt jar is not available; set maxent_jar in the niche configuration")
}

bio_files <- function(directory) {
  files <- list.files(directory, pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE)
  if (!length(files)) stop("No tif files in: ", directory)
  names_lower <- tolower(basename(files))
  hit <- regexec("bio[_-]?([0-9]+)", names_lower)
  pieces <- regmatches(names_lower, hit)
  ids <- vapply(pieces, function(x) if (length(x) >= 2) as.integer(x[[2]]) else NA_integer_, integer(1))
  if (anyNA(ids)) stop("Could not identify BIO number in: ", paste(basename(files[is.na(ids)]), collapse = ", "))
  if (anyDuplicated(ids)) stop("Duplicate BIO number in: ", directory)
  ord <- order(ids)
  setNames(files[ord], paste0("bio", ids[ord]))
}

mask_sf <- st_read(mask_path, quiet = TRUE)
if (is.na(st_crs(mask_sf))) {
  st_crs(mask_sf) <- 4326
} else if (!identical(st_crs(mask_sf)$epsg, 4326L)) {
  mask_sf <- st_transform(mask_sf, 4326)
}
if (any(!st_is_valid(mask_sf))) mask_sf <- st_make_valid(mask_sf)
mask_sp <- as(mask_sf, "Spatial")

read_env <- function(directory) {
  files <- bio_files(directory)
  env <- raster::stack(unname(files))
  names(env) <- names(files)
  if (is.na(raster::crs(env))) raster::crs(env) <- raster::crs("+proj=longlat +datum=WGS84 +no_defs")
  env <- raster::crop(env, raster::extent(mask_sp))
  env <- raster::mask(env, mask_sp)
  env
}

current_all <- read_env(current_env_dir)
current_ids <- names(current_all)
if (!length(current_ids)) stop("No current BIO layers remain after masking")

occ <- read.csv(occurrence_csv, stringsAsFactors = FALSE)
if (!all(c("lon", "lat") %in% names(occ))) stop("Occurrence CSV must contain lon and lat")
occ$lon <- as.numeric(occ$lon)
occ$lat <- as.numeric(occ$lat)
occ <- occ[is.finite(occ$lon) & is.finite(occ$lat), , drop = FALSE]
occ <- unique(occ)
bb <- st_bbox(mask_sf)
occ <- occ[occ$lon >= bb[["xmin"]] & occ$lon <= bb[["xmax"]] & occ$lat >= bb[["ymin"]] & occ$lat <= bb[["ymax"]], , drop = FALSE]
if (nrow(occ) < 5) stop("At least five occurrence points are required inside the mask")

occ_values <- raster::extract(current_all, occ[, c("lon", "lat")])
occ_ok <- complete.cases(occ_values)
occ_qc <- data.frame(input_rows = nrow(occ), used_rows = sum(occ_ok), dropped_rows = sum(!occ_ok), stringsAsFactors = FALSE)
write.table(occ_qc, file.path(out_dir, "tables", "occurrence_qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
occ <- occ[occ_ok, , drop = FALSE]
occ_values <- occ_values[occ_ok, , drop = FALSE]
if (nrow(occ) < 5) stop("Fewer than five occurrence points have complete BIO values inside the mask")

set.seed(as.integer(get_cfg("seed", "123")))
background_n <- as.integer(get_cfg("background_n", "10000"))
bg <- dismo::randomPoints(current_all[[1]], n = min(background_n, raster::ncell(current_all[[1]])), p = occ[, c("lon", "lat")], excludep = TRUE)
bg <- as.data.frame(bg)
names(bg) <- c("lon", "lat")
bg_values <- raster::extract(current_all, bg[, c("lon", "lat")])
bg <- bg[complete.cases(bg_values), , drop = FALSE]
bg_values <- bg_values[complete.cases(bg_values), , drop = FALSE]
if (nrow(bg) < 20) stop("Too few valid background points after masking")

scope <- tolower(get_cfg("correlation_scope", "occurrence"))
cor_values <- switch(scope,
  occurrence = occ_values,
  background = bg_values,
  mask = {
    n_mask <- min(10000L, raster::ncell(current_all[[1]]))
    pts <- dismo::randomPoints(current_all[[1]], n = n_mask)
    raster::extract(current_all, pts)
  },
  stop("Unsupported correlation_scope: ", scope)
)
cor_values <- as.data.frame(cor_values)
names(cor_values) <- current_ids
cor_values <- cor_values[, vapply(cor_values, function(x) sum(is.finite(x)) >= 3, logical(1)), drop = FALSE]
if (!ncol(cor_values)) stop("No usable variables for correlation filtering")

cor_method <- tolower(get_cfg("correlation_method", "pearson"))
cor_threshold <- as.numeric(get_cfg("correlation_threshold", "0.8"))
cor_matrix <- cor(cor_values, use = "pairwise.complete.obs", method = cor_method)
write.csv(cor_matrix, file.path(out_dir, "tables", "correlation_matrix.csv"), row.names = TRUE)

cor_filter <- function(values, threshold, method) {
  keep <- names(values)
  if (length(keep) < 2) return(keep)
  repeat {
    cm <- abs(cor(values[, keep, drop = FALSE], use = "pairwise.complete.obs", method = method))
    cm[!is.finite(cm)] <- 0
    diag(cm) <- 0
    max_cor <- max(cm)
    if (!is.finite(max_cor) || max_cor <= threshold || length(keep) < 2) break
    pair <- which(cm == max_cor, arr.ind = TRUE)[1, ]
    mean_cor <- colMeans(cm)
    drop_idx <- if (mean_cor[pair[[1]]] >= mean_cor[pair[[2]]]) pair[[1]] else pair[[2]]
    keep <- keep[-drop_idx]
  }
  keep
}

filtered_vars <- cor_filter(cor_values, cor_threshold, cor_method)
if (!length(filtered_vars)) stop("Correlation filtering removed every variable")
writeLines(filtered_vars, file.path(out_dir, "tables", "variables_after_correlation.txt"))

subset_sizes <- unique(as.integer(num_cfg("subset_sizes")))
subset_sizes <- subset_sizes[is.finite(subset_sizes) & subset_sizes > 0]
param_count <- length(csv_cfg("feature_classes")) * length(num_cfg("beta_multipliers"))
max_models <- as.integer(get_cfg("max_models", "500"))
max_sets <- max(1L, floor(max_models / max(1L, param_count)))

make_var_sets <- function(vars, sizes, max_sets, seed) {
  set.seed(seed)
  keys <- character()
  sets <- list()
  add_set <- function(x) {
    x <- sort(unique(x))
    key <- paste(x, collapse = ",")
    if (!(key %in% keys) && length(sets) < max_sets) {
      keys <<- c(keys, key)
      sets[[length(sets) + 1L]] <<- x
    }
  }
  if (length(vars) <= max(sizes)) add_set(vars)
  for (size in sizes) {
    size <- min(size, length(vars))
    if (size < 1 || length(sets) >= max_sets) next
    total <- choose(length(vars), size)
    if (is.finite(total) && total <= 10000) {
      combos <- combn(vars, size)
      order_idx <- sample(seq_len(ncol(combos)))
      for (j in order_idx) {
        add_set(combos[, j])
        if (length(sets) >= max_sets) break
      }
    } else {
      attempts <- 0L
      while (length(sets) < max_sets && attempts < max_sets * 50L) {
        add_set(sample(vars, size))
        attempts <- attempts + 1L
      }
    }
  }
  sets
}

var_sets <- make_var_sets(filtered_vars, subset_sizes, max_sets, as.integer(get_cfg("seed", "123")))
if (!length(var_sets)) var_sets <- list(filtered_vars)
features <- toupper(csv_cfg("feature_classes"))
betas <- as.numeric(num_cfg("beta_multipliers"))
grid <- expand.grid(var_set_id = seq_along(var_sets), feature_class = features, beta_multiplier = betas, stringsAsFactors = FALSE)
if (nrow(grid) > max_models) {
  set.seed(as.integer(get_cfg("seed", "123")))
  grid <- grid[sample(seq_len(nrow(grid)), max_models), , drop = FALSE]
}
grid$model_id <- seq_len(nrow(grid))
grid$n_variables <- vapply(grid$var_set_id, function(i) length(var_sets[[i]]), integer(1))

feature_flags <- function(feature_class) {
  chars <- strsplit(feature_class, "", fixed = TRUE)[[1]]
  c(linear = "L" %in% chars, quadratic = "Q" %in% chars, product = "P" %in% chars, threshold = "T" %in% chars, hinge = "H" %in% chars)
}
maxent_args <- function(feature_class, beta) {
  flags <- feature_flags(feature_class)
  c("responsecurves=false", "jackknife=false", "pictures=false",
    paste0("linear=", tolower(flags[["linear"]])), paste0("quadratic=", tolower(flags[["quadratic"]])),
    paste0("product=", tolower(flags[["product"]])), paste0("threshold=", tolower(flags[["threshold"]])),
    paste0("hinge=", tolower(flags[["hinge"]])), paste0("betamultiplier=", beta),
    "outputformat=cloglog", "extrapolate=false", "doclamp=true", "writeclampgrid=true")
}

tss_score <- function(pred_p, pred_a) {
  pred <- c(pred_p, pred_a)
  obs <- c(rep(1, length(pred_p)), rep(0, length(pred_a)))
  thresholds <- seq(0, 1, length.out = 201)
  scores <- vapply(thresholds, function(thr) {
    hit <- pred >= thr
    tp <- sum(hit & obs == 1); fn <- sum(!hit & obs == 1)
    fp <- sum(hit & obs == 0); tn <- sum(!hit & obs == 0)
    sens <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    spec <- if ((tn + fp) == 0) 0 else tn / (tn + fp)
    sens + spec - 1
  }, numeric(1))
  max(scores, na.rm = TRUE)
}

cv_folds <- min(as.integer(get_cfg("cv_folds", "5")), nrow(occ))
set.seed(as.integer(get_cfg("seed", "123")))
occ_fold <- dismo::kfold(occ[, c("lon", "lat")], k = cv_folds)
bg_fold <- sample(rep(seq_len(cv_folds), length.out = nrow(bg)))
metric_rows <- vector("list", nrow(grid))
selection_metric <- tolower(get_cfg("selection_metric", "auc"))

for (gi in seq_len(nrow(grid))) {
  g <- grid[gi, ]
  vars <- var_sets[[g$var_set_id]]
  fold_rows <- list()
  for (fold in seq_len(cv_folds)) {
    train_occ <- occ[occ_fold != fold, c("lon", "lat"), drop = FALSE]
    test_occ <- occ[occ_fold == fold, c("lon", "lat"), drop = FALSE]
    train_bg <- bg[bg_fold != fold, c("lon", "lat"), drop = FALSE]
    test_bg <- bg[bg_fold == fold, c("lon", "lat"), drop = FALSE]
    model_dir <- file.path(out_dir, "models", sprintf("tuning_%03d", g$model_id), sprintf("fold_%02d", fold))
    dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
    fit <- tryCatch(dismo::maxent(current_all[[vars]], p = train_occ, a = train_bg, args = maxent_args(g$feature_class, g$beta_multiplier), path = model_dir), error = function(e) NULL)
    if (is.null(fit)) next
    pred_r <- tryCatch(raster::predict(fit, current_all[[vars]], progress = ""), error = function(e) NULL)
    if (is.null(pred_r)) next
    pred_p <- raster::extract(pred_r, test_occ)
    pred_a <- raster::extract(pred_r, test_bg)
    pred_p <- pred_p[is.finite(pred_p)]; pred_a <- pred_a[is.finite(pred_a)]
    if (length(pred_p) < 1 || length(pred_a) < 1) next
    ev <- tryCatch(dismo::evaluate(p = pred_p, a = pred_a), error = function(e) NULL)
    fold_rows[[length(fold_rows) + 1L]] <- data.frame(model_id = g$model_id, fold = fold, auc = if (is.null(ev)) NA_real_ else ev@auc, tss = tss_score(pred_p, pred_a), stringsAsFactors = FALSE)
  }
  folds <- if (length(fold_rows)) do.call(rbind, fold_rows) else data.frame(auc = numeric(), tss = numeric())
  metric_rows[[gi]] <- data.frame(model_id = g$model_id, var_set_id = g$var_set_id, feature_class = g$feature_class, beta_multiplier = g$beta_multiplier, n_variables = g$n_variables, successful_folds = nrow(folds), auc_mean = if (nrow(folds)) mean(folds$auc, na.rm = TRUE) else NA_real_, auc_sd = if (nrow(folds) > 1) sd(folds$auc, na.rm = TRUE) else NA_real_, tss_mean = if (nrow(folds)) mean(folds$tss, na.rm = TRUE) else NA_real_, tss_sd = if (nrow(folds) > 1) sd(folds$tss, na.rm = TRUE) else NA_real_, stringsAsFactors = FALSE)
}

metrics <- do.call(rbind, metric_rows)
write.csv(metrics, file.path(out_dir, "tuning_metrics.csv"), row.names = FALSE)
score_col <- if (selection_metric == "tss") "tss_mean" else "auc_mean"
metrics_ok <- metrics[is.finite(metrics[[score_col]]) & metrics$successful_folds > 0, , drop = FALSE]
if (!nrow(metrics_ok)) stop("No MaxEnt tuning model completed successfully")
ord <- order(-metrics_ok[[score_col]], -metrics_ok$auc_mean, -metrics_ok$tss_mean, metrics_ok$n_variables, metrics_ok$beta_multiplier)
best <- metrics_ok[ord[[1]], , drop = FALSE]
best_vars <- var_sets[[best$var_set_id]]
write.table(data.frame(parameter = c("model_id", "feature_class", "beta_multiplier", "variables", "selection_metric", "selection_score", "auc_mean", "tss_mean"), value = c(best$model_id, best$feature_class, best$beta_multiplier, paste(best_vars, collapse = ","), selection_metric, best[[score_col]], best$auc_mean, best$tss_mean)), file.path(out_dir, "selected_model.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(best_vars, file.path(out_dir, "tables", "selected_variables.txt"))

best_dir <- file.path(out_dir, "models", "best_model")
dir.create(best_dir, recursive = TRUE, showWarnings = FALSE)
final_fit <- dismo::maxent(current_all[[best_vars]], p = occ[, c("lon", "lat")], a = bg[, c("lon", "lat")], args = maxent_args(best$feature_class, best$beta_multiplier), path = best_dir)
current_pred <- raster::predict(final_fit, current_all[[best_vars]], progress = "")
writeRaster(current_pred, file.path(out_dir, "rasters", "current_suitability.tif"), overwrite = TRUE)

periods <- csv_cfg("periods")
scenarios <- csv_cfg("scenarios")
gcms <- csv_cfg("gcms")
ensemble_method <- tolower(get_cfg("ensemble_method", "mean"))
future_preds <- list()
manifest_rows <- list()
for (period in periods) for (scenario in scenarios) for (gcm in gcms) {
  key <- paste(period, scenario, gcm, sep = "__")
  directory <- file.path(future_root, paste(period, scenario, gcm, sep = "-"))
  if (!dir.exists(directory)) stop("Future climate directory missing: ", directory)
  future_all <- read_env(directory)
  missing_vars <- setdiff(best_vars, names(future_all))
  if (length(missing_vars)) stop("Future directory missing selected variables: ", directory, " [", paste(missing_vars, collapse = ","), "]")
  future_env <- future_all[[best_vars]]
  current_env <- current_all[[best_vars]]
  same_grid <- raster::compareRaster(current_env, future_env, extent = TRUE, rowcol = TRUE, crs = TRUE, res = TRUE, orig = TRUE, rotation = TRUE, stopiffalse = FALSE)
  if (!same_grid) {
    future_env <- raster::resample(future_env, current_env, method = "bilinear")
    future_env <- raster::mask(future_env, mask_sp)
  }
  pred <- raster::predict(final_fit, future_env, progress = "")
  future_preds[[key]] <- pred
  out_name <- paste0("future_suitability_", period, "_", scenario, "_", gcm, ".tif")
  writeRaster(pred, file.path(out_dir, "rasters", out_name), overwrite = TRUE)
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(period = period, scenario = scenario, gcm = gcm, path = file.path("rasters", out_name), stringsAsFactors = FALSE)
}

for (period in periods) for (scenario in scenarios) {
  keys <- paste(period, scenario, gcms, sep = "__")
  if (!all(keys %in% names(future_preds))) stop("No complete GCM set for ", period, " ", scenario)
  stack_future <- raster::stack(future_preds[keys])
  ensemble <- if (ensemble_method == "median") raster::calc(stack_future, median, na.rm = TRUE) else raster::calc(stack_future, mean, na.rm = TRUE)
  ensemble_name <- paste0("future_suitability_", period, "_", scenario, "_ensemble_", ensemble_method, ".tif")
  writeRaster(ensemble, file.path(out_dir, "rasters", ensemble_name), overwrite = TRUE)
  delta <- current_pred - ensemble
  delta_name <- paste0("maladaptation_", period, "_", scenario, "_ensemble_", ensemble_method, ".tif")
  writeRaster(delta, file.path(out_dir, "rasters", delta_name), overwrite = TRUE)
  vals_current <- raster::extract(current_pred, occ[, c("lon", "lat")])
  vals_future <- raster::extract(ensemble, occ[, c("lon", "lat")])
  point_df <- data.frame(lon = occ$lon, lat = occ$lat, period = period, scenario = scenario, current_suitability = vals_current, future_suitability = vals_future, maladaptation = vals_current - vals_future)
  write.csv(point_df, file.path(out_dir, "tables", paste0("maladaptation_points_", period, "_", scenario, ".csv")), row.names = FALSE)
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(period = period, scenario = scenario, gcm = paste0("ensemble_", ensemble_method), path = file.path("rasters", ensemble_name), stringsAsFactors = FALSE)
  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(period = period, scenario = scenario, gcm = paste0("maladaptation_ensemble_", ensemble_method), path = file.path("rasters", delta_name), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, manifest_rows), file.path(out_dir, "tables", "projection_manifest.csv"), row.names = FALSE)
cat("MaxEnt niche pipeline complete\n")
cat("Selected variables:", paste(best_vars, collapse = ", "), "\n")
cat("Feature class:", best$feature_class, " beta:", best$beta_multiplier, "\n")
