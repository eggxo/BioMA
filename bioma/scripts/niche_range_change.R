#!/usr/bin/env Rscript

# Derive habitat transitions from the continuous MaxEnt rasters produced by
# niche_pipeline.R. Values are 1 = Contract, 2 = Invariant, and 3 = Expand;
# zero means unsuitable in both periods. The transition rasters are an
# optional downstream product and never replace the continuous suitability or
# maladaptation rasters.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: niche_range_change.R niche_config.tsv")

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

out_dir <- get_cfg("output_dir")
raster_dir <- file.path(out_dir, "rasters")
range_dir <- file.path(out_dir, "range_change")
range_raster_dir <- file.path(range_dir, "rasters")
table_dir <- file.path(range_dir, "tables")
dir.create(range_raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

periods <- csv_cfg("periods")
scenarios <- csv_cfg("scenarios")
gcms <- csv_cfg("gcms")
if (!length(periods) || !length(scenarios) || !length(gcms)) stop("The range-change step requires periods, scenarios, and gcms")

threshold_raw <- tolower(trimws(get_cfg("range_threshold", "auto")))
if (threshold_raw == "auto") {
  threshold_method <- "maximum_tss"
} else {
  threshold_value <- suppressWarnings(as.numeric(threshold_raw))
  if (!is.finite(threshold_value) || threshold_value <= 0 || threshold_value >= 1) {
    stop("range_threshold must be auto or a number between 0 and 1")
  }
  threshold_method <- "user_fixed"
}

consensus <- tolower(trimws(get_cfg("range_consensus", "majority")))
if (!(consensus %in% c("majority", "all", "any", "minimum"))) stop("Unsupported range_consensus: ", consensus)
minimum_raw <- tolower(trimws(get_cfg("range_minimum_gcm_agreement", "auto")))

current_path <- file.path(raster_dir, "current_suitability.tif")
if (!file.exists(current_path)) stop("Current suitability raster is missing: ", current_path)
current <- raster(current_path)

mask_path <- get_cfg("mask_shp", "")
if (nzchar(mask_path) && file.exists(mask_path)) {
  mask_sf <- st_read(mask_path, quiet = TRUE)
  if (is.na(st_crs(mask_sf))) st_crs(mask_sf) <- 4326
  if (!identical(st_crs(mask_sf)$epsg, 4326L)) mask_sf <- st_transform(mask_sf, 4326)
  if (any(!st_is_valid(mask_sf))) mask_sf <- st_make_valid(mask_sf)
  current <- mask(current, as(mask_sf, "Spatial"))
}

occurrence_csv <- get_cfg("occurrence_csv", "")
if (!nzchar(occurrence_csv) || !file.exists(occurrence_csv)) stop("Occurrence CSV is required for automatic threshold selection")
occ <- read.csv(occurrence_csv, stringsAsFactors = FALSE)
if (!all(c("lon", "lat") %in% names(occ))) stop("Occurrence CSV must contain lon and lat")
occ <- occ[is.finite(as.numeric(occ$lon)) & is.finite(as.numeric(occ$lat)), , drop = FALSE]
occ$lon <- as.numeric(occ$lon)
occ$lat <- as.numeric(occ$lat)

background_path <- file.path(out_dir, "tables", "background_points.csv")
if (file.exists(background_path)) {
  background <- read.csv(background_path, stringsAsFactors = FALSE)
} else {
  set.seed(as.integer(get_cfg("seed", "123")))
  sampled <- raster::sampleRandom(current, size = min(10000L, raster::ncell(current)), xy = TRUE, na.rm = TRUE)
  background <- data.frame(lon = sampled[, ncol(sampled) - 1L], lat = sampled[, ncol(sampled)], stringsAsFactors = FALSE)
  occurrence_keys <- paste(occ$lon, occ$lat, sep = "|")
  background <- background[!(paste(background$lon, background$lat, sep = "|") %in% occurrence_keys), , drop = FALSE]
}

extract_values <- function(r, points) {
  values <- raster::extract(r, points[, c("lon", "lat")])
  as.numeric(values)
}

threshold_from_tss <- function(present, absent) {
  present <- present[is.finite(present)]
  absent <- absent[is.finite(absent)]
  if (!length(present) || !length(absent)) stop("Cannot select a threshold without valid occurrence and background predictions")
  candidates <- sort(unique(c(0, present, absent, 1)))
  scores <- vapply(candidates, function(thr) {
    sensitivity <- mean(present >= thr)
    specificity <- mean(absent < thr)
    sensitivity + specificity - 1
  }, numeric(1))
  # First threshold wins ties, keeping the result deterministic and no more
  # liberal than an equivalent later threshold.
  candidates[which.max(scores)]
}

if (threshold_raw == "auto") {
  threshold_value <- threshold_from_tss(extract_values(current, occ), extract_values(current, background))
}
if (!is.finite(threshold_value) || threshold_value <= 0 || threshold_value >= 1) stop("The selected threshold is outside (0, 1): ", threshold_value)

current_binary <- current >= threshold_value
current_binary[!is.finite(getValues(current))] <- NA
writeRaster(current_binary, file.path(range_raster_dir, "current_binary.tif"), datatype = "INT1U", NAflag = 255, overwrite = TRUE)

n_gcm <- length(gcms)
if (minimum_raw == "auto") {
  minimum_agreement <- switch(consensus, majority = floor(n_gcm / 2) + 1L, all = n_gcm, any = 1L, minimum = floor(n_gcm / 2) + 1L)
} else {
  minimum_agreement <- suppressWarnings(as.integer(minimum_raw))
  if (!is.finite(minimum_agreement) || minimum_agreement < 1L) stop("range_minimum_gcm_agreement must be auto or a positive integer")
}
if (minimum_agreement > n_gcm) stop("minimum_gcm_agreement cannot exceed the number of configured GCMs")
if (consensus == "all") minimum_agreement <- n_gcm
if (consensus == "any") minimum_agreement <- 1L

transition_from_binary <- function(current_binary, future_binary) {
  values <- cbind(getValues(current_binary), getValues(future_binary))
  classes <- rep(NA_integer_, nrow(values))
  valid <- is.finite(values[, 1]) & is.finite(values[, 2])
  classes[valid & values[, 1] < 0.5 & values[, 2] < 0.5] <- 0L
  classes[valid & values[, 1] >= 0.5 & values[, 2] < 0.5] <- 1L
  classes[valid & values[, 1] >= 0.5 & values[, 2] >= 0.5] <- 2L
  classes[valid & values[, 1] < 0.5 & values[, 2] >= 0.5] <- 3L
  out <- raster(current_binary)
  values(out) <- classes
  out
}

area_summary <- list()
row_id <- 0L
for (period in periods) for (scenario in scenarios) {
  future_paths <- file.path(raster_dir, paste0("future_suitability_", period, "_", scenario, "_", gcms, ".tif"))
  if (!all(file.exists(future_paths))) stop("Future suitability raster is missing for ", period, " ", scenario, ": ", paste(future_paths[!file.exists(future_paths)], collapse = ", "))
  future_layers <- lapply(future_paths, raster)
  for (layer in future_layers) {
    if (!compareRaster(current, layer, extent = TRUE, rowcol = TRUE, crs = TRUE, res = TRUE, orig = TRUE, rotation = TRUE, stopiffalse = FALSE)) stop("Current and future rasters do not share a grid for ", period, " ", scenario)
  }
  future_values <- do.call(cbind, lapply(future_layers, getValues))
  future_binary_values <- rep(NA_real_, nrow(future_values))
  # A transition cell is comparable only when every selected GCM has a
  # prediction there. This prevents a partially masked GCM set from receiving
  # an accidental majority vote.
  valid <- complete.cases(future_values)
  future_binary_values[valid] <- apply(future_values[valid, , drop = FALSE], 1, function(v) sum(v >= threshold_value, na.rm = TRUE) >= minimum_agreement)
  future_binary <- raster(future_layers[[1]])
  values(future_binary) <- future_binary_values
  transition <- transition_from_binary(current_binary, future_binary)

  prefix <- paste0(period, "_", scenario)
  writeRaster(future_binary, file.path(range_raster_dir, paste0("future_binary_", prefix, ".tif")), datatype = "INT1U", NAflag = 255, overwrite = TRUE)
  writeRaster(transition, file.path(range_raster_dir, paste0("transition_class_", prefix, ".tif")), datatype = "INT1U", NAflag = 255, overwrite = TRUE)

  cell_area <- raster::area(current)
  transition_values <- getValues(transition)
  cell_area_values <- getValues(cell_area)
  area_for_class <- function(class) sum(cell_area_values[is.finite(cell_area_values) & transition_values == class], na.rm = TRUE)
  lost <- area_for_class(1L)
  invariant <- area_for_class(2L)
  expanded <- area_for_class(3L)
  current_area <- sum(cell_area_values[is.finite(cell_area_values) & getValues(current_binary) >= 0.5], na.rm = TRUE)
  future_area <- sum(cell_area_values[is.finite(cell_area_values) & future_binary_values >= 0.5], na.rm = TRUE)
  union_area <- lost + invariant + expanded
  row_id <- row_id + 1L
  area_summary[[row_id]] <- data.frame(
    period = period, scenario = scenario, gcm = paste0("consensus_", minimum_agreement, "_of_", n_gcm),
    threshold = threshold_value, threshold_method = threshold_method, consensus = consensus,
    minimum_gcm_agreement = minimum_agreement, gcm_count = n_gcm, gcms = paste(gcms, collapse = ","),
    lost_km2 = lost, invariant_km2 = invariant, expanded_km2 = expanded,
    current_suitable_km2 = current_area, future_suitable_km2 = future_area, union_suitable_km2 = union_area,
    lost_pct_current = 100 * lost / current_area, invariant_pct_current = 100 * invariant / current_area,
    expanded_pct_current = 100 * expanded / current_area, future_minus_current_pct = 100 * (future_area - current_area) / current_area,
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, area_summary)
write.csv(summary_table, file.path(range_dir, "range_change_area_summary.csv"), row.names = FALSE)
write.table(data.frame(parameter = c("threshold", "threshold_method", "consensus", "minimum_gcm_agreement", "gcm_count"), value = c(threshold_value, threshold_method, consensus, minimum_agreement, n_gcm)), file.path(range_dir, "range_change_parameters.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Range-change analysis complete. Threshold:", format(threshold_value, digits = 8), "Agreement:", minimum_agreement, "of", n_gcm, "\n")
