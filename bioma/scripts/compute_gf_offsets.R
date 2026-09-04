#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 11) {
  stop("usage: compute_gf_offsets.R MODEL CURRENT LOCAL SEARCH OUTDIR SCENARIO RADII INITIAL_K BATCH VERIFY_N TIE_TOL")
}

suppressPackageStartupMessages(library(gradientForest))
suppressPackageStartupMessages(library(FNN))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(geosphere))

started <- proc.time()[["elapsed"]]
model_path <- args[[1]]
current_path <- args[[2]]
local_path <- args[[3]]
search_path <- args[[4]]
output_dir <- args[[5]]
scenario_name <- args[[6]]
radii_km <- as.numeric(strsplit(args[[7]], ",", fixed = TRUE)[[1]])
initial_k <- as.integer(args[[8]])
batch_size <- as.integer(args[[9]])
verify_n <- as.integer(args[[10]])
tie_tolerance <- as.numeric(args[[11]])

scenario_match <- regexec("^([0-9]{4}-[0-9]{4})-ssp([0-9]+)-(.+)$", scenario_name)
scenario_parts <- regmatches(scenario_name, scenario_match)[[1]]
if (length(scenario_parts) != 4) {
  stop("scenario must follow PERIOD-sspNNN-MODEL naming")
}
period_name <- scenario_parts[[2]]
ssp_name <- paste0("ssp", scenario_parts[[3]])
model_name <- scenario_parts[[4]]

loaded <- load(model_path)
if (!"all_gfmod" %in% loaded) {
  stop("model file must contain an object named all_gfmod")
}
predictors <- colnames(all_gfmod$X)
if (is.null(predictors) || length(predictors) == 0) {
  stop("cannot determine predictor names from Gradient Forest model")
}

current <- fread(current_path, na.strings = c("NA", "NaN", ""))
future_local <- fread(local_path, na.strings = c("NA", "NaN", ""))
future_search <- fread(search_path, na.strings = c("NA", "NaN", ""))
required <- c("lon", "lat", predictors)
for (item in list(current = current, future_local = future_local, future_search = future_search)) {
  if (!all(required %in% names(item))) {
    stop("a climate table is missing coordinates or model predictor columns")
  }
}
if (nrow(current) != nrow(future_local)) {
  stop("current and future-local tables have different row counts")
}
if (!isTRUE(all.equal(current$lon, future_local$lon, tolerance = 0)) ||
    !isTRUE(all.equal(current$lat, future_local$lat, tolerance = 0))) {
  stop("current and future-local coordinates are not identical")
}

complete_current <- complete.cases(current[, ..predictors])
complete_local <- complete.cases(future_local[, ..predictors])
complete_search <- complete.cases(future_search[, ..predictors])
if (!all(complete_current)) {
  stop("current background contains incomplete climate rows")
}
if (!any(complete_search)) {
  stop("future-search background contains no complete climate rows")
}

transform_gf <- function(table, complete_rows) {
  output <- matrix(NA_real_, nrow = nrow(table), ncol = length(predictors))
  colnames(output) <- predictors
  if (any(complete_rows)) {
    output[complete_rows, ] <- as.matrix(
      predict(all_gfmod, table[complete_rows, ..predictors])
    )
  }
  output
}

current_gf <- transform_gf(current, complete_current)
local_gf <- transform_gf(future_local, complete_local)
search_gf_all <- transform_gf(future_search, complete_search)
search_valid_rows <- which(complete_search)
search_gf <- search_gf_all[search_valid_rows, , drop = FALSE]
current_coords <- as.matrix(current[, .(lon, lat)])
local_coords <- as.matrix(future_local[, .(lon, lat)])
search_coords_all <- as.matrix(future_search[, .(lon, lat)])
search_coords <- search_coords_all[search_valid_rows, , drop = FALSE]

local_offset <- rep(NA_real_, nrow(current))
local_offset[complete_local] <- sqrt(
  rowSums((current_gf[complete_local, , drop = FALSE] - local_gf[complete_local, , drop = FALSE])^2)
)
local_output <- data.table(
  lon = current$lon,
  lat = current$lat,
  local_offset = local_offset,
  status = ifelse(complete_local, "matched", "future_nodata")
)
fwrite(local_output, file.path(output_dir, "local_offset.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")

make_k_sequence <- function(n_data, first_k) {
  values <- integer()
  k <- min(first_k, n_data)
  repeat {
    values <- c(values, k)
    if (k == n_data) break
    k <- min(n_data, max(k + 1L, k * 4L))
  }
  unique(values)
}

forward_exact <- function(data_gf, query_gf, data_coords, query_coords, radii_km,
                          initial_k, batch_size, tie_tolerance) {
  n_query <- nrow(query_gf)
  n_data <- nrow(data_gf)
  n_radius <- length(radii_km)
  selected_index <- matrix(NA_integer_, nrow = n_query, ncol = n_radius)
  selected_gf <- matrix(NA_real_, nrow = n_query, ncol = n_radius)
  selected_geo <- matrix(NA_real_, nrow = n_query, ncol = n_radius)
  unresolved <- matrix(TRUE, nrow = n_query, ncol = n_radius)
  k_sequence <- make_k_sequence(n_data, initial_k)
  maximum_k_used <- 0L
  max_pairs <- 2000000L

  for (k in k_sequence) {
    active_rows <- which(rowSums(unresolved) > 0)
    if (length(active_rows) == 0) break
    maximum_k_used <- k
    chunk_size <- max(1L, min(batch_size, floor(max_pairs / k)))
    chunks <- split(active_rows, ceiling(seq_along(active_rows) / chunk_size))
    for (rows in chunks) {
      neighbours <- get.knnx(data_gf, query_gf[rows, , drop = FALSE], k = k, algorithm = "kd_tree")
      indexes <- neighbours$nn.index
      distances <- neighbours$nn.dist
      if (is.null(dim(indexes))) indexes <- matrix(indexes, nrow = length(rows), byrow = TRUE)
      if (is.null(dim(distances))) distances <- matrix(distances, nrow = length(rows), byrow = TRUE)
      target_coords <- data_coords[as.vector(t(indexes)), , drop = FALSE]
      source_coords <- query_coords[rep(rows, each = k), , drop = FALSE]
      geographic <- matrix(
        distGeo(source_coords, target_coords),
        nrow = length(rows),
        ncol = k,
        byrow = TRUE
      )

      for (position in seq_along(rows)) {
        query_row <- rows[[position]]
        for (radius_index in which(unresolved[query_row, ])) {
          radius_m <- radii_km[[radius_index]] * 1000
          eligible <- if (is.infinite(radius_m)) {
            seq_len(k)
          } else {
            which(geographic[position, ] <= radius_m)
          }
          if (length(eligible) == 0) next
          best_distance <- min(distances[position, eligible])
          tied <- eligible[abs(distances[position, eligible] - best_distance) <= tie_tolerance]
          tied_geo <- geographic[position, tied]
          best_geo <- min(tied_geo)
          tied <- tied[abs(tied_geo - best_geo) <= 1e-7]
          best_position <- tied[which.min(indexes[position, tied])]

          boundary_clear <- k == n_data || distances[position, k] > best_distance + tie_tolerance
          if (boundary_clear) {
            selected_index[query_row, radius_index] <- indexes[position, best_position]
            selected_gf[query_row, radius_index] <- distances[position, best_position]
            selected_geo[query_row, radius_index] <- geographic[position, best_position]
            unresolved[query_row, radius_index] <- FALSE
          }
        }
      }
    }
  }
  list(index = selected_index, gf = selected_gf, geo = selected_geo, maximum_k = maximum_k_used)
}

forward <- forward_exact(
  search_gf,
  current_gf,
  search_coords,
  current_coords,
  radii_km,
  initial_k,
  batch_size,
  tie_tolerance
)

pair_bearing <- function(source, destination) {
  output <- rep(NA_real_, nrow(source))
  valid <- complete.cases(source) & complete.cases(destination)
  if (any(valid)) output[valid] <- bearing(source[valid, , drop = FALSE], destination[valid, , drop = FALSE])
  output
}

forward_tables <- vector("list", length(radii_km))
for (radius_index in seq_along(radii_km)) {
  selected <- forward$index[, radius_index]
  destination <- matrix(NA_real_, nrow = nrow(current), ncol = 2)
  matched <- !is.na(selected)
  destination[matched, ] <- search_coords[selected[matched], , drop = FALSE]
  forward_tables[[radius_index]] <- data.table(
    source_lon = current$lon,
    source_lat = current$lat,
    radius_km = if (is.infinite(radii_km[[radius_index]])) "unlimited" else format(radii_km[[radius_index]], trim = TRUE, scientific = FALSE),
    local_offset = local_offset,
    local_status = ifelse(complete_local, "matched", "future_nodata"),
    forward_offset = forward$gf[, radius_index],
    pred_distance_m = forward$geo[, radius_index],
    bearing_deg = pair_bearing(current_coords, destination),
    destination_lon = destination[, 1],
    destination_lat = destination[, 2],
    status = ifelse(matched, "matched", "no_candidate")
  )
}
forward_output <- rbindlist(forward_tables)
fwrite(forward_output, file.path(output_dir, "forward_offset.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")

reverse <- forward_exact(
  current_gf,
  search_gf,
  current_coords,
  search_coords,
  Inf,
  initial_k,
  batch_size,
  tie_tolerance
)
reverse_index <- reverse$index[, 1]
reverse_gf <- reverse$gf[, 1]
reverse_destination <- current_coords[reverse_index, , drop = FALSE]
reverse_geo <- reverse$geo[, 1]
reverse_bearing <- bearing(search_coords, reverse_destination)
reverse_output <- data.table(
  source_lon = future_search$lon,
  source_lat = future_search$lat,
  reverse_offset = NA_real_,
  pred_distance_m = NA_real_,
  bearing_deg = NA_real_,
  destination_lon = NA_real_,
  destination_lat = NA_real_,
  status = ifelse(complete_search, "matched", "future_nodata")
)
reverse_output[search_valid_rows, `:=`(
  reverse_offset = reverse_gf,
  pred_distance_m = reverse_geo,
  bearing_deg = reverse_bearing,
  destination_lon = reverse_destination[, 1],
  destination_lat = reverse_destination[, 2]
)]
fwrite(reverse_output, file.path(output_dir, "reverse_offset.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")

# One plot-ready final table. Reverse-only cells are retained when a future
# search mask is larger than the current mask, and are repeated per radius.
radius_table <- data.table(radius_km = vapply(
  radii_km,
  function(value) if (is.infinite(value)) "unlimited" else format(value, trim = TRUE, scientific = FALSE),
  character(1)
))
reverse_row_count <- nrow(reverse_output)
reverse_expanded <- reverse_output[rep(seq_len(reverse_row_count), each = nrow(radius_table))]
reverse_expanded[, radius_km := rep(radius_table$radius_km, times = reverse_row_count)]
setnames(
  reverse_expanded,
  c("source_lon", "source_lat", "pred_distance_m", "bearing_deg", "destination_lon", "destination_lat", "status"),
  c("lon", "lat", "reverse_pred_distance_m", "reverse_bearing_deg", "reverse_destination_lon", "reverse_destination_lat", "reverse_status")
)
forward_final <- copy(forward_output)
setnames(
  forward_final,
  c("source_lon", "source_lat", "pred_distance_m", "bearing_deg", "destination_lon", "destination_lat", "status"),
  c("lon", "lat", "forward_pred_distance_m", "forward_bearing_deg", "forward_destination_lon", "forward_destination_lat", "forward_status")
)
all_offsets <- merge(
  forward_final,
  reverse_expanded,
  by = c("lon", "lat", "radius_km"),
  all = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)
all_offsets[, `:=`(
  scenario = scenario_name,
  period = period_name,
  ssp = ssp_name,
  model = model_name
)]
setcolorder(
  all_offsets,
  c(
    "scenario", "period", "ssp", "model", "radius_km", "lon", "lat",
    "local_offset", "local_status",
    "forward_offset", "forward_pred_distance_m", "forward_bearing_deg",
    "forward_destination_lon", "forward_destination_lat", "forward_status",
    "reverse_offset", "reverse_pred_distance_m", "reverse_bearing_deg",
    "reverse_destination_lon", "reverse_destination_lat", "reverse_status"
  )
)
setorder(all_offsets, radius_km, lat, lon)
fwrite(all_offsets, file.path(output_dir, "all_offsets.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")

write_transformed <- function(table, transformed, filename) {
  output <- data.table(lon = table$lon, lat = table$lat)
  transformed_table <- as.data.table(transformed)
  setnames(transformed_table, paste0("gf_", predictors))
  fwrite(cbind(output, transformed_table), file.path(output_dir, filename), sep = "\t", compress = "gzip", na = "NA")
}
write_transformed(current, current_gf, "current_transformed.tsv.gz")
write_transformed(future_local, local_gf, "future_local_transformed.tsv.gz")
write_transformed(future_search, search_gf_all, "future_search_transformed.tsv.gz")

verification_cases <- 0L
if (verify_n > 0) {
  verify_rows <- unique(round(seq(1, nrow(current), length.out = min(verify_n, nrow(current)))))
  for (query_row in verify_rows) {
    brute_gf <- sqrt(rowSums((search_gf - matrix(current_gf[query_row, ], nrow(search_gf), ncol(search_gf), byrow = TRUE))^2))
    brute_geo <- distGeo(
      matrix(current_coords[query_row, ], nrow(search_coords), 2, byrow = TRUE),
      search_coords
    )
    for (radius_index in seq_along(radii_km)) {
      eligible <- if (is.infinite(radii_km[[radius_index]])) seq_len(nrow(search_gf)) else which(brute_geo <= radii_km[[radius_index]] * 1000)
      expected <- NA_integer_
      if (length(eligible)) {
        best_gf <- min(brute_gf[eligible])
        tied <- eligible[abs(brute_gf[eligible] - best_gf) <= tie_tolerance]
        best_geo <- min(brute_geo[tied])
        tied <- tied[abs(brute_geo[tied] - best_geo) <= 1e-7]
        expected <- min(tied)
      }
      observed <- forward$index[query_row, radius_index]
      if (!identical(as.integer(expected), as.integer(observed))) {
        stop(paste("forward nearest-neighbour verification failed at row", query_row, "radius", radii_km[[radius_index]]))
      }
      verification_cases <- verification_cases + 1L
    }
  }

  verify_search_rows <- unique(round(seq(1, nrow(search_gf), length.out = min(verify_n, nrow(search_gf)))))
  for (query_row in verify_search_rows) {
    brute_gf <- sqrt(rowSums((current_gf - matrix(search_gf[query_row, ], nrow(current_gf), ncol(current_gf), byrow = TRUE))^2))
    brute_geo <- distGeo(
      matrix(search_coords[query_row, ], nrow(current_coords), 2, byrow = TRUE),
      current_coords
    )
    best_gf <- min(brute_gf)
    tied <- which(abs(brute_gf - best_gf) <= tie_tolerance)
    best_geo <- min(brute_geo[tied])
    tied <- tied[abs(brute_geo[tied] - best_geo) <= 1e-7]
    expected <- min(tied)
    observed <- reverse$index[query_row, 1]
    if (!identical(as.integer(expected), as.integer(observed))) {
      stop(paste("reverse nearest-neighbour verification failed at valid-search row", query_row))
    }
    verification_cases <- verification_cases + 1L
  }
}

summarize_values <- function(type, radius, values) {
  valid <- values[is.finite(values)]
  data.table(
    offset_type = type,
    radius_km = radius,
    total_rows = length(values),
    valid_rows = length(valid),
    missing_rows = sum(!is.finite(values)),
    minimum = if (length(valid)) min(valid) else NA_real_,
    mean = if (length(valid)) mean(valid) else NA_real_,
    median = if (length(valid)) median(valid) else NA_real_,
    maximum = if (length(valid)) max(valid) else NA_real_
  )
}
summary_tables <- list(summarize_values("local", "same_cell", local_offset))
for (radius_index in seq_along(radii_km)) {
  label <- if (is.infinite(radii_km[[radius_index]])) "unlimited" else format(radii_km[[radius_index]], trim = TRUE, scientific = FALSE)
  summary_tables[[length(summary_tables) + 1L]] <- summarize_values("forward", label, forward$gf[, radius_index])
}
summary_tables[[length(summary_tables) + 1L]] <- summarize_values("reverse", "unlimited", reverse_output$reverse_offset)
fwrite(rbindlist(summary_tables), file.path(output_dir, "scenario_summary.tsv"), sep = "\t", na = "NA")

capture.output(sessionInfo(), file = file.path(output_dir, "r_session_info.txt"))
elapsed <- proc.time()[["elapsed"]] - started
summary <- data.table(
  summary_key = c(
    "r_version", "gradientforest_version", "fnn_version", "current_cells",
    "future_local_cells", "future_search_cells", "future_search_valid_cells", "combined_rows",
    "local_missing", "reverse_missing", "forward_unmatched", "maximum_k_used",
    "verification_cases", "elapsed_seconds"
  ),
  value = c(
    as.character(getRversion()), as.character(packageVersion("gradientForest")),
    as.character(packageVersion("FNN")), nrow(current), nrow(future_local),
    nrow(future_search), sum(complete_search), nrow(all_offsets), sum(!complete_local),
    sum(!complete_search), sum(is.na(forward$index)), max(forward$maximum_k, reverse$maximum_k),
    verification_cases, elapsed
  )
)
setnames(summary, "summary_key", "key")
fwrite(summary, file.path(output_dir, ".r_summary.tsv"), sep = "\t")
