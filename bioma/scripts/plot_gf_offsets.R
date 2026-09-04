#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("usage: plot_gf_offsets.R OUTDIR RADIUS MODELS PERIOD SSP POPULATIONS BOUNDARY MIN_MODELS INPUT...")
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ragg))

output_dir <- args[[1]]
selected_radius <- args[[2]]
requested_models <- if (nzchar(args[[3]])) strsplit(args[[3]], ",", fixed = TRUE)[[1]] else character()
selected_period <- args[[4]]
selected_ssp <- args[[5]]
population_path <- args[[6]]
boundary_path <- args[[7]]
minimum_models <- as.integer(args[[8]])
input_paths <- args[9:length(args)]

fread_auto <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  magic <- readBin(connection, what = "raw", n = 2)
  is_gzip <- length(magic) == 2 && identical(as.integer(magic), c(31L, 139L))
  if (is_gzip) {
    return(fread(cmd = paste("gzip -dc --", shQuote(path)), na.strings = c("NA", "NaN", "")))
  }
  fread(path, na.strings = c("NA", "NaN", ""))
}

tables <- lapply(input_paths, fread_auto)
all_data <- rbindlist(tables, use.names = TRUE, fill = TRUE)
input_rows <- nrow(all_data)
required <- c("period", "ssp", "model", "radius_km", "lon", "lat", "local_offset", "forward_offset", "reverse_offset")
if (!all(required %in% names(all_data))) stop("an input final table is missing plotting columns")

all_data[, radius_km := as.character(radius_km)]
selected <- all_data[radius_km == selected_radius]
if (nzchar(selected_period)) selected <- selected[period == selected_period]
if (nzchar(selected_ssp)) selected <- selected[ssp == selected_ssp]
if (length(requested_models)) selected <- selected[model %in% requested_models]
if (nrow(selected) == 0) stop("no rows match the requested radius/period/SSP/models")

periods <- unique(selected$period)
ssps <- unique(selected$ssp)
if (length(periods) != 1) stop("plot selection must contain exactly one period; use --period")
if (length(ssps) != 1) stop("plot selection must contain exactly one SSP; use --ssp")
present_models <- sort(unique(selected$model))
if (length(requested_models) && !all(requested_models %in% present_models)) {
  stop(paste("requested models are missing:", paste(setdiff(requested_models, present_models), collapse = ",")))
}
if (length(present_models) < minimum_models) {
  stop(paste0(
    "ensemble plotting requires at least ", minimum_models,
    " distinct models; found ", length(present_models)
  ))
}
duplicate_rows <- selected[, .N, by = .(period, ssp, radius_km, model, lon, lat)][N > 1]
if (nrow(duplicate_rows)) {
  stop("duplicate model-coordinate rows would give a model extra weight in the ensemble mean")
}
coordinate_model_counts <- selected[, .(n_models = uniqueN(model)), by = .(period, ssp, radius_km, lon, lat)]
if (any(coordinate_model_counts$n_models != length(present_models))) {
  stop("selected models do not have an identical coordinate grid")
}

strict_ensemble_mean <- function(values) {
  if (sum(!is.na(values)) != length(present_models)) return(NA_real_)
  mean(values)
}
ensemble_data <- selected[, .(
  local_offset_mean = strict_ensemble_mean(local_offset),
  forward_offset_mean = strict_ensemble_mean(forward_offset),
  reverse_offset_mean = strict_ensemble_mean(reverse_offset),
  n_models_local = sum(!is.na(local_offset)),
  n_models_forward = sum(!is.na(forward_offset)),
  n_models_reverse = sum(!is.na(reverse_offset))
), by = .(period, ssp, radius_km, lon, lat)]
ensemble_data[, `:=`(
  ensemble_models = paste(present_models, collapse = ","),
  n_models_requested = length(present_models)
)]
setcolorder(
  ensemble_data,
  c(
    "period", "ssp", "radius_km", "lon", "lat", "ensemble_models", "n_models_requested",
    "local_offset_mean", "n_models_local", "forward_offset_mean", "n_models_forward",
    "reverse_offset_mean", "n_models_reverse"
  )
)
fwrite(
  ensemble_data,
  file.path(output_dir, "ensemble_mean_offsets.tsv.gz"),
  sep = "\t", compress = "gzip", na = "NA"
)
plot_data <- copy(ensemble_data)

normalize01 <- function(values) {
  output <- rep(NA_real_, length(values))
  valid <- is.finite(values)
  if (!any(valid)) return(output)
  limits <- range(values[valid])
  if (diff(limits) == 0) {
    output[valid] <- 0
  } else {
    output[valid] <- (values[valid] - limits[[1]]) / diff(limits)
  }
  output
}
offset_bin <- function(values) {
  cut(
    values,
    breaks = c(-Inf, 0.25, 0.5, 0.75, Inf),
    labels = c("0–0.25", "0.25–0.50", "0.50–0.75", "0.75–1.00"),
    include.lowest = TRUE,
    right = TRUE
  )
}
plot_data[, `:=`(
  local_norm = normalize01(local_offset_mean),
  forward_norm = normalize01(forward_offset_mean),
  reverse_norm = normalize01(reverse_offset_mean)
)]
plot_data[, `:=`(
  local_bin = offset_bin(local_norm),
  forward_bin = offset_bin(forward_norm),
  reverse_bin = offset_bin(reverse_norm)
)]
complete_rgb <- complete.cases(plot_data[, .(local_norm, forward_norm, reverse_norm)])
plot_data[, color_rgb := NA_character_]
plot_data[complete_rgb, color_rgb := rgb(local_norm, forward_norm, reverse_norm)]
fwrite(plot_data, file.path(output_dir, "plot_data.tsv.gz"), sep = "\t", compress = "gzip", na = "NA")

boundary <- NULL
if (nzchar(boundary_path)) {
  boundary <- suppressMessages(st_read(boundary_path, quiet = TRUE))
  if (is.na(st_crs(boundary))) stop("boundary has no coordinate reference system")
  boundary <- st_transform(boundary, 4326)
}
pops <- NULL
pop_labels <- NULL
if (nzchar(population_path)) {
  pops <- fread(population_path)
  lon_name <- intersect(c("lon", "longitude"), names(pops))[[1]]
  lat_name <- intersect(c("lat", "latitude"), names(pops))[[1]]
  label_name <- intersect(c("group", "pop", "population_id", "ID"), names(pops))[[1]]
  setnames(pops, c(lon_name, lat_name, label_name), c("lon", "lat", "label"))
  pop_labels <- pops[, .(lon = mean(lon), lat = mean(lat)), by = label]
}

x_limits <- range(plot_data$lon, na.rm = TRUE) + c(-0.25, 0.25)
y_limits <- range(plot_data$lat, na.rm = TRUE) + c(-0.25, 0.25)
palette <- c(
  "0–0.25" = "#C1E0BA",
  "0.25–0.50" = "#6EA5C9",
  "0.50–0.75" = "#F7E4AC",
  "0.75–1.00" = "#D8747D"
)
context <- paste0(
  periods[[1]], " | ", ssps[[1]], " | ", selected_radius,
  "\nModels: ", paste(present_models, collapse = " + ")
)

map_base <- function() {
  plot <- ggplot()
  plot <- plot + coord_sf(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "#4B5563"),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.1, "cm")
    ) + labs(x = "Longitude", y = "Latitude")
  plot
}

add_map_context <- function(plot) {
  if (!is.null(boundary)) plot <- plot + geom_sf(data = boundary, fill = NA, color = "#222222", linewidth = 0.35)
  if (!is.null(pops)) {
    plot <- plot + geom_point(data = pops, aes(lon, lat), shape = 21, fill = "white", color = "black", size = 1.8, stroke = 0.35)
    plot <- plot + geom_text(data = pop_labels, aes(lon, lat, label = label), size = 2.5, color = "#222222", check_overlap = TRUE)
  }
  plot
}

composite_data <- plot_data[complete_rgb]
p1 <- map_base() +
  geom_point(data = composite_data, aes(lon, lat, color = color_rgb), shape = 15, size = 0.75) +
  scale_color_identity() +
  labs(
    title = "Ensemble-mean GF offsets",
    subtitle = context,
    caption = paste(
      "Means require complete model coverage.",
      "RGB: local = red, forward = green, reverse = blue.",
      "Channels are min–max normalized.",
      sep = "\n"
    )
  )
p1 <- add_map_context(p1)

offset_map <- function(value_column, bin_column, title) {
  plot <- map_base() +
    geom_point(
      data = plot_data[!is.na(get(value_column))],
      aes(x = lon, y = lat, color = get(bin_column)),
      shape = 15,
      size = 0.75
    ) +
    scale_color_manual(values = palette, drop = FALSE, name = "Normalized offset") +
    labs(title = title, subtitle = paste0(context, "\nScale: min–max normalized"))
  add_map_context(plot)
}
p2 <- offset_map("local_offset_mean", "local_bin", "Ensemble-mean local offset")
p3 <- offset_map("forward_offset_mean", "forward_bin", "Ensemble-mean forward offset")
p4 <- offset_map("reverse_offset_mean", "reverse_bin", "Ensemble-mean reverse offset")

map_grid <- plot_grid(p1, p2, p3, p4, ncol = 2, labels = c("A", "B", "C", "D"), align = "hv")
ggsave(file.path(output_dir, "all_offsets_maps.pdf"), map_grid, width = 15, height = 11, device = cairo_pdf)
ggsave(file.path(output_dir, "all_offsets_maps.png"), map_grid, width = 15, height = 11, dpi = 300, device = ragg::agg_png, background = "white")

relationship_plot <- function(x_column, y_column, x_label, y_label) {
  data <- plot_data[complete_rgb & complete.cases(plot_data[, c(x_column, y_column), with = FALSE])]
  limits <- range(c(data[[x_column]], data[[y_column]]), na.rm = TRUE)
  ggplot(data, aes(x = get(x_column), y = get(y_column))) +
    geom_point(aes(color = color_rgb), size = 0.6, alpha = 0.45) +
    scale_color_identity() +
    geom_abline(intercept = 0, slope = 1, color = "#222222", linewidth = 0.45, linetype = 2) +
    coord_equal(xlim = limits, ylim = limits) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "none") +
    labs(x = x_label, y = y_label, subtitle = paste0("n = ", nrow(data), " cells; dashed line = 1:1"))
}
p5 <- relationship_plot("local_offset_mean", "forward_offset_mean", "Mean local offset", "Mean forward offset")
p6 <- relationship_plot("local_offset_mean", "reverse_offset_mean", "Mean local offset", "Mean reverse offset")
p7 <- relationship_plot("forward_offset_mean", "reverse_offset_mean", "Mean forward offset", "Mean reverse offset")
relationship_grid <- plot_grid(p5, p6, p7, ncol = 3, labels = c("E", "F", "G"), align = "hv")
ggsave(file.path(output_dir, "offset_relationships.pdf"), relationship_grid, width = 15, height = 5, device = cairo_pdf)
ggsave(file.path(output_dir, "offset_relationships.png"), relationship_grid, width = 15, height = 5, dpi = 300, device = ragg::agg_png, background = "white")

summary <- data.frame(
  key = c(
    "r_version", "ggplot2_version", "sf_version", "input_rows", "selected_rows", "plot_cells", "models",
    "full_local_coverage_cells", "full_forward_coverage_cells", "full_reverse_coverage_cells",
    "period", "ssp", "radius", "model_names"
  ),
  value = c(
    as.character(getRversion()), as.character(packageVersion("ggplot2")), as.character(packageVersion("sf")),
    input_rows, nrow(selected), nrow(plot_data), length(present_models),
    sum(plot_data$n_models_local == length(present_models)),
    sum(plot_data$n_models_forward == length(present_models)),
    sum(plot_data$n_models_reverse == length(present_models)),
    periods[[1]], ssps[[1]], selected_radius,
    paste(present_models, collapse = ",")
  ),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(output_dir, ".plot_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "r_session_info.txt"))
