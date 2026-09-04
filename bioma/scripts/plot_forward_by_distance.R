#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("usage: plot_forward_by_distance.R OUTDIR MODELS PERIOD SSPS MIN_MODELS INPUT...")
}

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ragg))

output_dir <- args[[1]]
requested_models <- if (nzchar(args[[2]])) strsplit(args[[2]], ",", fixed = TRUE)[[1]] else character()
selected_period <- args[[3]]
requested_ssps <- if (nzchar(args[[4]])) strsplit(args[[4]], ",", fixed = TRUE)[[1]] else character()
minimum_models <- as.integer(args[[5]])
input_paths <- args[6:length(args)]

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
required <- c("period", "ssp", "model", "radius_km", "lon", "lat", "forward_offset")
if (!all(required %in% names(all_data))) stop("an input final table is missing forward-distance columns")

all_data[, ssp := tolower(as.character(ssp))]
all_data[grepl("^[0-9]+$", ssp), ssp := paste0("ssp", ssp)]
all_data[, radius_key := tolower(as.character(radius_km))]
all_data[radius_key %in% c("inf", "infinity", "all", "unlimit"), radius_key := "unlimited"]
all_data[, radius_key := sub("km$", "", radius_key)]

selected <- all_data
if (nzchar(selected_period)) selected <- selected[period == selected_period]
if (length(requested_ssps)) selected <- selected[ssp %in% requested_ssps]
if (length(requested_models)) selected <- selected[model %in% requested_models]
if (nrow(selected) == 0) stop("no rows match the requested period/SSPs/models")
finite_radius_rows <- selected$radius_key != "unlimited"
finite_radius_values <- suppressWarnings(as.numeric(selected$radius_key[finite_radius_rows]))
if (any(!is.finite(finite_radius_values)) || any(finite_radius_values <= 0)) {
  stop("migration radii must be positive numbers or unlimited")
}
selected[finite_radius_rows, radius_key := format(
  finite_radius_values, scientific = FALSE, trim = TRUE
)]

periods <- unique(selected$period)
if (length(periods) != 1) stop("forward-distance plot requires exactly one period; use --period")
present_models <- sort(unique(selected$model))
if (length(requested_models) && !all(requested_models %in% present_models)) {
  stop(paste("requested models are missing:", paste(setdiff(requested_models, present_models), collapse = ",")))
}
if (length(present_models) < minimum_models) {
  stop(paste0(
    "forward-distance plotting requires at least ", minimum_models,
    " distinct models; found ", length(present_models)
  ))
}
present_ssps <- unique(selected$ssp)
ssp_numbers <- suppressWarnings(as.numeric(sub("^ssp", "", present_ssps)))
present_ssps <- present_ssps[order(is.na(ssp_numbers), ssp_numbers, present_ssps)]
if (length(requested_ssps) && !all(requested_ssps %in% present_ssps)) {
  stop(paste("requested SSPs are missing:", paste(setdiff(requested_ssps, present_ssps), collapse = ",")))
}

duplicate_rows <- selected[, .N, by = .(period, ssp, radius_key, model, lon, lat)][N > 1]
if (nrow(duplicate_rows)) {
  stop("duplicate model-coordinate-radius rows would give a model extra weight")
}
coordinate_model_counts <- selected[, .(n_models = uniqueN(model)), by = .(period, ssp, radius_key, lon, lat)]
if (any(coordinate_model_counts$n_models != length(present_models))) {
  stop("selected models do not have identical coordinate grids for every SSP and radius")
}

strict_ensemble_mean <- function(values) {
  if (sum(!is.na(values)) != length(present_models)) return(NA_real_)
  mean(values)
}
ensemble_grid <- selected[, .(
  forward_offset_mean = strict_ensemble_mean(forward_offset),
  n_models_forward = sum(!is.na(forward_offset))
), by = .(period, ssp, radius_key, lon, lat)]
ensemble_grid[, `:=`(
  ensemble_models = paste(present_models, collapse = ","),
  n_models_requested = length(present_models)
)]
setcolorder(
  ensemble_grid,
  c(
    "period", "ssp", "radius_key", "lon", "lat", "ensemble_models", "n_models_requested",
    "forward_offset_mean", "n_models_forward"
  )
)
fwrite(
  ensemble_grid,
  file.path(output_dir, "forward_offset_ensemble_grid.tsv.gz"),
  sep = "\t", compress = "gzip", na = "NA"
)

summary_data <- ensemble_grid[!is.na(forward_offset_mean), .(
  median_offset = median(forward_offset_mean),
  q25 = as.numeric(quantile(forward_offset_mean, 0.25, names = FALSE)),
  q75 = as.numeric(quantile(forward_offset_mean, 0.75, names = FALSE)),
  n_cells = .N
), by = .(period, ssp, radius_key)]
if (nrow(summary_data) == 0) stop("no complete-model forward offsets are available to summarize")

radius_keys <- unique(summary_data$radius_key)
has_unlimited <- "unlimited" %in% radius_keys
finite_values <- suppressWarnings(as.numeric(setdiff(radius_keys, "unlimited")))
finite_values <- sort(unique(finite_values))
if (!length(finite_values)) stop("at least one finite migration radius is required")
radius_layout <- data.table(
  radius_key = format(finite_values, scientific = FALSE, trim = TRUE),
  x_pos = seq_along(finite_values),
  radius_label = format(finite_values, scientific = FALSE, trim = TRUE)
)
if (has_unlimited) {
  radius_layout <- rbind(
    radius_layout,
    data.table(
      radius_key = "unlimited",
      x_pos = length(finite_values) + 1,
      radius_label = "Unlimited"
    )
  )
}
summary_data <- radius_layout[summary_data, on = "radius_key"]
if (any(is.na(summary_data$x_pos))) stop("failed to position one or more migration radii")
summary_data[, ssp_label := toupper(ssp)]
setorder(summary_data, ssp, x_pos)
fwrite(summary_data, file.path(output_dir, "forward_offset_by_distance.tsv"), sep = "\t", na = "NA")

if (length(present_ssps) > 4) stop("at most four SSP series can be shown legibly in one panel")
known_colors <- c(
  "SSP126" = "#4C78A8",
  "SSP245" = "#2878D0",
  "SSP370" = "#D49A24",
  "SSP585" = "#C44700"
)
fallback_colors <- c("#4C78A8", "#2878D0", "#D49A24", "#C44700")
series_labels <- toupper(present_ssps)
series_colors <- setNames(fallback_colors[seq_along(series_labels)], series_labels)
known <- intersect(series_labels, names(known_colors))
series_colors[known] <- known_colors[known]
series_linetypes <- setNames(c("solid", "dashed", "dotdash", "longdash")[seq_along(series_labels)], series_labels)
series_shapes <- setNames(c(21, 22, 24, 23)[seq_along(series_labels)], series_labels)

y_limits <- range(c(summary_data$q25, summary_data$q75), na.rm = TRUE)
y_span <- diff(y_limits)
if (!is.finite(y_span) || y_span == 0) y_span <- max(abs(y_limits), 1) * 0.1
y_plot_limits <- c(max(0, y_limits[[1]] - 0.08 * y_span), y_limits[[2]] + 0.08 * y_span)
cell_counts <- range(summary_data$n_cells)
cell_text <- if (cell_counts[[1]] == cell_counts[[2]]) {
  paste0(format(cell_counts[[1]], big.mark = ","), " cells per point")
} else {
  paste0(format(cell_counts[[1]], big.mark = ","), "–", format(cell_counts[[2]], big.mark = ","), " cells per point")
}
subtitle <- paste0(
  periods[[1]], " | ", length(present_models), "-model ensemble mean",
  " | median and IQR across ", cell_text
)

p <- ggplot(
  summary_data,
  aes(x = x_pos, y = median_offset, color = ssp_label, group = ssp_label)
) +
  geom_ribbon(aes(ymin = q25, ymax = q75, fill = ssp_label), alpha = 0.14, color = NA) +
  geom_line(aes(linetype = ssp_label), linewidth = 1.2) +
  geom_point(aes(shape = ssp_label), size = 3.5, fill = "white", stroke = 0.9) +
  scale_color_manual(values = series_colors, breaks = series_labels, name = NULL) +
  scale_fill_manual(values = series_colors, breaks = series_labels, guide = "none") +
  scale_linetype_manual(values = series_linetypes, breaks = series_labels, name = NULL) +
  scale_shape_manual(values = series_shapes, breaks = series_labels, name = NULL) +
  scale_x_continuous(
    breaks = radius_layout$x_pos,
    labels = radius_layout$radius_label,
    expand = expansion(mult = c(0.03, 0.04))
  ) +
  coord_cartesian(ylim = y_plot_limits, clip = "off") +
  labs(
    title = "Forward offset across migration-distance limits",
    subtitle = subtitle,
    x = "Maximum migration-distance limit",
    y = "Forward offset",
    caption = "Line = median; ribbon = IQR. Distance limits are shown as ordered categories."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
    panel.border = element_rect(color = "#30343B", fill = NA, linewidth = 0.6),
    axis.title = element_text(color = "#20242A"),
    axis.text = element_text(color = "#30343B"),
    plot.title = element_text(face = "bold", color = "#20242A"),
    plot.subtitle = element_text(color = "#4B5563"),
    plot.caption = element_text(color = "#5F6773", hjust = 0),
    legend.position = "top",
    legend.justification = "right",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(1.1, "cm"),
    plot.margin = margin(8, 14, 8, 8)
  )

ggsave(
  file.path(output_dir, "forward_offset_by_distance.pdf"),
  p, width = 7.2, height = 5.3, device = cairo_pdf
)
ggsave(
  file.path(output_dir, "forward_offset_by_distance.png"),
  p, width = 7.2, height = 5.3, dpi = 600,
  device = ragg::agg_png, background = "white"
)

summary <- data.frame(
  key = c(
    "r_version", "ggplot2_version", "input_rows", "selected_rows", "ensemble_grid_rows",
    "valid_ensemble_grid_rows", "summary_rows", "models", "ssps", "radii",
    "period", "model_names", "ssp_names", "radius_names"
  ),
  value = c(
    as.character(getRversion()), as.character(packageVersion("ggplot2")), input_rows, nrow(selected),
    nrow(ensemble_grid), sum(!is.na(ensemble_grid$forward_offset_mean)), nrow(summary_data),
    length(present_models), length(present_ssps), nrow(radius_layout), periods[[1]],
    paste(present_models, collapse = ","), paste(present_ssps, collapse = ","),
    paste(radius_layout$radius_label, collapse = ",")
  ),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(output_dir, ".plot_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "r_session_info.txt"))
