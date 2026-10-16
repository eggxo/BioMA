#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: niche_plot.R output_dir occurrence_csv mask_shp")

suppressPackageStartupMessages({
  library(raster)
  library(sf)
  library(ggplot2)
})
root <- normalizePath(args[[1]], mustWork = TRUE)
occurrence_csv <- args[[2]]
mask_path <- args[[3]]
raster_dir <- file.path(root, "rasters")
fig_dir <- file.path(root, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

mask_sf <- st_read(mask_path, quiet = TRUE)
if (is.na(st_crs(mask_sf))) {
  st_crs(mask_sf) <- 4326
} else if (!identical(st_crs(mask_sf)$epsg, 4326L)) {
  mask_sf <- st_transform(mask_sf, 4326)
}
if (any(!st_is_valid(mask_sf))) mask_sf <- st_make_valid(mask_sf)
mask_polys <- st_cast(st_geometry(mask_sf), "POLYGON")
mask_outline <- do.call(rbind, lapply(seq_along(mask_polys), function(i) {
  xy <- st_coordinates(mask_polys[[i]])
  data.frame(x = xy[, "X"], y = xy[, "Y"], group = paste(i, xy[, "L1"], sep = "_"))
}))

occ <- read.csv(occurrence_csv, stringsAsFactors = FALSE)
occ <- occ[is.finite(as.numeric(occ$lon)) & is.finite(as.numeric(occ$lat)), , drop = FALSE]
occ$lon <- as.numeric(occ$lon)
occ$lat <- as.numeric(occ$lat)
group_palette <- function(values) {
  groups <- sort(unique(as.character(values)))
  groups <- groups[nzchar(groups) & !is.na(groups)]
  if (!length(groups)) return(setNames(character(), character()))
  setNames(grDevices::hcl.colors(length(groups), palette = "Dark 3"), groups)
}

bb <- st_bbox(mask_sf)
xlim <- c(as.numeric(bb[["xmin"]]) - 0.5, as.numeric(bb[["xmax"]]) + 0.5)
ylim <- c(as.numeric(bb[["ymin"]]) - 0.5, as.numeric(bb[["ymax"]]) + 0.5)

raster_df <- function(path) {
  r <- raster::raster(path)
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (!nrow(d)) return(NULL)
  names(d)[3] <- "value"
  d
}

plot_map <- function(path, kind, title) {
  d <- raster_df(path)
  if (is.null(d)) return(NULL)
  p <- ggplot() +
    geom_raster(data = d, aes(x = x, y = y, fill = value)) +
    geom_path(data = mask_outline, aes(x = x, y = y, group = group), color = "black", linewidth = 0.35) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, x = "Longitude", y = "Latitude")
  if (kind == "maladaptation") {
    p <- p + scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, name = "Current - future")
  } else {
    p <- p + scale_fill_viridis_c(option = "plasma", name = "Suitability")
  }
  if ("cluster" %in% names(occ)) {
    occ$cluster <- as.character(occ$cluster)
    cols <- group_palette(occ$cluster)
    p <- p + geom_point(data = occ, aes(lon, lat, color = cluster), shape = 21, fill = "white", size = 2.2, stroke = 0.35) + scale_color_manual(values = cols, guide = "none")
  } else {
    p <- p + geom_point(data = occ, aes(lon, lat), shape = 21, fill = "white", color = "black", size = 2.2, stroke = 0.35)
  }
  p + theme_bw() + theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold", size = 12), legend.position = "right")
}

save_plot <- function(p, stem, width = 7, height = 5.5) {
  if (is.null(p)) return(invisible(NULL))
  ggsave(file.path(fig_dir, paste0(stem, ".png")), p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), p, width = width, height = height, bg = "white")
}

current_path <- file.path(raster_dir, "current_suitability.tif")
if (file.exists(current_path)) save_plot(plot_map(current_path, "suitability", "Current MaxEnt suitability"), "current_suitability")

future_files <- list.files(raster_dir, pattern = "^future_suitability_.*_ensemble_(mean|median)\\.tif$", full.names = TRUE)
future_plots <- list()
for (path in sort(future_files)) {
  label <- sub("^future_suitability_", "", sub("_ensemble_(mean|median)\\.tif$", "", basename(path)))
  stem <- paste0("future_suitability_", gsub("[^A-Za-z0-9]+", "_", label))
  p <- plot_map(path, "suitability", paste0("Future MaxEnt suitability (", gsub("_", " | ", label), ")"))
  save_plot(p, stem)
  future_plots[[length(future_plots) + 1L]] <- p
}

delta_files <- list.files(raster_dir, pattern = "^maladaptation_.*_ensemble_(mean|median)\\.tif$", full.names = TRUE)
delta_plots <- list()
for (path in sort(delta_files)) {
  label <- sub("^maladaptation_", "", sub("_ensemble_(mean|median)\\.tif$", "", basename(path)))
  stem <- paste0("maladaptation_", gsub("[^A-Za-z0-9]+", "_", label))
  p <- plot_map(path, "maladaptation", paste0("Maladaptation (current - future; ", gsub("_", " | ", label), ")"))
  save_plot(p, stem)
  delta_plots[[length(delta_plots) + 1L]] <- p
}

if (requireNamespace("cowplot", quietly = TRUE)) {
  if (length(delta_plots)) {
    combo <- cowplot::plot_grid(plotlist = delta_plots, ncol = 2, labels = LETTERS[seq_along(delta_plots)])
    ggsave(file.path(fig_dir, "Figure_Niche_Maladaptation_Maps.png"), combo, width = 13, height = 5.8 * ceiling(length(delta_plots) / 2), dpi = 300, bg = "white")
    ggsave(file.path(fig_dir, "Figure_Niche_Maladaptation_Maps.pdf"), combo, width = 13, height = 5.8 * ceiling(length(delta_plots) / 2), bg = "white")
  }
  if (length(future_plots)) {
    combo <- cowplot::plot_grid(plotlist = future_plots, ncol = 2, labels = LETTERS[seq_along(future_plots)])
    ggsave(file.path(fig_dir, "Figure_Niche_Future_Suitability_Maps.png"), combo, width = 13, height = 5.8 * ceiling(length(future_plots) / 2), dpi = 300, bg = "white")
    ggsave(file.path(fig_dir, "Figure_Niche_Future_Suitability_Maps.pdf"), combo, width = 13, height = 5.8 * ceiling(length(future_plots) / 2), bg = "white")
  }
}

metrics_path <- file.path(root, "tuning_metrics.csv")
if (file.exists(metrics_path)) {
  metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
  metrics <- metrics[is.finite(metrics$auc_mean) | is.finite(metrics$tss_mean), , drop = FALSE]
  if (nrow(metrics)) {
    p_auc <- ggplot(metrics, aes(x = beta_multiplier, y = auc_mean, color = feature_class, shape = factor(n_variables))) +
      geom_point(size = 2.8, alpha = 0.85) +
      geom_line(aes(group = interaction(feature_class, n_variables)), alpha = 0.45) +
      labs(x = "Beta multiplier", y = "Cross-validated AUC", color = "Feature class", shape = "Number of variables", title = "MaxEnt variable and parameter tuning") +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_plot(p_auc, "Figure_Niche_Tuning", width = 9, height = 5.5)
  }
}

# Optional binary habitat transition products generated by niche_range_change.R.
# These maps are kept separate from continuous suitability and maladaptation.
range_dir <- file.path(root, "range_change")
range_raster_dir <- file.path(range_dir, "rasters")
range_files <- if (dir.exists(range_raster_dir)) list.files(range_raster_dir, pattern = "^transition_class_.*\\.tif$", full.names = TRUE) else character()
if (length(range_files)) {
  transition_labels <- c("0" = "Unsuitable in both", "1" = "Contract", "2" = "Invariant", "3" = "Expand")
  transition_colors <- c("0" = "#FFFFFF", "1" = "#C97979", "2" = "#8A6F9E", "3" = "#4F91C2")
  transition_rows <- lapply(range_files, function(path) {
    d <- raster_df(path)
    if (is.null(d)) return(NULL)
    d$class <- factor(as.character(round(d$value)), levels = names(transition_colors))
    d$scenario <- sub("^transition_class_", "", sub("\\.tif$", "", basename(path)))
    d
  })
  transition_data <- do.call(rbind, transition_rows[!vapply(transition_rows, is.null, logical(1))])
  if (nrow(transition_data)) {
    transition_plot <- ggplot(transition_data, aes(x = x, y = y, fill = class)) +
      geom_raster(interpolate = FALSE) +
      geom_path(data = mask_outline, aes(x = x, y = y, group = group), inherit.aes = FALSE, color = "black", linewidth = 0.35) +
      facet_wrap(~scenario, ncol = 2) +
      scale_fill_manual(values = transition_colors, labels = transition_labels, drop = FALSE, name = "Habitat transition") +
      coord_equal(expand = FALSE) +
      labs(title = "Habitat contraction, persistence, and expansion", x = "Longitude", y = "Latitude") +
      theme_bw() +
      theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right")
    save_plot(transition_plot, "Figure_Niche_Range_Change_Maps", width = 12, height = 5.8 * ceiling(length(range_files) / 2))
  }

  area_path <- file.path(range_dir, "range_change_area_summary.csv")
  if (file.exists(area_path)) {
    area <- read.csv(area_path, stringsAsFactors = FALSE)
    if (nrow(area)) {
      area_long <- rbind(
        data.frame(period = area$period, scenario = area$scenario, category = "Contract", percent = area$lost_pct_current),
        data.frame(period = area$period, scenario = area$scenario, category = "Invariant", percent = area$invariant_pct_current),
        data.frame(period = area$period, scenario = area$scenario, category = "Expand", percent = area$expanded_pct_current)
      )
      area_long$scenario <- paste(area_long$period, area_long$scenario)
      area_plot <- ggplot(area_long, aes(x = scenario, y = percent, fill = category)) +
        geom_col(position = "dodge", width = 0.72) +
        scale_fill_manual(values = c(Contract = "#C97979", Invariant = "#8A6F9E", Expand = "#4F91C2")) +
        labs(title = "Habitat transition area relative to current suitable area", x = NULL, y = "Area (%)", fill = NULL) +
        theme_bw() +
        theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(hjust = 0.5, face = "bold"))
      save_plot(area_plot, "Figure_Niche_Range_Change_Area", width = 10, height = 5.5)
    }
  }
}
cat("Niche figures written to ", fig_dir, "\n")
