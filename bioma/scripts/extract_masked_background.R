#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 22 || !args[[1]] %in% c("mask", "extract")) {
  stop("usage: extract_masked_background.R mask|extract MASK_OR_COORDS OUTPUT.tsv.gz BIO1.tif ... BIO19.tif")
}

suppressPackageStartupMessages(library(raster))

mode <- args[[1]]
input_path <- args[[2]]
output_path <- args[[3]]
raster_paths <- args[4:22]

climate <- stack(raster_paths)
names(climate) <- paste0("bio", seq_len(19))
if (mode == "mask") {
  mask_layer <- shapefile(input_path)
  if (is.na(crs(mask_layer))) {
    stop("mask has no coordinate reference system")
  }
  if (!compareCRS(climate, mask_layer)) {
    stop("mask and climate rasters use different coordinate reference systems")
  }
  masked <- mask(crop(climate, extent(mask_layer)), mask_layer)
  table <- as.data.frame(masked, xy = TRUE, na.rm = TRUE)
  names(table)[1:2] <- c("lon", "lat")
} else {
  reference <- read.delim(input_path, check.names = FALSE)
  if (!all(c("lon", "lat") %in% names(reference))) {
    stop("coordinate reference must contain lon and lat columns")
  }
  extracted <- extract(climate, reference[, c("lon", "lat")])
  table <- data.frame(reference[, c("lon", "lat")], extracted, check.names = FALSE)
  names(table) <- c("lon", "lat", paste0("bio", seq_len(19)))
}

connection <- gzfile(output_path, open = "wt")
write.table(
  table,
  file = connection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
close(connection)
