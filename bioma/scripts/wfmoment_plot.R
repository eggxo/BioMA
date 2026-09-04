#!/usr/bin/env Rscript

# Parameterised publication-style plot for the two-dimensional WFmoment
# habitat-loss workflow.  The script deliberately keeps input values intact:
# any interpolation or optional terminal tail is recorded as display-only data.

options(stringsAsFactors = FALSE)

required_packages <- c("ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "),
       ". Install them before running wfmoment_plot.R.", call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript wfmoment_plot.R --curve-file CURVE.tsv [options]\n\n",
    "Required:\n",
    "  --curve-file PATH       WFmoment curve summary (TSV/CSV, optionally .gz)\n",
    "  --point-file PATH       scenario lookup points (TSV/CSV, optionally .gz)\n\n",
    "Optional inputs:\n",
    "  --short-file PATH       shrt_time_pi_loss.csv (used when GDAR values are absent)\n",
    "  --scenario-label-file PATH  two-column scenario,label mapping\n",
    "  --meta-file PATH        metadata table containing z_gdar (optional)\n\n",
    "Output and display:\n",
    "  --outdir DIR            output directory (default: curve directory)\n",
    "  --prefix NAME           output stem (default: wfmoment_2d)\n",
    "  --include-equilibrium  add equilibrium line/points when available\n",
    "  --include-zero-deme-endpoint  keep occupied_demes=0 endpoint in the main line\n",
    "  --all-wf-points        plot WF3/WF5 points for every scenario (default: period-matched)\n",
    "  --tail-mode MODE        auto, none, smooth, zero, or linear; display only\n",
    "  --tail-cutoff-demes N   start the display tail at the last state with N demes (default: 5)\n",
    "  --min-display-demes N   minimum extant demes in a continuous main line (default: 5)\n",
    "  --tail-end NUMBER      display tail endpoint (default: x-max)\n",
    "  --tail-exponent NUMBER curvature power for smooth terminal tail (default: 1.5)\n",
    "  --gdar-display-epsilon NUMBER  stop GDAR theory at 100*(1-epsilon)% (default: 0.001)\n",
    "  --display-monotone MODE auto, none, or cummin; auto applies to edge loss\n",
    "  --line-method MODE      raw, linear, or spline (default: spline; raw when --plot-raw=true)\n",
    "  --n-grid INTEGER        display interpolation points (default: 1200)\n",
    "  --x-max NUMBER          x-axis maximum (default: 100)\n",
    "  --y-max NUMBER          y-axis maximum (default: auto, at least 105)\n",
    "  --label-series NAME     GDAR, WF3, WF5, or equilibrium (default: GDAR)\n",
    "  --no-labels             suppress scenario labels\n",
    "  --point-x-column NAME   point x column (default: auto)\n",
    "  --gdar-column NAME      GDAR curve/point column (default: auto)\n",
    "  --wf3-column NAME       WF 3-generation curve/point column (default: auto)\n",
    "  --wf5-column NAME       WF 5-generation curve/point column (default: auto)\n",
    "  --equilibrium-column NAME  equilibrium column (default: auto)\n",
    "  --input-scale MODE      auto, fraction, or percent (default: auto)\n",
    "  --help                  show this message\n",
    sep = ""
  )
}

parse_args <- function(x) {
  if (!length(x)) {
    usage()
    stop("--curve-file is required", call. = FALSE)
  }
  out <- list()
  i <- 1L
  flags <- c("include-equilibrium", "include-zero-deme-endpoint", "all-wf-points", "no-labels", "help")
  while (i <= length(x)) {
    token <- x[[i]]
    if (token %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0L)
    }
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token, call. = FALSE)
    }
    key <- sub("^--", "", token)
    if (key %in% flags) {
      out[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(x) || startsWith(x[[i + 1L]], "--")) {
      stop("Option --", key, " requires a value", call. = FALSE)
    }
    out[[key]] <- x[[i + 1L]]
    i <- i + 2L
  }
  out
}

raw_command_args <- commandArgs(trailingOnly = TRUE)
# The BioMA Python wrapper uses a compact positional interface.  Keep it
# alongside the more expressive named interface so the script remains useful
# as a standalone command.
legacy_positional <- length(raw_command_args) > 0L && !startsWith(raw_command_args[[1L]], "--")
if (legacy_positional) {
  if (length(raw_command_args) < 3L) {
    usage()
    stop("Positional mode requires curve_summary, scenario_points, and outdir", call. = FALSE)
  }
  args <- list(
    `curve-file` = raw_command_args[[1L]],
    `point-file` = raw_command_args[[2L]],
    outdir = raw_command_args[[3L]],
    species = if (length(raw_command_args) >= 4L) raw_command_args[[4L]] else "",
    title = if (length(raw_command_args) >= 5L) raw_command_args[[5L]] else "",
    `plot-equilibrium` = if (length(raw_command_args) >= 6L) raw_command_args[[6L]] else "false",
    `plot-raw` = if (length(raw_command_args) >= 7L) raw_command_args[[7L]] else "false"
  )
} else {
  args <- parse_args(raw_command_args)
}
get_arg <- function(name, default = NULL) {
  if (!is.null(args[[name]])) as.character(args[[name]]) else default
}
has_flag <- function(name) isTRUE(args[[name]])

curve_file <- get_arg("curve-file")
if (is.null(curve_file) || !nzchar(curve_file)) {
  usage()
  stop("--curve-file is required", call. = FALSE)
}
point_file <- get_arg("point-file")
short_file <- get_arg("short-file")
meta_file <- get_arg("meta-file")
label_file <- get_arg("scenario-label-file")
outdir <- get_arg("outdir", dirname(normalizePath(curve_file, mustWork = FALSE)))
prefix <- get_arg("prefix", "wfmoment_2d")
if (legacy_positional) {
  prefix <- "Figure_wfmoment_2D_deme"
}
as_bool <- function(value, default = FALSE) {
  if (is.null(value) || !length(value) || !nzchar(as.character(value))) return(default)
  tolower(as.character(value)) %in% c("1", "true", "yes", "y", "on")
}
include_equilibrium <- if (legacy_positional) as_bool(get_arg("plot-equilibrium"), FALSE) else has_flag("include-equilibrium")
plot_raw <- if (legacy_positional) as_bool(get_arg("plot-raw"), FALSE) else as_bool(get_arg("plot-raw", "false"), FALSE)
exclude_zero_deme_endpoint <- !has_flag("include-zero-deme-endpoint")
all_wf_points <- has_flag("all-wf-points")
no_labels <- has_flag("no-labels")
tail_mode <- tolower(get_arg("tail-mode", "auto"))
tail_cutoff_demes <- suppressWarnings(as.numeric(get_arg("tail-cutoff-demes", "5")))
min_display_demes <- suppressWarnings(as.numeric(get_arg("min-display-demes", "5")))
tail_end_arg <- get_arg("tail-end", "auto")
tail_end_auto <- tolower(tail_end_arg) %in% c("auto", "x-max")
tail_exponent <- suppressWarnings(as.numeric(get_arg("tail-exponent", "1.5")))
gdar_display_epsilon <- suppressWarnings(as.numeric(get_arg("gdar-display-epsilon", "0.001")))
display_monotone <- tolower(get_arg("display-monotone", "auto"))
line_method_arg <- tolower(get_arg("line-method", "auto"))
line_method <- if (line_method_arg == "auto") if (plot_raw) "raw" else "spline" else line_method_arg
n_grid <- suppressWarnings(as.integer(get_arg("n-grid", "1200")))
x_max <- suppressWarnings(as.numeric(get_arg("x-max", "100")))
y_max_arg <- get_arg("y-max", "auto")
label_series <- toupper(get_arg("label-series", "GDAR"))
plot_title <- get_arg("title", "")
plot_species <- get_arg("species", "")
# BioMA compute writes metadata.tsv beside curve_summary.tsv.  Discover it by
# default so the positional wrapper can reconstruct the GDAR theory line without
# adding another positional argument.
if (is.null(meta_file)) {
  metadata_candidate <- file.path(dirname(curve_file), "metadata.tsv")
  if (file.exists(metadata_candidate)) meta_file <- metadata_candidate
}
point_x_column <- get_arg("point-x-column", "auto")
gdar_column_arg <- get_arg("gdar-column", "auto")
wf3_column_arg <- get_arg("wf3-column", "auto")
wf5_column_arg <- get_arg("wf5-column", "auto")
eq_column_arg <- get_arg("equilibrium-column", "auto")
input_scale <- tolower(get_arg("input-scale", "auto"))

if (!file.exists(curve_file)) stop("curve file does not exist: ", curve_file, call. = FALSE)
if (!is.null(point_file) && !file.exists(point_file)) stop("point file does not exist: ", point_file, call. = FALSE)
if (!is.null(short_file) && !file.exists(short_file)) stop("short file does not exist: ", short_file, call. = FALSE)
if (!is.null(meta_file) && !file.exists(meta_file)) stop("meta file does not exist: ", meta_file, call. = FALSE)
if (!is.null(label_file) && !file.exists(label_file)) stop("scenario label file does not exist: ", label_file, call. = FALSE)
if (!tail_mode %in% c("auto", "none", "smooth", "zero", "linear")) stop("--tail-mode must be auto, none, smooth, zero, or linear", call. = FALSE)
if (!is.finite(tail_cutoff_demes) || tail_cutoff_demes < 1) stop("--tail-cutoff-demes must be at least 1", call. = FALSE)
if (!is.finite(min_display_demes) || min_display_demes < 1) stop("--min-display-demes must be at least 1", call. = FALSE)
if (tail_end_auto) {
  # The display-only bridge terminates on the requested 100% boundary.  The
  # complete-extinction row is still retained in the raw export, but is not
  # used as an observed WF state.
  tail_end <- x_max
} else {
  tail_end <- suppressWarnings(as.numeric(tail_end_arg))
}
if (!is.finite(tail_end) || tail_end < x_max) stop("--tail-end must be numeric and >= --x-max", call. = FALSE)
if (!is.finite(tail_exponent) || tail_exponent <= 0) stop("--tail-exponent must be positive", call. = FALSE)
if (!is.finite(gdar_display_epsilon) || gdar_display_epsilon < 0 || gdar_display_epsilon >= 1) {
  stop("--gdar-display-epsilon must be in [0, 1)", call. = FALSE)
}
if (!display_monotone %in% c("auto", "none", "cummin")) stop("--display-monotone must be auto, none, or cummin", call. = FALSE)
if (!line_method %in% c("raw", "linear", "spline")) stop("--line-method must be raw, linear, or spline", call. = FALSE)
if (!input_scale %in% c("auto", "fraction", "percent")) stop("--input-scale must be auto, fraction, or percent", call. = FALSE)
if (!is.finite(n_grid) || n_grid < 20L) stop("--n-grid must be at least 20", call. = FALSE)
if (!is.finite(x_max) || x_max <= 0) stop("--x-max must be positive", call. = FALSE)
if (!label_series %in% c("GDAR", "WF3", "WF5", "EQUILIBRIUM")) {
  stop("--label-series must be GDAR, WF3, WF5, or equilibrium", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_table_auto <- function(path) {
  con1 <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  first <- tryCatch(readLines(con1, n = 1L, warn = FALSE, encoding = "UTF-8"),
                    finally = close(con1))
  if (!length(first)) stop("empty table: ", path, call. = FALSE)
  sep <- if (grepl("\\t", first[[1]])) "\t" else if (grepl(",", first[[1]])) "," else ""
  # Re-open because the first line was consumed while detecting the delimiter.
  con2 <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
  tryCatch(
    utils::read.table(con2, header = TRUE, sep = sep, quote = "\"", comment.char = "",
                      check.names = FALSE, na.strings = c("", "NA", "NaN", "NULL", "."),
                      stringsAsFactors = FALSE),
    finally = close(con2)
  )
}

as_numeric <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  z <- trimws(as.character(x))
  z[z %in% c("", "NA", "NaN", "NULL", ".", "-")] <- NA_character_
  suppressWarnings(as.numeric(gsub(",", "", z, fixed = TRUE)))
}

normalise_name <- function(x) tolower(gsub("[^a-z0-9]", "", x))

find_column <- function(df, candidates, explicit = NULL, required = FALSE) {
  nms <- names(df)
  if (!is.null(explicit) && nzchar(explicit) && tolower(explicit) != "auto") {
    if (explicit %in% nms) return(explicit)
    j <- which(normalise_name(nms) == normalise_name(explicit))
    if (length(j)) return(nms[[j[[1]]]])
    if (required) stop("Column '", explicit, "' not found. Available columns: ", paste(nms, collapse = ", "), call. = FALSE)
    return(NULL)
  }
  if (!length(candidates)) {
    if (required) stop("No candidate column names supplied", call. = FALSE)
    return(NULL)
  }
  j <- match(normalise_name(candidates), normalise_name(nms), nomatch = 0L)
  if (any(j > 0L)) return(nms[[j[[which(j > 0L)[[1]]]]]])
  if (required) stop("None of these columns found: ", paste(candidates, collapse = ", "),
                     ". Available columns: ", paste(nms, collapse = ", "), call. = FALSE)
  NULL
}

finite_range <- function(x) {
  x <- as_numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) c(NA_real_, NA_real_) else range(x)
}

normalise_value_scale <- function(x) {
  z <- tolower(trimws(as.character(x)))
  if (any(z %in% c("percent", "pct", "percentage", "%"), na.rm = TRUE)) return("percent")
  if (any(z %in% c("fraction", "frac", "proportion", "prop"), na.rm = TRUE)) return("fraction")
  NULL
}

infer_value_scale <- function(x, column_name = "", hint = NULL) {
  if (!is.null(hint) && length(hint) && hint[[1L]] %in% c("fraction", "percent")) return(hint[[1L]])
  if (input_scale %in% c("fraction", "percent")) return(input_scale)
  finite <- as_numeric(x)
  finite <- finite[is.finite(finite)]
  if (!length(finite)) return("fraction")
  # Explicit percent-like names take precedence for values in [0, 1].
  is_pct_name <- grepl("pct|percent|percentage", tolower(column_name))
  if (is_pct_name || any(abs(finite) > 1.5)) "percent" else "fraction"
}

value_to_percent <- function(x, scale = input_scale, column_name = "") {
  z <- as_numeric(x)
  resolved <- infer_value_scale(z, column_name = column_name, hint = scale)
  if (resolved == "fraction") z * 100 else z
}

loss_to_percent <- function(df, explicit = "auto") {
  xcol <- find_column(df,
                      c("habitat_loss_pct", "habitat_loss_percent", "habitat_loss",
                        "loss_pct", "loss_percent", "x"),
                      explicit = if (tolower(explicit) == "auto") NULL else explicit)
  if (!is.null(xcol)) {
    x <- as_numeric(df[[xcol]])
    if (grepl("pct|percent", tolower(xcol)) || (length(x[is.finite(x)]) && any(abs(x[is.finite(x)]) > 1.5))) {
      return(list(x = x, column = xcol, source = "loss"))
    }
    return(list(x = x * 100, column = xcol, source = "loss_fraction"))
  }
  acol <- find_column(df, c("A_remaining", "A_remain", "frac_remain", "fraction_remaining", "area_remaining"))
  if (!is.null(acol)) {
    a <- as_numeric(df[[acol]])
    return(list(x = (1 - a) * 100, column = acol, source = "A_remaining"))
  }
  stop("Could not identify habitat-loss column. Expected habitat_loss_pct or A_remaining.", call. = FALSE)
}

scenario_id <- function(df) {
  c <- find_column(df, c("scenario", "scenario_raw", "scenario_id", "future_scenario"))
  if (is.null(c)) return(rep("scenario", nrow(df)))
  as.character(df[[c]])
}

pretty_scenario <- function(x) {
  x <- as.character(x)
  out <- x
  for (i in seq_along(x)) {
    z <- x[[i]]
    m <- regexec("^([0-9]{4})[_-]([0-9]{4})[_-](?:ssp)?([0-9]+)$", tolower(z))
    hit <- regmatches(tolower(z), m)[[1]]
    if (length(hit) == 4L) {
      out[[i]] <- paste0(hit[[2]], "-", hit[[3]], "\nSSP", hit[[4]])
    } else {
      out[[i]] <- gsub("[_-]+", " ", z)
    }
  }
  out
}

scenario_period_class <- function(x) {
  z <- tolower(as.character(x))
  out <- rep(NA_character_, length(z))
  out[grepl("2061[^0-9]*2080", z, perl = TRUE)] <- "early"
  out[grepl("2081[^0-9]*2100", z, perl = TRUE)] <- "late"
  out
}

read_labels <- function(path) {
  if (is.null(path)) return(setNames(character(), character()))
  d <- read_table_auto(path)
  if (ncol(d) < 2L) stop("scenario label file must have at least two columns", call. = FALSE)
  setNames(as.character(d[[2L]]), as.character(d[[1L]]))
}

label_map <- read_labels(label_file)
apply_labels <- function(ids, df = NULL) {
  if (!is.null(df)) {
    c <- find_column(df, c("scenario_label", "label", "scenario_name"))
    if (!is.null(c)) {
      vals <- as.character(df[[c]])
      vals[!is.na(vals) & nzchar(vals)] <- gsub("[_-]", "-", vals[!is.na(vals) & nzchar(vals)])
      return(ifelse(!is.na(vals) & nzchar(vals), vals, pretty_scenario(ids)))
    }
  }
  mapped <- unname(label_map[ids])
  mapped[is.na(mapped) | !nzchar(mapped)] <- pretty_scenario(ids)[is.na(mapped) | !nzchar(mapped)]
  mapped
}

curve_raw <- read_table_auto(curve_file)
if (!nrow(curve_raw)) stop("curve file has no rows: ", curve_file, call. = FALSE)

# Locate the x axis before determining whether the curve is a wide or long table.
curve_x <- loss_to_percent(curve_raw)
curve_raw$.wfmoment_x_pct <- curve_x$x
curve_raw$.wfmoment_scenario <- scenario_id(curve_raw)
curve_occupied_col <- find_column(curve_raw, c("occupied_demes", "n_occupied_demes", "num_demes", "n_demes"))
curve_mode_col <- find_column(curve_raw, c("loss_mode", "mode", "loss_type"))
curve_mode_values <- if (is.null(curve_mode_col)) character() else tolower(trimws(as.character(curve_raw[[curve_mode_col]])))
metadata_mode <- character()
if (!is.null(meta_file) && file.exists(meta_file)) {
  meta_probe <- tryCatch(read_table_auto(meta_file), error = function(e) NULL)
  if (!is.null(meta_probe)) {
    mcol <- find_column(meta_probe, c("loss_mode", "mode", "loss_type"))
    if (!is.null(mcol)) metadata_mode <- tolower(trimws(as.character(meta_probe[[mcol]])))
  }
}
edge_mode <- any(c(curve_mode_values, metadata_mode) == "edge", na.rm = TRUE)
complete_loss_endpoint <- FALSE
if (is.finite(x_max)) {
  endpoint_hit <- is.finite(curve_x$x) & curve_x$x >= 100 - 1e-8
  if (any(endpoint_hit)) {
    if (is.null(curve_occupied_col)) {
      # Tables without an occupied-deme column are accepted when they expose a
      # finite x=100 endpoint; the value itself is retained as raw input.
      complete_loss_endpoint <- TRUE
    } else {
      occupied_probe <- as_numeric(curve_raw[[curve_occupied_col]])
      complete_loss_endpoint <- any(endpoint_hit & is.finite(occupied_probe) & occupied_probe <= 0)
    }
  }
}
if (tail_end_auto && complete_loss_endpoint) {
  # Treat 100% as the biological endpoint even when the visible x-axis is a
  # zoomed subset.  make_display_line clips the synthetic tail to x_max.
  tail_end <- 100
}
monotone_wf <- if (display_monotone == "cummin") TRUE else if (display_monotone == "none") FALSE else edge_mode
# A complete-extinction row is a sentinel, not a regular point on a continuous
# trajectory.  For edge contraction with a complete-loss endpoint, the default
# display therefore uses a display-only tail that reaches the x-axis boundary
# continuously.
# Raw tables remain untouched; users who want a strict truncation can request
# --tail-mode none, while --include-zero-deme-endpoint is retained for QA.
effective_tail_mode <- if (tail_mode == "auto") {
  if (edge_mode && complete_loss_endpoint) "smooth" else "none"
} else tail_mode
if (has_flag("include-zero-deme-endpoint") && tail_mode == "auto") effective_tail_mode <- "none"

curve_type_col <- find_column(curve_raw, c("curve_type", "curve", "time_group", "metric"))
curve_value_candidates <- c("mean_pi_remaining", "pi_remaining", "value", "y", "predicted_pi_remaining")

series_specs <- list(
  GDAR = list(label = "GDAR short-term", color = "#4F5BD5", linetype = "solid", shape = 16L),
  WF3 = list(label = "Wfmoment 3 generation", color = "#D9500B", linetype = "dotted", shape = 15L),
  WF5 = list(label = "Wfmoment 5 generation", color = "#31B95B", linetype = "solid", shape = 17L),
  EQUILIBRIUM = list(label = "Wfmoment equilibrium", color = "#8064A2", linetype = "dotdash", shape = 18L)
)

pick_long_series <- function(df, type_col, patterns) {
  if (is.null(type_col)) return(NULL)
  vals <- tolower(as.character(df[[type_col]]))
  hit <- rep(FALSE, length(vals))
  for (p in patterns) hit <- hit | grepl(p, vals, perl = TRUE)
  if (!any(hit)) return(NULL)
  vcol <- find_column(df, curve_value_candidates)
  if (is.null(vcol)) return(NULL)
  occ_col <- find_column(df, c("occupied_demes", "n_occupied_demes", "num_demes", "n_demes"))
  scale_col <- find_column(df, c("value_scale", "scale", "unit", "units"))
  scale_hint <- if (is.null(scale_col)) NULL else normalise_value_scale(df[[scale_col]][hit])
  list(x = df$.wfmoment_x_pct[hit], raw = as_numeric(df[[vcol]][hit]),
       occupied = if (is.null(occ_col)) NULL else as_numeric(df[[occ_col]][hit]),
       source_column = vcol, source_type = type_col, scale = scale_hint,
       index = which(hit))
}

extract_wide_series <- function(df, explicit, candidates) {
  c <- find_column(df, candidates, explicit = explicit)
  if (is.null(c)) return(NULL)
  occ_col <- find_column(df, c("occupied_demes", "n_occupied_demes", "num_demes", "n_demes"))
  scale_col <- find_column(df, c("value_scale", "scale", "unit", "units"))
  scale_hint <- if (is.null(scale_col)) NULL else normalise_value_scale(df[[scale_col]])
  list(x = df$.wfmoment_x_pct, raw = as_numeric(df[[c]]), source_column = c,
       occupied = if (is.null(occ_col)) NULL else as_numeric(df[[occ_col]]),
       source_type = "wide", scale = scale_hint)
}

curve_series <- list()

# GDAR: prefer a directly reported gdar_remaining column.  A short-term column
# is accepted for legacy tables; if neither exists, the formula is reconstructed
# from A_remaining and z_gdar below.
gdar_curve_obj <- extract_wide_series(
  curve_raw, gdar_column_arg,
  c("gdar_remaining", "gdar_remaining_pct", "GDAR_short_term", "gdar_short_term",
    "pi_remaining_short_from_GDAR", "pi_remaining_short", "pi_remaining_short_lookup")
)
if (is.null(gdar_curve_obj)) {
  # WF_Immediate is a model output, not the GDAR theory line.  Only explicit
  # GDAR/short-term labels are accepted here to avoid silently mislabelling it.
  gdar_curve_obj <- pick_long_series(curve_raw, curve_type_col, c("gdar", "short.*term"))
}

wf3_curve_obj <- extract_wide_series(
  curve_raw, wf3_column_arg,
  c("pi_remaining_3gen", "pi_remaining_3_gen", "WF_3_gen", "WF_3gen", "pi_remaining_mid_3gen")
)
if (is.null(wf3_curve_obj)) {
  wf3_curve_obj <- pick_long_series(curve_raw, curve_type_col, c("3[ _-]*gen", "gen[ _-]*3", "3generation"))
}

wf5_curve_obj <- extract_wide_series(
  curve_raw, wf5_column_arg,
  c("pi_remaining_5gen", "pi_remaining_5_gen", "WF_5_gen", "WF_5gen", "pi_remaining_mid_5gen")
)
if (is.null(wf5_curve_obj)) {
  wf5_curve_obj <- pick_long_series(curve_raw, curve_type_col, c("5[ _-]*gen", "gen[ _-]*5", "5generation"))
}

eq_curve_obj <- extract_wide_series(
  curve_raw, eq_column_arg,
  c("pi_remaining_equilibrium", "WF_equilibrium", "equilibrium", "pi_remaining_eq")
)
if (is.null(eq_curve_obj)) {
  eq_curve_obj <- pick_long_series(curve_raw, curve_type_col, c("equilibrium", "long[ _-]*term", "eq"))
}

read_scalar_z <- function(df, species = NULL) {
  if (!is.null(species) && nzchar(species) && "species" %in% names(df)) {
    selected <- as.character(df$species) == species
    if (any(selected, na.rm = TRUE)) df <- df[selected, , drop = FALSE]
  }
  c <- find_column(df, c("z_gdar", "zGDAR", "gdar_z"))
  if (is.null(c)) return(NA_real_)
  z <- as_numeric(df[[c]])
  z <- z[is.finite(z)]
  if (length(z)) z[[1L]] else NA_real_
}

z_gdar <- read_scalar_z(curve_raw, plot_species)
if (!is.finite(z_gdar) && !is.null(short_file)) z_gdar <- read_scalar_z(read_table_auto(short_file), plot_species)
if (!is.finite(z_gdar) && !is.null(meta_file)) z_gdar <- read_scalar_z(read_table_auto(meta_file), plot_species)

if (is.null(gdar_curve_obj) && is.finite(z_gdar)) {
  # Formula-derived GDAR is explicitly marked as derived in the export table.
  gdar_curve_obj <- list(x = curve_raw$.wfmoment_x_pct,
                         raw = (pmax(0, 1 - curve_raw$.wfmoment_x_pct / 100) ^ z_gdar),
                          occupied = if (is.null(curve_occupied_col)) NULL else as_numeric(curve_raw[[curve_occupied_col]]),
                          source_column = "derived:A_remaining^z_gdar",
                          source_type = "derived", scale = "fraction")
}
if (is.null(wf3_curve_obj)) warning("No WF 3-generation curve found; it will be omitted.")
if (is.null(wf5_curve_obj)) warning("No WF 5-generation curve found; it will be omitted.")
if (is.null(gdar_curve_obj)) warning("No GDAR curve found; GDAR line will be omitted.")
if (include_equilibrium && is.null(eq_curve_obj)) warning("Equilibrium requested but no equilibrium curve found; it will be omitted.")

collapse_curve <- function(obj, series_name) {
  if (is.null(obj)) return(NULL)
  raw <- as_numeric(obj$raw)
  x <- as_numeric(obj$x)
  occupied <- obj$occupied
  if (is.null(occupied) || length(occupied) != length(x)) occupied <- rep(NA_real_, length(x))
  raw_df_all <- data.frame(habitat_loss_pct = x, value_raw = raw,
                           occupied_demes = as_numeric(occupied),
                           series = series_name, source_column = obj$source_column,
                            stringsAsFactors = FALSE)
  # Resolve the unit once for the whole series.  Inferring it separately for
  # each x-group misclassifies small percent values near extinction as
  # fractions (for example 1.46% becomes 146%).
  series_scale <- infer_value_scale(raw, column_name = obj$source_column,
                                    hint = obj$scale %||% NULL)
  raw_df_all$value_pct <- value_to_percent(raw_df_all$value_raw, scale = series_scale,
                                           column_name = obj$source_column)
  keep_x <- is.finite(x)
  if (!any(keep_x & is.finite(raw))) return(NULL)

  # Collapse duplicate x values for display while retaining all original rows
  # (including NA values) in raw_df_all. A wholly missing x-group remains NA,
  # allowing the display builder to break rather than interpolate across it.
  ord <- order(x, seq_along(x), na.last = TRUE)
  x_ord <- x[ord]; raw_ord <- raw[ord]; occ_ord <- occupied[ord]
  idx_by_x <- split(which(is.finite(x_ord)), x_ord[is.finite(x_ord)], drop = TRUE)
  x_unique <- as.numeric(names(idx_by_x))
  y_unique <- vapply(idx_by_x, function(ii) {
    vals <- raw_ord[ii]; vals <- vals[is.finite(vals)]
    if (!length(vals)) NA_real_ else mean(value_to_percent(vals, scale = series_scale,
                                                            column_name = obj$source_column))
  }, numeric(1))
  occ_unique <- vapply(idx_by_x, function(ii) {
    vals <- occ_ord[ii]; vals <- vals[is.finite(vals)]
    if (!length(vals)) NA_real_ else max(vals)
  }, numeric(1))
  display <- data.frame(habitat_loss_pct = x_unique, value_pct = y_unique,
                        occupied_demes = occ_unique, series = series_name,
                        display_source = "raw", stringsAsFactors = FALSE)
  display <- display[order(display$habitat_loss_pct), , drop = FALSE]
  # The all-zero row represents complete extinction, not a regular deme state:
  # with no extant deme, species-wide pi is not a regular continuous state.
  # Keep the sentinel in the raw export, but exclude it from the default
  # display.  A theoretical GDAR line is rebuilt below from z_gdar.
  if (exclude_zero_deme_endpoint && nrow(display) >= 2L) {
    terminal_zero <- is.finite(display$occupied_demes) & display$occupied_demes <= 0
    candidate <- !terminal_zero
    if (sum(candidate & is.finite(display$value_pct)) >= 2L) display <- display[candidate, , drop = FALSE]
  }
  # A one-deme value remains a valid raw diagnostic, but joining it to the
  # multi-deme trajectory creates a near-vertical segment at the right edge.
  # WF publication lines therefore stop at the last state that still has the
  # configured minimum number of demes. GDAR is area-theoretical and is handled
  # separately, so it must not inherit this deme cutoff.
  if (series_name != "GDAR" && nrow(display) >= 2L && any(is.finite(display$occupied_demes))) {
    candidate <- !is.finite(display$occupied_demes) | display$occupied_demes >= min_display_demes
    # The explicit diagnostic flag deliberately opts back into the zero-deme
    # sentinel; it is never included by the default display.
    if (!exclude_zero_deme_endpoint) {
      candidate <- candidate | (is.finite(display$occupied_demes) & display$occupied_demes <= 0)
    }
    if (sum(candidate & is.finite(display$value_pct)) >= 2L) display <- display[candidate, , drop = FALSE]
  }
  if (monotone_wf && series_name %in% c("WF3", "WF5") && nrow(display) >= 2L) {
    # Display-only envelope for the connected edge-contraction presentation.
    # Once demes become disconnected, species-wide pi can increase because
    # between-deme divergence grows. Preserve every raw value, but use a
    # bounded cumulative minimum for this publication-style line so that the
    # edge trajectory remains a non-increasing area-loss reference. Random or
    # fragmentation runs can disable this with --display-monotone none.
    bounded <- pmin(100, pmax(0, display$value_pct))
    finite_y <- is.finite(bounded)
    if (any(finite_y)) {
      runs_y <- rle(finite_y)
      ends_y <- cumsum(runs_y$lengths)
      starts_y <- c(1L, head(ends_y, -1L) + 1L)
      for (run_idx in which(runs_y$values)) {
        ii <- starts_y[[run_idx]]:ends_y[[run_idx]]
        bounded[ii] <- cummin(bounded[ii])
      }
    }
    display$value_pct <- bounded
    display$display_source <- "edge_cummin_display"
  }
  # raw_df_all intentionally retains NA rows; they are useful diagnostics for
  # disconnected equilibrium states and are excluded only from the display.
  list(raw = raw_df_all, display = display, source_column = obj$source_column)
}

curve_objs <- list(
  GDAR = collapse_curve(gdar_curve_obj, "GDAR"),
  WF3 = collapse_curve(wf3_curve_obj, "WF3"),
  WF5 = collapse_curve(wf5_curve_obj, "WF5"),
  EQUILIBRIUM = collapse_curve(eq_curve_obj, "EQUILIBRIUM")
)
if (!include_equilibrium) curve_objs$EQUILIBRIUM <- NULL
curve_objs <- curve_objs[!vapply(curve_objs, is.null, logical(1))]
if (!length(curve_objs)) stop("No finite curve series could be extracted", call. = FALSE)

# GDAR is an area-only theory line, not a deme-by-deme WF trajectory.  When
# z_gdar is available, use a dense theoretical display grid and stop at a tiny
# positive remaining-area fraction.  This keeps the line attached to the 100%
# x-axis boundary without treating the A_remaining=0 sentinel as an observed
# pi value.  The original GDAR rows remain in curve_raw_export.
gdar_display_theory <- FALSE
gdar_display_max_x <- NA_real_
if (is.finite(z_gdar) && !is.null(curve_objs$GDAR) && exclude_zero_deme_endpoint) {
  gdar_display_max_x <- min(x_max, 100 * (1 - gdar_display_epsilon))
  if (is.finite(gdar_display_max_x) && gdar_display_max_x > 0) {
    gx <- seq(0, gdar_display_max_x, length.out = max(200L, n_grid))
    curve_objs$GDAR$display <- data.frame(
      habitat_loss_pct = gx,
      value_pct = 100 * pmax(0, 1 - gx / 100) ^ z_gdar,
      occupied_demes = NA_real_,
      series = "GDAR",
      display_source = "theory_A_remaining^z_gdar",
      stringsAsFactors = FALSE
    )
    gdar_display_theory <- TRUE
  }
}

append_tail <- function(display, mode, x_end, cutoff_demes = tail_cutoff_demes,
                        tail_anchor = tail_end, exponent = tail_exponent) {
  if (mode == "none" || !nrow(display)) return(display)
  display <- display[order(display$habitat_loss_pct), , drop = FALSE]
  series_name <- unique(display$series)[[1L]]
  # GDAR has its own epsilon-truncated theoretical grid; never append a WF
  # extinction tail to it.
  if (identical(series_name, "GDAR")) return(display)
  # For edge contraction, start the display-only tail at the last state that
  # still has a small connected core rather than at the one-deme endpoint.
  if ("occupied_demes" %in% names(display) && is.finite(cutoff_demes)) {
    eligible <- is.finite(display$occupied_demes) & display$occupied_demes >= cutoff_demes &
      is.finite(display$habitat_loss_pct) & is.finite(display$value_pct)
    if (any(eligible)) {
      cutoff_x <- max(display$habitat_loss_pct[eligible])
      display <- display[display$habitat_loss_pct <= cutoff_x, , drop = FALSE]
    }
  }
  finite_idx <- which(is.finite(display$habitat_loss_pct) & is.finite(display$value_pct))
  if (!length(finite_idx)) return(display)
  last_idx <- tail(finite_idx, 1L)
  last_x <- display$habitat_loss_pct[[last_idx]]
  last_y <- display$value_pct[[last_idx]]
  if (!is.finite(last_x)) return(display)
  if (mode == "smooth") {
    # Use a bounded, monotone display-only bridge to the complete-loss
    # boundary.  1 - u^q is deliberately used instead of joining the raw
    # one/zero-deme sentinels: it follows the reference shape (a gentle start
    # followed by accelerated loss), reaches exactly zero at x_end, and never
    # changes the raw calculation tables.
    if (!is.finite(tail_anchor) || tail_anchor <= last_x) return(display)
    span <- max(0.01, (tail_anchor - last_x) / max(tail_anchor, 1))
    n_tail <- max(80L, min(600L, as.integer(round(n_grid * min(1, span)))))
    tx <- sort(unique(c(seq(last_x, tail_anchor, length.out = n_tail), x_end)))
    u <- pmin(1, pmax(0, (tx - last_x) / (tail_anchor - last_x)))
    ty <- pmax(0, last_y) * pmax(0, 1 - u ^ exponent)
    add <- data.frame(habitat_loss_pct = tx[-1L], value_pct = ty[-1L],
                      occupied_demes = NA_real_, series = series_name,
                      display_source = "smooth_terminal_tail", stringsAsFactors = FALSE)
    return(rbind(display, add))
  }
  if (!is.finite(x_end) || last_x >= x_end) return(display)
  if (mode == "zero") {
    add <- data.frame(habitat_loss_pct = x_end, value_pct = 0,
                      occupied_demes = NA_real_, series = series_name,
                      display_source = "appended_zero_tail", stringsAsFactors = FALSE)
  } else {
    # A linear continuation to zero is retained as an explicit legacy option.
    add <- data.frame(habitat_loss_pct = x_end, value_pct = 0,
                      occupied_demes = NA_real_, series = series_name,
                      display_source = "appended_linear_zero_tail", stringsAsFactors = FALSE)
  }
  rbind(display, add)
}

make_display_line <- function(display, method, x_end, n_grid, tail_mode) {
  display <- append_tail(display, tail_mode, x_end)
  display <- display[order(display$habitat_loss_pct), , drop = FALSE]
  # Keep the exported display table on the requested axis.  A custom tail
  # endpoint beyond x_max is supported for exploratory plots, but is clipped
  # from the default visible export.
  display <- display[!is.finite(display$habitat_loss_pct) |
                       display$habitat_loss_pct <= x_end + 1e-10, , drop = FALSE]
  if (method == "raw") return(display)
  finite <- is.finite(display$habitat_loss_pct) & is.finite(display$value_pct)
  if (!any(finite)) return(display[FALSE, , drop = FALSE])
  runs <- rle(finite)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  segments <- list()
  for (k in which(runs$values)) {
    idx <- starts[[k]]:ends[[k]]
    seg <- display[idx, , drop = FALSE]
    seg <- seg[!duplicated(seg$habitat_loss_pct), , drop = FALSE]
    x <- seg$habitat_loss_pct
    y <- seg$value_pct
    if (length(x) < 3L) {
      out <- seg[, c("habitat_loss_pct", "value_pct", "series", "display_source"), drop = FALSE]
      out$display_source <- method
    } else {
      grid <- seq(min(x), max(x), length.out = max(20L, n_grid))
      if (method == "spline" && length(x) >= 4L) {
        fn <- tryCatch(stats::splinefun(x, y, method = "hyman"), error = function(e) NULL)
        yy <- if (is.null(fn)) stats::approx(x, y, xout = grid, rule = 1)$y else fn(grid)
      } else {
        yy <- stats::approx(x, y, xout = grid, rule = 1)$y
      }
      has_tail <- any(grepl("smooth_terminal_tail", seg$display_source, fixed = TRUE))
      source_name <- if (has_tail) paste0(method, "+smooth_tail") else method
      out <- data.frame(habitat_loss_pct = grid, value_pct = yy,
                        series = unique(display$series)[[1L]], display_source = source_name,
                        stringsAsFactors = FALSE)
    }
    if (length(segments)) {
      segments[[length(segments) + 1L]] <- data.frame(
        habitat_loss_pct = NA_real_, value_pct = NA_real_,
        series = unique(display$series)[[1L]], display_source = "NA_break",
        stringsAsFactors = FALSE
      )
    }
    segments[[length(segments) + 1L]] <- out
  }
  if (length(segments)) do.call(rbind, segments) else display[FALSE, , drop = FALSE]
}

curve_raw_export <- do.call(rbind, lapply(curve_objs, `[[`, "raw"))
curve_display <- do.call(rbind, lapply(curve_objs, function(obj) {
  make_display_line(obj$display, line_method, x_max, n_grid, effective_tail_mode)
}))

# ------------------------- scenario points -------------------------
point_objs <- list()
if (!is.null(point_file)) {
  points_raw <- read_table_auto(point_file)
  if (!nrow(points_raw)) warning("point file has no rows: ", point_file)
} else {
  points_raw <- NULL
}

short_raw <- if (!is.null(short_file)) read_table_auto(short_file) else NULL

extract_point_series <- function(df, series, explicit = "auto") {
  if (is.null(df) || !nrow(df)) return(NULL)
  xobj <- loss_to_percent(df, explicit = point_x_column)
  df$.wfmoment_x_pct <- xobj$x
  ids <- scenario_id(df)
  labels <- apply_labels(ids, df)
  time_col <- find_column(df, c("time_forward", "generation", "generations"))
  time_vals <- if (is.null(time_col)) rep(NA_real_, nrow(df)) else as_numeric(df[[time_col]])
  type_col <- find_column(df, c("curve_type", "curve", "time_group", "metric"))
  if (series == "GDAR") {
    obj <- extract_wide_series(df, explicit,
      c("GDAR_short_term", "gdar_remaining", "gdar_remaining_pct", "pi_remaining_short",
        "pi_remaining_short_lookup", "pi_remaining_short_from_GDAR"))
    if (is.null(obj) && !is.null(type_col)) obj <- pick_long_series(df, type_col, c("gdar", "short.*term"))
  } else if (series == "WF3") {
    obj <- extract_wide_series(df, explicit,
      c("WF_3_gen", "WF_3gen", "pi_remaining_3gen", "pi_remaining_3_gen", "pi_remaining_3gen_lookup"))
    if (is.null(obj) && !is.null(type_col)) obj <- pick_long_series(df, type_col, c("3[ _-]*gen", "gen[ _-]*3", "3generation"))
    # Some phase-2 tables are long by time_forward rather than curve_type.
    if (is.null(obj)) {
      tcol <- find_column(df, c("time_forward", "generation", "generations"))
      vcol <- find_column(df, c("pi_remaining_mid", "pi_remaining", "mean_pi_remaining", "value"))
      if (!is.null(tcol) && !is.null(vcol)) {
        tt <- as_numeric(df[[tcol]]); hit <- is.finite(tt) & abs(tt - 3) < 1e-8
        if (any(hit)) obj <- list(x = xobj$x[hit], raw = as_numeric(df[[vcol]][hit]), source_column = vcol, source_type = tcol,
                                  ids = ids[hit], labels = labels[hit], time_forward = tt[hit])
      }
    }
  } else if (series == "WF5") {
    obj <- extract_wide_series(df, explicit,
      c("WF_5_gen", "WF_5gen", "pi_remaining_5gen", "pi_remaining_5_gen", "pi_remaining_5gen_lookup"))
    if (is.null(obj) && !is.null(type_col)) obj <- pick_long_series(df, type_col, c("5[ _-]*gen", "gen[ _-]*5", "5generation"))
    if (is.null(obj)) {
      tcol <- find_column(df, c("time_forward", "generation", "generations"))
      vcol <- find_column(df, c("pi_remaining_mid", "pi_remaining", "mean_pi_remaining", "value"))
      if (!is.null(tcol) && !is.null(vcol)) {
        tt <- as_numeric(df[[tcol]]); hit <- is.finite(tt) & abs(tt - 5) < 1e-8
        if (any(hit)) obj <- list(x = xobj$x[hit], raw = as_numeric(df[[vcol]][hit]), source_column = vcol, source_type = tcol,
                                  ids = ids[hit], labels = labels[hit], time_forward = tt[hit])
      }
    }
  } else {
    obj <- extract_wide_series(df, explicit,
      c("WF_equilibrium", "pi_remaining_equilibrium", "equilibrium", "pi_remaining_eq"))
    if (is.null(obj) && !is.null(type_col)) obj <- pick_long_series(df, type_col, c("equilibrium", "long[ _-]*term", "eq"))
  }
  if (is.null(obj)) return(NULL)
  # extract_wide_series does not carry IDs; attach the full-row metadata here.
  if (is.null(obj$ids)) {
    idx <- obj$index %||% seq_len(nrow(df))
    obj$ids <- ids[idx]
    obj$labels <- labels[idx]
    obj$x <- xobj$x[idx]
    obj$time_forward <- time_vals[idx]
  }
  raw <- as_numeric(obj$raw)
  if (is.null(obj$time_forward)) obj$time_forward <- rep(NA_real_, length(raw))
  n <- min(length(raw), length(obj$x), length(obj$ids), length(obj$labels), length(obj$time_forward))
  point_scale <- infer_value_scale(raw, column_name = obj$source_column,
                                   hint = obj$scale %||% NULL)
  data.frame(scenario = as.character(obj$ids[seq_len(n)]),
             scenario_label = as.character(obj$labels[seq_len(n)]),
             habitat_loss_pct = as_numeric(obj$x[seq_len(n)]),
             time_forward = as_numeric(obj$time_forward[seq_len(n)]),
             value_raw = raw[seq_len(n)],
             value_pct = value_to_percent(raw[seq_len(n)], scale = point_scale,
                                          column_name = obj$source_column),
             series = series, source_column = obj$source_column,
             stringsAsFactors = FALSE)
}

if (!is.null(points_raw)) {
  point_objs$GDAR <- extract_point_series(points_raw, "GDAR", explicit = gdar_column_arg)
  point_objs$WF3 <- extract_point_series(points_raw, "WF3", explicit = wf3_column_arg)
  point_objs$WF5 <- extract_point_series(points_raw, "WF5", explicit = wf5_column_arg)
  point_objs$EQUILIBRIUM <- extract_point_series(points_raw, "EQUILIBRIUM", explicit = eq_column_arg)
}

# Legacy short_file can supply GDAR points when lookup points only contain WF columns.
if (is.null(point_objs$GDAR) && !is.null(short_raw)) {
  point_objs$GDAR <- extract_point_series(short_raw, "GDAR", explicit = gdar_column_arg)
}

# If no direct GDAR point values exist, use the same z_gdar formula as the line.
if (is.null(point_objs$GDAR) && !is.null(points_raw) && is.finite(z_gdar)) {
  xobj <- loss_to_percent(points_raw, explicit = point_x_column)
  ids <- scenario_id(points_raw)
  point_objs$GDAR <- data.frame(
    scenario = ids, scenario_label = apply_labels(ids, points_raw),
    habitat_loss_pct = xobj$x,
    time_forward = {
      tc <- find_column(points_raw, c("time_forward", "generation", "generations"))
      if (is.null(tc)) rep(NA_real_, nrow(points_raw)) else as_numeric(points_raw[[tc]])
    },
    value_raw = pmax(0, 1 - xobj$x / 100) ^ z_gdar,
    value_pct = value_to_percent(pmax(0, 1 - xobj$x / 100) ^ z_gdar, scale = "fraction", column_name = "derived"),
    series = "GDAR", source_column = "derived:A_remaining^z_gdar", stringsAsFactors = FALSE
  )
}

if (!include_equilibrium) point_objs$EQUILIBRIUM <- NULL
point_objs <- point_objs[!vapply(point_objs, is.null, logical(1))]
point_objs <- point_objs[vapply(point_objs, nrow, integer(1)) > 0L]
point_export <- if (length(point_objs)) do.call(rbind, point_objs) else data.frame()
if (nrow(point_export)) {
  point_export$plot_include <- TRUE
  if (!all_wf_points) {
    period <- scenario_period_class(point_export$scenario)
    # The reference figure places the short/mid-term marker at the matching
    # future period: WF3 for 2061-2080 and WF5 for 2081-2100.  For tables with
    # no parseable period, use time_forward when two distinct times are given;
    # otherwise retain the point rather than guessing.
    wf3 <- point_export$series == "WF3"
    wf5 <- point_export$series == "WF5"
    known_period <- !is.na(period)
    point_export$plot_include[wf3 & known_period] <- period[wf3 & known_period] == "early"
    point_export$plot_include[wf5 & known_period] <- period[wf5 & known_period] == "late"
    unknown <- !(known_period) & (wf3 | wf5) & is.finite(point_export$time_forward)
    wf_times <- sort(unique(point_export$time_forward[(wf3 | wf5) & is.finite(point_export$time_forward)]))
    if (length(wf_times) >= 2L) {
      early_time <- wf_times[[1L]]
      late_time <- tail(wf_times, 1L)
      point_export$plot_include[unknown & wf3] <- abs(point_export$time_forward[unknown & wf3] - early_time) < 1e-8
      point_export$plot_include[unknown & wf5] <- abs(point_export$time_forward[unknown & wf5] - late_time) < 1e-8
    }
  }
  # Keep NA rows in the exported point table; create a finite-only view for
  # plotting below. Missing equilibrium values are expected after disconnects.
  point_plot <- point_export[point_export$plot_include & is.finite(point_export$habitat_loss_pct) & is.finite(point_export$value_raw), , drop = FALSE]
} else {
  point_plot <- point_export
}

# Ensure plot only includes series with both a finite line and finite points.
series_present <- unique(curve_display$series)
if (nrow(point_export)) {
  point_export <- point_export[point_export$series %in% series_present, , drop = FALSE]
  point_plot <- point_plot[point_plot$series %in% series_present, , drop = FALSE]
}

# ------------------------- output tables -------------------------
curve_export_path <- file.path(outdir, paste0(prefix, "_curve_data.tsv"))
display_export_path <- file.path(outdir, paste0(prefix, "_curve_display.tsv"))
point_export_path <- file.path(outdir, paste0(prefix, "_point_data.tsv"))
diagnostic_path <- file.path(outdir, paste0(prefix, "_diagnostics.tsv"))
utils::write.table(curve_raw_export, curve_export_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
utils::write.table(curve_display, display_export_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
if (nrow(point_export)) {
  # Scenario labels intentionally contain line breaks for the figure. Quote
  # character fields here so the exported TSV remains parseable.
  utils::write.table(point_export, point_export_path, sep = "\t", row.names = FALSE, quote = TRUE, na = "")
} else {
  utils::write.table(data.frame(), point_export_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

diag_rows <- data.frame(
  item = c("curve_file", "point_file", "short_file", "meta_file", "line_method", "tail_mode",
           "effective_tail_mode", "tail_end_auto", "complete_loss_endpoint", "tail_cutoff_demes", "min_display_demes", "tail_end", "tail_exponent",
           "gdar_display_epsilon", "gdar_display_theory", "gdar_display_max_x", "display_monotone", "edge_mode",
           "exclude_zero_deme_endpoint", "all_wf_points", "n_curve_raw", "n_curve_display", "n_point_raw",
           "n_point_plotted", "z_gdar", "include_equilibrium", "x_max"),
  value = c(curve_file, point_file %||% "", short_file %||% "", meta_file %||% "", line_method, tail_mode,
            effective_tail_mode, as.character(tail_end_auto), as.character(complete_loss_endpoint),
            format(tail_cutoff_demes, digits = 12), format(min_display_demes, digits = 12),
            format(tail_end, digits = 12), format(tail_exponent, digits = 12), format(gdar_display_epsilon, digits = 12),
            as.character(gdar_display_theory), ifelse(is.finite(gdar_display_max_x), format(gdar_display_max_x, digits = 12), ""),
            as.character(monotone_wf), as.character(edge_mode), as.character(exclude_zero_deme_endpoint), as.character(all_wf_points),
            nrow(curve_raw_export), nrow(curve_display), nrow(point_export), nrow(point_plot),
            ifelse(is.finite(z_gdar), format(z_gdar, digits = 12), ""),
            as.character(include_equilibrium), format(x_max, digits = 12)),
  stringsAsFactors = FALSE
)
utils::write.table(diag_rows, diagnostic_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# ------------------------- plot -------------------------
series_levels <- intersect(c("GDAR", "WF3", "WF5", "EQUILIBRIUM"), unique(curve_display$series))
series_labels <- vapply(series_levels, function(s) series_specs[[s]]$label, character(1))
series_colors <- vapply(series_levels, function(s) series_specs[[s]]$color, character(1))
series_linetypes <- vapply(series_levels, function(s) series_specs[[s]]$linetype, character(1))
names(series_labels) <- names(series_colors) <- names(series_linetypes) <- series_levels

finite_y <- c(curve_display$value_pct, if (nrow(point_plot)) point_plot$value_pct else numeric())
finite_y <- finite_y[is.finite(finite_y)]
if (tolower(y_max_arg) == "auto") {
  # Leave a small, stable headroom for the two-line scenario labels.
  y_max <- if (length(finite_y)) max(115, ceiling(max(finite_y) / 5) * 5) else 115
} else {
  y_max <- suppressWarnings(as.numeric(y_max_arg))
  if (!is.finite(y_max) || y_max <= 0) stop("--y-max must be positive or auto", call. = FALSE)
}
x_breaks <- seq(0, x_max, length.out = 5L)
y_breaks <- if (y_max >= 100) seq(0, 100, by = 25) else seq(0, y_max, length.out = 5L)

line_data <- curve_display[curve_display$series %in% series_levels, , drop = FALSE]
point_data <- point_plot[point_plot$series %in% series_levels, , drop = FALSE]
if (nrow(point_data)) {
  point_period <- scenario_period_class(point_data$scenario)
  # GDAR points are circles; WF points use triangle/square to distinguish the
  # early (2061-2080) and late (2081-2100) scenarios as in the reference figure.
  point_data$point_shape <- ifelse(point_data$series == "GDAR", 16L,
                                   ifelse(point_period == "early", 17L,
                                          ifelse(point_period == "late", 15L,
                                                 ifelse(point_data$series == "WF3", 15L,
                                                        ifelse(point_data$series == "WF5", 17L, 18L)))))
}

p <- ggplot()
if (isTRUE(plot_raw) && line_method != "raw") {
  raw_line_plot <- curve_raw_export[
    is.finite(curve_raw_export$habitat_loss_pct) &
      is.finite(curve_raw_export$value_pct) &
      curve_raw_export$series %in% series_levels &
      (!exclude_zero_deme_endpoint |
         !is.finite(curve_raw_export$occupied_demes) |
         curve_raw_export$occupied_demes > 0), , drop = FALSE
  ]
  p <- p + geom_line(
    data = raw_line_plot,
    aes(x = habitat_loss_pct, y = value_pct, colour = series),
    linewidth = 0.5, alpha = 0.28, show.legend = FALSE, na.rm = TRUE
  )
}
p <- p + geom_line(data = line_data,
            aes(x = habitat_loss_pct, y = value_pct, colour = series, linetype = series),
            linewidth = 1.15, na.rm = TRUE)
if (nrow(point_data)) {
  p <- p + geom_point(data = point_data,
                      aes(x = habitat_loss_pct, y = value_pct, shape = point_shape, colour = series, fill = series),
                      size = 4.2, stroke = 0.7, show.legend = FALSE, na.rm = TRUE)
}

if (!no_labels && nrow(point_data) && label_series %in% series_levels) {
    labels_df <- point_data[point_data$series == label_series &
                              !is.na(point_data$scenario_label) &
                              nzchar(point_data$scenario_label), , drop = FALSE]
  if (nrow(labels_df)) {
    # Keep labels tied to their point and in the same left-to-right order as
    # the reference panel. Alternating two fixed rows avoids ggrepel moving a
    # label over a neighbouring scenario when losses are close together.
    labels_df <- labels_df[order(labels_df$habitat_loss_pct, labels_df$scenario), , drop = FALSE]
    labels_df$label_y <- ifelse(seq_len(nrow(labels_df)) %% 2L == 1L, 103, 109)
    labels_df$label_y <- pmin(labels_df$label_y, y_max - 2)
    p <- p + geom_text(data = labels_df,
                       aes(x = habitat_loss_pct, y = label_y, label = scenario_label),
                       vjust = 0, hjust = 0.5, size = 3.25, lineheight = 0.9,
                       show.legend = FALSE)
  }
}

p <- p +
  scale_colour_manual(values = series_colors, breaks = series_levels, labels = series_labels, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = series_linetypes, breaks = series_levels, labels = series_labels, name = NULL, drop = FALSE) +
  scale_shape_identity(guide = "none") +
  scale_fill_manual(values = series_colors, guide = "none") +
  scale_x_continuous(breaks = x_breaks, expand = c(0, 0)) +
  scale_y_continuous(breaks = y_breaks, expand = c(0, 0)) +
  labs(x = "Habitat loss (%)", y = "Predicted genetic diversity remaining (%)",
       title = if (nzchar(plot_title)) plot_title else NULL,
       subtitle = NULL) +
  coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max), expand = FALSE, clip = "on") +
  theme_bw(base_size = 15) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.7, colour = "black"),
    legend.position = c(0.20, 0.20),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 11),
    plot.margin = margin(10, 16, 10, 10)
  )

png_path <- file.path(outdir, paste0(prefix, ".png"))
pdf_path <- file.path(outdir, paste0(prefix, ".pdf"))
# Explicitly select the base R device.  Some R 4.5 installations register an
# ``agg`` default that fails on network-mounted paths even though PNG support is
# available through grDevices::png.
ggsave(png_path, p, width = 10.5, height = 7.8, dpi = 400,
       device = grDevices::png, bg = "white")
ggsave(pdf_path, p, width = 10.5, height = 7.8, bg = "white")

cat("WFmoment plot written:\n")
cat(" - ", png_path, "\n", sep = "")
cat(" - ", pdf_path, "\n", sep = "")
cat(" - ", curve_export_path, " (raw curve values)\n", sep = "")
cat(" - ", display_export_path, " (display line values)\n", sep = "")
cat(" - ", point_export_path, " (raw scenario values)\n", sep = "")
cat("Series: ", paste(series_levels, collapse = ", "), "\n", sep = "")
cat("GDAR source: ", ifelse(is.finite(z_gdar) && (is.null(gdar_curve_obj) || grepl("derived", curve_objs$GDAR$source_column)),
                              paste0("derived (z_gdar=", format(z_gdar, digits = 8), ")"),
                              if ("GDAR" %in% names(curve_objs)) curve_objs$GDAR$source_column else "not available"), "\n", sep = "")
