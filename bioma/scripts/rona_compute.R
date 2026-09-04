#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) stop("usage: rona_compute.R alt_frequency unld_dir environment future_root output_dir models ssps periods")
alt_path <- args[[1]]; unld_dir <- args[[2]]; env_path <- args[[3]]; future_root <- args[[4]]; out_root <- args[[5]]
models <- if (nzchar(args[[6]])) strsplit(args[[6]], ",", fixed=TRUE)[[1]] else character()
ssps <- if (nzchar(args[[7]])) strsplit(args[[7]], ",", fixed=TRUE)[[1]] else character()
periods <- if (nzchar(args[[8]])) strsplit(args[[8]], ",", fixed=TRUE)[[1]] else character()

suppressPackageStartupMessages({ library(data.table); library(plotrix); library(terra) })
dir.create(file.path(out_root, "weighted"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out_root, "SE"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out_root, "weighted_SE"), recursive=TRUE, showWarnings=FALSE)

freq <- fread(alt_path, data.table=FALSE, check.names=TRUE)
if (!"Pop" %in% names(freq)) names(freq)[1] <- "Pop"
env <- fread(env_path, data.table=FALSE, check.names=FALSE)
if (!all(c("ID", "pop", "lon", "lat") %in% names(env))) stop("environment must contain ID,pop,lon,lat")
env$ID <- as.character(env$ID); freq$Pop <- as.character(freq$Pop)
keep <- env$ID[env$ID %in% freq$Pop]
if (length(keep) < 3) stop("fewer than three overlapping populations")
env <- env[match(keep, env$ID), , drop=FALSE]; freq <- freq[match(keep, freq$Pop), , drop=FALSE]
popname <- env[, c("ID", "pop", "lon", "lat"), drop=FALSE]

rona_pred <- function(gen, present, future) {
  gen <- as.matrix(gen); storage.mode(gen) <- "double"; n <- nrow(gen); m <- ncol(gen)
  rona <- matrix(NA_real_, n, m); rsq <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    ok <- is.finite(gen[,j]) & is.finite(present)
    if (sum(ok) < 3 || length(unique(present[ok])) < 2) next
    fit <- lm(gen[ok,j] ~ present[ok]); cf <- coef(fit); rsq[j] <- summary(fit)$r.squared
    rona[,j] <- abs(cf[[2]] * future + cf[[1]] - gen[,j])
  }
  avg <- se <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    x <- rona[i,]; ok <- is.finite(x)
    if (any(ok)) { avg[i] <- weighted.mean(x[ok], rsq[ok], na.rm=TRUE); se[i] <- if (sum(ok)>1) sd(x[ok])/sqrt(sum(ok)) else NA_real_ }
  }
  list(avg=avg, se=se)
}

scenario_dirs <- list.dirs(future_root, recursive=TRUE, full.names=TRUE)
scenario_dirs <- scenario_dirs[grepl("^[0-9]{4}-[0-9]{4}-ssp[0-9]+-.+", basename(scenario_dirs), ignore.case=TRUE)]
for (directory in scenario_dirs) {
  scenario <- trimws(basename(directory)); parts <- regmatches(scenario, regexec("^([0-9]{4}-[0-9]{4})-(ssp[0-9]+)-(.+)$", scenario, ignore.case=TRUE))[[1]]
  if (length(parts) != 4) next
  period <- parts[[2]]; ssp <- tolower(parts[[3]]); model <- parts[[4]]
  if (length(models) && !(model %in% models)) next
  if (length(ssps) && !(sub("^ssp", "", ssp) %in% ssps)) next
  if (length(periods) && !(period %in% periods)) next
  tif_files <- list.files(directory, pattern="bio[0-9]+\\.cut\\.tif$", full.names=TRUE, ignore.case=TRUE)
  if (length(tif_files) < 19) { message("skip incomplete scenario: ", scenario); next }
  future_env <- matrix(NA_real_, nrow=nrow(env), ncol=19, dimnames=list(env$ID, paste0("BIO",1:19)))
  points <- vect(env[, c("lon", "lat")], geom=c("lon", "lat"), crs="EPSG:4326")
  for (j in 1:19) {
    hit <- tif_files[grepl(paste0("bio", j, "\\.cut\\.tif$"), tif_files, ignore.case=TRUE)]
    if (length(hit) != 1) stop("BIO", j, " not found or duplicated in ", directory)
    future_env[,j] <- terra::extract(rast(hit), points)[,2]
  }
  result_w <- popname; result_se <- popname; result_ws <- popname
  for (j in 1:19) {
    prune <- file.path(unld_dir, paste0("LD_BIO", j, ".prune.in")); if (!file.exists(prune)) stop("missing prune file: ", prune)
    ids <- gsub(":", ".", trimws(readLines(prune, warn=FALSE))); ids <- ids[nzchar(ids)]
    loci <- intersect(ids, names(freq)[-1]); if (!length(loci)) stop("no frequency columns overlap ", prune)
    fit <- rona_pred(freq[, loci, drop=FALSE], env[[paste0("bio",j)]], future_env[,j]); key <- paste0(ssp, "_BIO", j)
    result_w[[key]] <- fit$avg; result_se[[key]] <- fit$se
    result_ws[[paste0(key, "(SE)")]] <- ifelse(is.finite(fit$avg), paste0(round(fit$avg,4), "±(", round(fit$se,4), ")"), NA)
  }
  dir.create(file.path(out_root, "weighted", model), recursive=TRUE, showWarnings=FALSE); dir.create(file.path(out_root, "SE", model), recursive=TRUE, showWarnings=FALSE); dir.create(file.path(out_root, "weighted_SE", model), recursive=TRUE, showWarnings=FALSE)
  fwrite(result_w, file.path(out_root, "weighted", model, paste0(period,"_",model,"_",ssp,"_RONA_weighted.csv")), quote=FALSE)
  fwrite(result_se, file.path(out_root, "SE", model, paste0(period,"_",model,"_",ssp,"_RONA_SE.csv")), quote=FALSE)
  fwrite(result_ws, file.path(out_root, "weighted_SE", model, paste0(period,"_",model,"_",ssp,"_RONA_weightedSE.csv")), quote=FALSE)
  message("completed ", scenario)
}
