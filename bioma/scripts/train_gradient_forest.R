#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("usage: train_gradient_forest.R TRAINING.tsv OUTDIR NTREE NBIN CORR MAXLEVEL SEED PREDICTORS")
}

suppressPackageStartupMessages(library(gradientForest))

training_path <- args[[1]]
output_dir <- args[[2]]
ntree <- as.integer(args[[3]])
nbin <- as.integer(args[[4]])
corr_threshold <- as.numeric(args[[5]])
max_level <- as.numeric(args[[6]])
seed <- as.integer(args[[7]])
predictors <- strsplit(args[[8]], ",", fixed = TRUE)[[1]]

training <- read.delim(training_path, check.names = FALSE)
required <- c("population_id", predictors)
if (!all(required %in% names(training))) {
  stop("training table is missing population_id or predictor columns")
}
responses <- setdiff(names(training), required)
if (length(responses) == 0) {
  stop("training table contains no response loci")
}
model_data <- training[, c(predictors, responses), drop = FALSE]

set.seed(seed)
all_gfmod <- gradientForest(
  data = model_data,
  predictor.vars = predictors,
  response.vars = responses,
  ntree = ntree,
  compact = TRUE,
  nbin = nbin,
  maxLevel = max_level,
  trace = TRUE,
  corr.threshold = corr_threshold,
  check.names = FALSE
)

save(all_gfmod, file = file.path(output_dir, "all_gfmod.data"))

importance <- data.frame(
  predictor = predictors,
  importance = as.numeric(all_gfmod$overall.imp[predictors]),
  stringsAsFactors = FALSE
)
write.table(
  importance,
  file = file.path(output_dir, "predictor_importance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

positive_names <- names(all_gfmod$result)
positive_r2 <- rep(NA_real_, length(responses))
names(positive_r2) <- responses
positive_r2[positive_names] <- as.numeric(all_gfmod$result)
performance <- data.frame(
  response = responses,
  positive_oob_r2 = as.numeric(positive_r2),
  retained_positive_r2 = responses %in% positive_names,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.table(
  performance,
  file = file.path(output_dir, "response_performance.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

transformed <- predict(all_gfmod, training[, predictors, drop = FALSE])
transformed <- data.frame(
  population_id = training$population_id,
  transformed,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write.table(
  transformed,
  file = file.path(output_dir, "training_transformed_environment.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

capture.output(sessionInfo(), file = file.path(output_dir, "r_session_info.txt"))
summary <- data.frame(
  key = c("r_version", "gradientforest_version", "positive_responses"),
  value = c(
    as.character(getRversion()),
    as.character(packageVersion("gradientForest")),
    as.character(length(positive_names))
  ),
  stringsAsFactors = FALSE
)
write.table(
  summary,
  file = file.path(output_dir, ".r_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
