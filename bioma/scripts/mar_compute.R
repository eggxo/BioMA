#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 11) stop("usage: mar_compute.R name workdir vcf lonlat scheme nrep xfrac quorum randseed maxsnps marsteps")
name <- args[[1]]; workdir <- args[[2]]; vcf <- args[[3]]; lonlat <- args[[4]]
scheme <- tolower(args[[5]]); nrep <- as.integer(args[[6]]); xfrac <- as.numeric(args[[7]])
quorum <- toupper(args[[8]]) == "TRUE"; randseed <- as.integer(args[[9]]); maxsnps <- as.integer(args[[10]])
marsteps <- strsplit(args[[11]], ",", fixed=TRUE)[[1]]
suppressPackageStartupMessages(library(mar))
dir.create(workdir, recursive=TRUE, showWarnings=FALSE)
cat("MAR version:", as.character(packageVersion("mar")), "\n")
cat("name:", name, "\nworkdir:", workdir, "\nvcf:", vcf, "\nlonlat:", lonlat, "\nscheme:", scheme, "\nnrep:", nrep, "\nxfrac:", xfrac, "\nmaxsnps:", maxsnps, "\n")
set.seed(randseed)
MARPIPELINE(
  name=name,
  workdir=workdir,
  genofile=vcf,
  lonlatfile=lonlat,
  filetype="vcf",
  randseed=randseed,
  option_geno=list(ploidy=2, maxsnps=maxsnps),
  option_map=list(mapres=NULL, mapcrs="OGC:CRS84"),
  option_marext=list(scheme=c(scheme), nrep=nrep, xfrac=xfrac, quorum=quorum, animate=FALSE),
  marsteps=marsteps,
  saveobj=TRUE
)
rda <- list.files(workdir, pattern="^extdflist.*\\.rda$", recursive=TRUE, full.names=TRUE)
if (!length(rda)) stop("MAR pipeline finished without an extdflist RDA")
writeLines(rda[[which.max(file.info(rda)$mtime)]], file.path(workdir, "mar_extdflist_path.txt"))
cat("extdflist:", rda[[which.max(file.info(rda)$mtime)]], "\n")
