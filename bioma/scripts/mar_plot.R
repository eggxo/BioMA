#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop("usage: mar_plot.R output_dir name [scenario_file] [geom_id]")
root <- args[[1]]; name <- args[[2]]; scenario_file <- if (length(args)>=3) args[[3]] else ""; geom_id <- if (length(args)>=4) as.integer(args[[4]]) else 7
suppressPackageStartupMessages({ library(mar); library(data.table); library(ggplot2); library(dplyr); library(tidyr); library(tibble) })
rda_marker <- file.path(root,"mar_extdflist_path.txt")
rda <- if (file.exists(rda_marker)) readLines(rda_marker,n=1) else list.files(root,pattern="^extdflist.*\\.rda$",recursive=TRUE,full.names=TRUE)[1]
if (!length(rda) || !file.exists(rda)) stop("Cannot locate extdflist RDA under ",root)
load(rda); if (!exists("extdflist")) stop("RDA does not contain extdflist")
scheme_name <- if (!is.null(extdflist$random)) "random" else names(extdflist)[1]; extdf <- extdflist[[scheme_name]]
if (is.null(extdf) || !nrow(extdf)) stop("Selected MAR extdflist is empty")
set.seed(123)

extract_coef <- function(fit) {
  txt <- capture.output(print(fit)); start <- grep("^Coefficients:",txt)
  if (!length(start)) stop("MARcalc output has no Coefficients block")
  lines <- txt[(start[1]+1):length(txt)]; nums <- unlist(regmatches(lines,gregexpr("[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?",lines,perl=TRUE))); nums <- suppressWarnings(as.numeric(nums)); nums <- nums[is.finite(nums)]
  if (length(nums)<2) stop("Unable to extract MAR c and z"); tibble(c_coef=nums[1],z=nums[2])
}
fits <- list(M=MARcalc(extdf,Mtype="M",Atype="A"), E=MARcalc(extdf,Mtype="E",Atype="A"), thetaW=MARcalc(extdf,Mtype="thetaw",Atype="A"), thetaPi=MARcalc(extdf,Mtype="thetapi",Atype="A"))
z_tbl <- rbindlist(lapply(names(fits), function(m) cbind(metric=m,as.data.frame(extract_coef(fits[[m]])))),fill=TRUE); z_tbl <- as.data.frame(z_tbl); z_tbl$c_coef <- as.numeric(z_tbl$c_coef); z_tbl$z <- as.numeric(z_tbl$z); z_tbl$model <- "MARcalc_power"; z_tbl$formula_raw <- "S = c * A^z"; z_tbl$formula_relative <- "D_remaining = (1 - habitat_loss)^z"
fwrite(z_tbl,file.path(root,"MAR_official_fit_params.tsv"),sep="\t")
curve_one <- function(metric,z) data.frame(metric=metric,habitat_loss_prop=seq(0,.999,length.out=1500),habitat_loss_pct=seq(0,.999,length.out=1500)*100,A_remaining_prop=1-seq(0,.999,length.out=1500),A_remaining_pct=(1-seq(0,.999,length.out=1500))*100,diversity_remaining_pct=((1-seq(0,.999,length.out=1500))^z)*100)
curve <- rbindlist(lapply(seq_len(nrow(z_tbl)),function(i) curve_one(z_tbl$metric[i],z_tbl$z[i]))); fwrite(curve,file.path(root,"MAR_display_curves.tsv"),sep="\t")

scenario <- NULL
if (nzchar(scenario_file) && file.exists(scenario_file)) {
  s <- fread(scenario_file)
  if (all(c("scenario","A_remaining","geom_job_id") %in% names(s))) {
    info <- data.frame(scenario=c("2061_2080_ssp245","2061_2080_ssp585","2081_2100_ssp245","2081_2100_ssp585"),period=c("2061-2080","2061-2080","2081-2100","2081-2100"),ssp=c("SSP245","SSP585","SSP245","SSP585"),scenario_label=c("2061-2080 SSP245","2061-2080 SSP585","2081-2100 SSP245","2081-2100 SSP585"))
    scenario <- merge(s[geom_job_id==geom_id,.(A_remaining=first(A_remaining)),by=scenario],info,by="scenario",all.x=TRUE)
    scenario <- scenario[!is.na(scenario_label)]
    mx <- max(scenario$A_remaining,na.rm=TRUE); scenario$A_remaining_prop <- if(mx<=1.5) scenario$A_remaining else if(mx<=100) scenario$A_remaining/100 else scenario$A_remaining/mx
    scenario$habitat_loss_prop <- 1-scenario$A_remaining_prop; scenario$habitat_loss_pct <- scenario$habitat_loss_prop*100
    points <- merge(as.data.table(scenario),z_tbl[,.(metric,z)],allow.cartesian=TRUE); points$diversity_remaining_pct <- ((1-points$habitat_loss_prop)^points$z)*100
    fwrite(points,file.path(root,"MAR_scenario_points.tsv"),sep="\t")
  }
}

curve$metric <- factor(curve$metric,levels=c("M","E","thetaW","thetaPi")); p <- ggplot(curve,aes(habitat_loss_pct,diversity_remaining_pct,color=metric))+geom_line(linewidth=1.2)+scale_x_continuous(limits=c(0,100),breaks=seq(0,100,25),expand=c(0,0))+scale_y_continuous(limits=c(0,105),breaks=seq(0,100,25),expand=c(0,0))+scale_color_manual(values=c(M="grey35",E="#D95F02",thetaW="#1B9E77",thetaPi="#7570B3"))+labs(x="Habitat loss (%)",y="Predicted genetic diversity remaining (%)",color=NULL)+theme_bw(base_size=15)+theme(panel.grid=element_blank(),axis.title=element_text(face="bold"),legend.position=c(.18,.2))
if (!is.null(scenario)) { points$metric <- factor(points$metric,levels=levels(curve$metric)); p <- p+geom_point(data=points[period=="2061-2080"],aes(fill=metric),shape=24,size=3.5,show.legend=FALSE)+geom_point(data=points[period=="2081-2100"],aes(fill=metric),shape=22,size=3.8,show.legend=FALSE) }
ggsave(file.path(root,"MAR_habitat_loss_curves.png"),p,width=10.5,height=7.8,dpi=400,bg="white"); ggsave(file.path(root,"MAR_habitat_loss_curves.pdf"),p,width=10.5,height=7.8,bg="white")
