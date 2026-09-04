#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) stop("usage: load_plot.R output_dir mask targets")
root <- args[[1]]; mask_path <- args[[2]]; targets <- strsplit(args[[3]], ",", fixed=TRUE)[[1]]
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(mgcv); library(raster); library(sf); library(cowplot); library(patchwork)})
if (!file.exists(mask_path)) stop("species mask does not exist: ", mask_path)
mask_sf <- st_read(mask_path, quiet=TRUE)
# Some legacy masks have no CRS metadata; project coordinates are lon/lat.
mask_crs <- st_crs(mask_sf)
if (is.na(mask_crs)) {
  mask_sf <- st_set_crs(mask_sf, 4326)
} else if (mask_crs$epsg != 4326) {
  mask_sf <- st_transform(mask_sf, 4326)
}
mask_sp <- as(mask_sf, "Spatial")
mask_df <- fortify(mask_sp)
group_palette <- function(values) {
  groups <- sort(unique(as.character(values)))
  groups <- groups[nzchar(groups) & !is.na(groups)]
  if (!length(groups)) return(setNames(character(), character()))
  cols <- grDevices::hcl.colors(length(groups), palette = "Dark 3")
  setNames(cols, groups)
}
pred_files <- list.files(file.path(root,"future_predictions"), pattern="^future_env_ssp[0-9]+_[0-9]{4}-[0-9]{4}_mean\\.(mean_load[MD]_relax)\\.csv$", full.names=TRUE)
if (!length(pred_files)) stop("No future prediction files found")
parse_meta <- function(f) { n <- basename(f); m <- regexec("^future_env_(ssp[0-9]+)_([0-9]{4}-[0-9]{4})_mean\\.(mean_load[MD]_relax)\\.csv$",n); z <- regmatches(n,m)[[1]]; if(length(z)!=4) return(NULL); list(ssp=z[2],period=z[3],target=z[4]) }
make_map <- function(f, target) {
  meta <- parse_meta(f); d <- fread(f, data.table=FALSE); if (!all(c("lon","lat","cluster","predicted_value") %in% names(d))) stop("prediction columns missing in ",f)
  d$lon <- as.numeric(d$lon); d$lat <- as.numeric(d$lat); d$value <- as.numeric(d$predicted_value); d <- d[is.finite(d$value) & is.finite(d$lon) & is.finite(d$lat),]
  if (nrow(d)<10) return(NULL)
  bb <- st_bbox(mask_sf); xlim <- c(as.numeric(bb["xmin"])-.5, as.numeric(bb["xmax"])+.5); ylim <- c(as.numeric(bb["ymin"])-.5, as.numeric(bb["ymax"])+.5)
  fit <- gam(value ~ s(lon,lat,k=min(10,nrow(d)-1)), data=d, method="REML")
  grd <- expand.grid(lon=seq(xlim[1],xlim[2],by=.05), lat=seq(ylim[1],ylim[2],by=.05)); grd$value <- as.numeric(predict(fit,newdata=grd)); rr <- rasterFromXYZ(grd[,c("lon","lat","value")]); rr <- tryCatch(mask(crop(rr,mask_sp),mask_sp),error=function(e) rr); tile <- as.data.frame(rr,xy=TRUE); names(tile)[3] <- "value"; tile <- tile[is.finite(tile$value),]
  if (nrow(tile)) { r <- range(tile$value); if (r[1]==r[2]) r[2] <- r[1]+1e-8; tile$class <- cut(tile$value,breaks=seq(r[1],r[2],length.out=6),include.lowest=TRUE,dig.lab=6) }
  p <- ggplot() + geom_polygon(data=mask_df,aes(long,lat,group=group),fill="#f2f2f2",color="black",linewidth=.3) + geom_tile(data=tile,aes(x,y,fill=class),width=.05,height=.05) + scale_fill_manual(values=c("#E6DCEB","#C9ADD8","#A26BBE","#6A3D9A","#3F007D"),name=ifelse(target=="mean_loadM_relax","Predicted loadM","Predicted loadD"),drop=FALSE) + geom_point(data=d,aes(lon,lat,color=cluster),size=2.3,shape=21,fill="white",stroke=.55) + scale_color_manual(values=group_palette(d$cluster),drop=FALSE) + coord_sf(xlim=xlim,ylim=ylim,expand=FALSE) + labs(title=paste0(ifelse(target=="mean_loadM_relax","Predicted loadM","Predicted loadM(DEL)")," (",gsub("ssp","SSP",toupper(meta$ssp)),", ",meta$period,")"),x="Longitude",y="Latitude",color="Population") + theme_bw() + theme(panel.grid=element_blank(),plot.title=element_text(hjust=.5,face="bold",size=12),legend.position="right",text=element_text(family="serif"))
  p
}
map_dir <- file.path(root,"maps"); dir.create(map_dir,recursive=TRUE,showWarnings=FALSE)
plots <- list(); labels <- character()
for (f in pred_files) { meta <- parse_meta(f); if (is.null(meta) || !(meta$target %in% targets)) next; p <- make_map(f,meta$target); if (is.null(p)) next; stub <- file.path(map_dir,paste0(meta$target,"_",meta$ssp,"_",meta$period)); ggsave(paste0(stub,".png"),p,width=6.4,height=5.8,dpi=300,bg="white"); ggsave(paste0(stub,".pdf"),p,width=6.4,height=5.8,bg="white"); plots[[length(plots)+1]] <- p; labels <- c(labels,paste(meta$target,meta$ssp,meta$period)) }
if (length(plots)) { combo <- plot_grid(plotlist=plots,ncol=2,labels=LETTERS[seq_along(plots)]); ggsave(file.path(root,"Figure_S25_load_maps.png"),combo,width=13,height=5.8*ceiling(length(plots)/2)/1.2,dpi=300,bg="white"); ggsave(file.path(root,"Figure_S25_load_maps.pdf"),combo,width=13,height=5.8*ceiling(length(plots)/2)/1.2,bg="white") }

# S24 tuning panels are promoted from the RF script's reference-style plots by
# the Python orchestrator, one file for each target. This plotting script is
# responsible only for the S25 spatial maps.
