#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 9) stop("usage: rona_plot.R output_dir environment mask models ssps periods grid_step grid_k interpolate")
root <- args[[1]]; env_path <- args[[2]]; mask_path <- args[[3]]
models <- if (nzchar(args[[4]])) strsplit(args[[4]], ",", fixed=TRUE)[[1]] else character()
ssps <- if (nzchar(args[[5]])) strsplit(args[[5]], ",", fixed=TRUE)[[1]] else character()
periods <- if (nzchar(args[[6]])) strsplit(args[[6]], ",", fixed=TRUE)[[1]] else character()
grid_step <- as.numeric(args[[7]]); grid_k <- as.integer(args[[8]]); interpolate <- args[[9]] != "0"
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(mgcv); library(raster) })
population_palette <- function(values) {
  groups <- sort(unique(as.character(values)))
  groups <- groups[nzchar(groups) & !is.na(groups)]
  if (!length(groups)) return(setNames(character(), character()))
  setNames(grDevices::hcl.colors(length(groups), palette = "Dark 3"), groups)
}

env <- fread(env_path, data.table=FALSE)
if (!all(c("ID","pop","lon","lat") %in% names(env))) stop("environment must contain ID,pop,lon,lat")
if (!length(models)) models <- list.dirs(file.path(root,"weighted"), full.names=FALSE, recursive=FALSE)
if (!length(ssps)) ssps <- unique(sub("^ssp", "", sub(".*_(ssp[0-9]+)_RONA.*", "\\1", basename(list.files(file.path(root,"weighted"), recursive=TRUE, full.names=FALSE)))))
if (!length(periods)) periods <- unique(sub("_ssp.*", "", basename(list.files(file.path(root,"weighted"), recursive=TRUE, full.names=FALSE))))
models <- models[nzchar(models)]; ssps <- ssps[nzchar(ssps)]; periods <- periods[nzchar(periods)]

read_ensemble <- function(period, ssp, bio) {
  values <- list(); errors <- list()
  for (model in models) {
    wf <- file.path(root,"weighted",model,paste0(period,"_",model,"_ssp",ssp,"_RONA_weighted.csv"))
    sf <- file.path(root,"SE",model,paste0(period,"_",model,"_ssp",ssp,"_RONA_SE.csv"))
    if (!file.exists(wf) || !file.exists(sf)) next
    w <- fread(wf, data.table=FALSE); e <- fread(sf, data.table=FALSE)
    key <- paste0("ssp",ssp,"_BIO",bio)
    if (!(key %in% names(w)) || !(key %in% names(e))) next
    values[[model]] <- w[,c("ID",key)]; errors[[model]] <- e[,c("ID",key)]
    names(values[[model]])[2] <- model; names(errors[[model]])[2] <- model
  }
  # Ensemble plots require every selected model. This prevents a missing model
  # from changing the spatial mean or receiving less weight silently.
  if (length(values) != length(models) || length(values) < 2) return(NULL)
  joined <- Reduce(function(x,y) merge(x,y,by="ID",all=TRUE), values)
  joined_e <- Reduce(function(x,y) merge(x,y,by="ID",all=TRUE), errors)
  val <- as.matrix(joined[,-1,drop=FALSE]); err <- as.matrix(joined_e[,-1,drop=FALSE]);
  out <- env[match(joined$ID, env$ID),c("ID","pop","lon","lat")]
  out$rona <- rowMeans(val, na.rm=TRUE); out$se <- rowMeans(err, na.rm=TRUE)
  out$n_models <- rowSums(is.finite(val)); out <- out[out$n_models == length(models) & is.finite(out$rona),]
  out$period <- period; out$ssp <- paste0("ssp",ssp); out
}

mask_sp <- NULL; mask_df <- NULL
if (nzchar(mask_path) && file.exists(mask_path)) {
  mask_sp <- tryCatch(raster::shapefile(mask_path), error=function(e) NULL)
  if (!is.null(mask_sp)) mask_df <- ggplot2::fortify(mask_sp)
}
map_one <- function(dat, bio, out_base) {
  xlim <- range(dat$lon, na.rm=TRUE); ylim <- range(dat$lat, na.rm=TRUE)
  if (!is.null(mask_sp)) { bb <- bbox(mask_sp); xlim <- c(bb[1,1],bb[1,2]); ylim <- c(bb[2,1],bb[2,2]) }
  tile <- NULL
  if (interpolate && nrow(dat) >= 5) {
    fit <- tryCatch(gam(rona ~ s(lon,lat,k=min(grid_k,nrow(dat)-1)), data=dat, method="REML"), error=function(e) NULL)
    if (!is.null(fit)) {
      grd <- expand.grid(lon=seq(xlim[1],xlim[2],by=grid_step), lat=seq(ylim[1],ylim[2],by=grid_step)); grd$rona <- as.numeric(predict(fit,newdata=grd))
      rr <- rasterFromXYZ(grd[,c("lon","lat","rona")]);
      if (!is.null(mask_sp)) rr <- tryCatch(mask(crop(rr,mask_sp),mask_sp), error=function(e) rr)
      tile <- as.data.frame(rr,xy=TRUE); names(tile)[3] <- "value"; tile <- tile[is.finite(tile$value),]
      if (nrow(tile)) { lo <- min(tile$value); hi <- max(tile$value); if (lo==hi) tile$class <- factor("value") else tile$class <- cut(tile$value,breaks=seq(lo,hi,length.out=6),include.lowest=TRUE) }
    }
  }
  p <- ggplot()
  if (!is.null(mask_df)) p <- p + geom_polygon(data=mask_df,aes(x=long,y=lat,group=group),fill="#f2f2f2",color="black",linewidth=.3)
  if (!is.null(tile) && nrow(tile)) p <- p + geom_tile(data=tile,aes(x=x,y=y,fill=class),width=grid_step,height=grid_step) + scale_fill_manual(values=c("#E3EEF8","#BFD7EA","#73ADD2","#3E7BB6","#1F4E8C"),name="RONA",drop=FALSE)
  p <- p + geom_point(data=dat,aes(x=lon,y=lat,color=pop),size=2.4,shape=21,fill="white",stroke=.6) + scale_color_manual(values=population_palette(dat$pop),drop=FALSE) + coord_cartesian(xlim=xlim,ylim=ylim,expand=FALSE) + labs(title=paste0("RONA BIO",bio," ",dat$period[1]," ",toupper(dat$ssp[1])),x="Longitude",y="Latitude",color="Population") + theme_bw() + theme(panel.grid=element_blank(),text=element_text(family="serif"),plot.title=element_text(hjust=.5),legend.position="right")
  ggsave(paste0(out_base,".png"),p,width=6.67,height=5.33,dpi=300,bg="white"); ggsave(paste0(out_base,".pdf"),p,width=6.67,height=5.33,bg="white")
}

dir.create(file.path(root,"maps"),recursive=TRUE,showWarnings=FALSE); dir.create(file.path(root,"boxplots"),recursive=TRUE,showWarnings=FALSE); dir.create(file.path(root,"ensemble_mean"),recursive=TRUE,showWarnings=FALSE)

# Save the exact coordinate-level model means used by both plot types.
for (period in periods) for (ssp in ssps) {
  merged <- NULL
  for (bio in 1:19) {
    d <- read_ensemble(period,ssp,bio); if (is.null(d)) next
    x <- d[,c("ID","pop","lon","lat","rona","se","n_models")]
    names(x)[5:6] <- c(paste0("BIO",bio,"_RONA"), paste0("BIO",bio,"_SE"))
    merged <- if (is.null(merged)) x else merge(merged,x,by=c("ID","pop","lon","lat","n_models"),all=TRUE)
  }
  if (!is.null(merged)) fwrite(merged,file.path(root,"ensemble_mean",paste0(period,"_ssp",ssp,"_RONA_ensemble_mean.tsv")),sep="\t")
}

for (bio in 1:19) {
  all <- list()
  for (period in periods) for (ssp in ssps) { d <- read_ensemble(period,ssp,bio); if (!is.null(d)) { all[[length(all)+1]] <- d; map_one(d,bio,file.path(root,"maps",paste0("bio",bio,".",period,".ssp",ssp,".rona"))) } }
  if (!length(all)) next
  box <- rbindlist(all,fill=TRUE); latest <- periods[length(periods)]; box <- box[box$period==latest & box$ssp %in% paste0("ssp",c("245","585")),]
  if (!nrow(box)) box <- rbindlist(all,fill=TRUE)
  p <- ggplot(box,aes(x=factor(ID,levels=unique(ID)),y=rona,color=ssp)) + geom_point(position=position_dodge(width=.5),size=2.4) + geom_errorbar(aes(ymin=pmax(0,rona-se),ymax=rona+se),position=position_dodge(width=.5),width=.35) + scale_color_manual(values=c(ssp245="#F8766D",ssp585="#31C0C4"),drop=FALSE) + labs(x="Population",y=paste0("BIO",bio,"_RONA (mean +/- SE)"),color="SSP") + theme_bw() + theme(text=element_text(family="serif"),axis.text.x=element_text(angle=90,hjust=1,size=7),panel.grid=element_blank(),legend.position="top")
  ggsave(file.path(root,"boxplots",paste0("bio",bio,".rona.boxplot.png")),p,width=6,height=6,dpi=300,bg="white"); ggsave(file.path(root,"boxplots",paste0("bio",bio,".rona.boxplot.pdf")),p,width=6,height=6,bg="white")
}
