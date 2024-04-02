reticulate::use_condaenv("torch-gpu")
dyn.load('/opt/R/arm64/gfortran/lib/libgcc_s.2.dylib')
dyn.load('/opt/R/arm64/gfortran/lib/libgfortran.5.dylib')
library(SPATA2)
library(Seurat)
library(EBImage)
library(tidyverse)
library(reticulate)
library(tidyverse)
library(EBImage)
library(tidyverse)
list2array = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}
getPatches <- function(img, patch=c(160,160)){
  
  nr_patches <- c(c(img %>% dim())[1:2]/patch[1]) %>% round(digits = 0)
  
  img <- EBImage::resize(img, w=patch[1]*nr_patches[1], h=patch[1]*nr_patches[2])
  img <- img[,,1:2]
  
  img_list <- list()
  index=1
  for(x in 1:nr_patches[1]){
    start.x <- (x*patch[1]-patch[1])+1
    end.x <- x*patch[1]
    
    for(y in 1:nr_patches[2]){
      start.y <- (y*patch[1]-patch[1])+1
      end.y <- y*patch[1]
      out <- list(image=img[start.x:end.x , start.y:end.y ,], coords=data.frame(start.x,end.x,start.y,end.y))
      
      #message(index)
      img_list[[index]] <- out
      index = index+1
    }
  }
  
  img_array <- list2array(map(img_list, ~.x$image)) #%>% aperm(., c(4,1,2,3))
  img_array %>% dim()
  
  coords_df <- map_dfr(1:length(img_list), ~img_list[[.x]]$coords %>% mutate(image=.x))
  
  
  
  return(list(image=img_array, coords=coords_df))
  
}
library(reticulate)
reticulate::use_condaenv("torch-gpu")
library(tidyverse)
np <- reticulate::import("numpy")
pydicom <- reticulate::import("pydicom")
source_python("/path/ALA Project_NYU/R_Project/Script/Run_val_fun.py")




object.list <- readRDS("object.list_ImageCluster.RDS")
colors <- readRDS("colors_cell_deconv.RDS")
names_col <- colors$colors; names(names_col) <- colors$annotation_level_4

hull <-   
  getCoordsDf(object.list[[1]]) %>% 
  ggforce::geom_mark_hull(mapping = ggplot2::aes(x = x, y = y), alpha = 1, color = "black", size = 1.5, expand = unit(3, "mm"))

#p1 <- plotSurface(object.list[[1]] , color_by = "ALA_cluster", pt_alpha = 0)+scale_color_manual(values=RColorBrewer::brewer.pal(6,"Set3"))
p2 <- 
  plotSurface(object.list[[1]], color_by = "ALA_cluster", pt_alpha = 1, display_image = F)+
  scale_color_manual(values=RColorBrewer::brewer.pal(6,"Set3"))+
  hull+
  SPATA2::ggpLayerScaleBarSI(object.list[[1]])
p1+p2


plotSurface(object.list[[1]], color_by = "Cells", alpha_by = "Cells", display_image = F, smooth = T)+
  #scale_color_manual(values=RColorBrewer::brewer.pal(6,"Set3"))+
  scale_color_gradientn(colors=c("#FFFFFF", RColorBrewer::brewer.pal(9, "Greens")),
                       limits=c(4,15), oob=scales::squish)+
  scale_alpha_continuous("mean", range = c(0.1, 1), limits=c(4,15), oob=scales::squish)+
  hull+
  SPATA2::ggpLayerScaleBarSI(object.list[[1]])


myeloid <- c("TAM_BDM_anti_infl","TAM_BDM_hypoxia_MES","TAM_BDM_INF","TAM_BDM_MHC","TAM_MG_aging_sig",
             "TAM_MG_pro_infl_I","TAM_MG_pro_infl_II","TAM_MG_prolif")


myeloid_v <- object.list[[1]] %>% joinWith(features = myeloid) %>% select({{myeloid}}) %>% rowMeans()
bc <- object.list[[1]] %>% joinWith(features = myeloid) %>% pull(barcodes)
object.list[[1]] <- object.list[[1]] %>% addFeatures(data.frame(barcodes=bc, Global_myeloid=myeloid_v))

plotSurface(object.list[[1]], color_by = "Global_myeloid", alpha_by = "Global_myeloid", display_image = F, smooth = T)+
  #scale_color_manual(values=RColorBrewer::brewer.pal(6,"Set3"))+
  scale_color_gradientn(colors=c("#FFFFFF", RColorBrewer::brewer.pal(9, "Reds")),
                        limits=c(0,0.01), oob=scales::squish)+
  scale_alpha_continuous("mean", range = c(0, 1), limits=c(0,0.01), oob=scales::squish)+
  hull+
  SPATA2::ggpLayerScaleBarSI(object.list[[1]])


#PpIX Intensity

runSegmentfromImage<- function(object, Image, spot_extension=0, multicore=T, workers = 16){
  
  #message fine correct spot size
  
  SPATA2::check_method(object)
  if (spot_extension > 0.7) stop("The spot extension will cause overlap in segmentation ")
  
  # Get spot radius
  getSpotRadius <- function(object){
    of_sample <- SPATA2::getSampleNames(object)
    coords <- SPATA2::getCoordsDf(object)
    bc_origin <- coords$barcodes
    bc_destination <- coords$barcodes
    d <- 
      tidyr::expand_grid(bc_origin, bc_destination) %>% 
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_origin = barcodes, xo = x, yo = y), by = "bc_origin") %>% 
      dplyr::left_join(x = ., y = dplyr::select(coords, bc_destination = barcodes, xd = x, yd = y), by = "bc_destination") %>% 
      dplyr::mutate(distance = base::round(base::sqrt((xd - xo)^2 + (yd - yo)^2), digits = 0)) %>% 
      dplyr::filter(distance!=0) %>% 
      dplyr::pull(distance) %>% 
      min()
    
    r = (d * c(55/100))/2
    
    return(r)
  }
  
  r <- getSpotRadius(object)
  if (!is.null(spot_extension)) {r = r + (r * spot_extension)}
  grid.plot <- SPATA2::getCoordsDf(object)
  
  dim(Image)
  if(dim(Image)[3]!=1){ Image <- Reduce(`+`, map(.x=1:dim(Image)[3], ~Image[,,.x])) }
  
  sc_dat <- reshape2::melt(Image)
  head(sc_dat)
  
  names(sc_dat)=c("x","y","var")
  
  yrange <- range(sc_dat$y)
  sc_dat$y <- (yrange[2] - sc_dat$y) + yrange[1]
  
  
  if (multicore == T) {
    base::options(future.fork.enable = TRUE)
    future::plan("multiprocess", workers = workers)
    future::supportsMulticore()
    base::options(future.globals.maxSize = 600 * 1024^2)
    message("... Run multicore ... ")
    
    #plot(sc_dat$x, sc_dat$y, pch = ".")
    segments <- furrr::future_map_dfr(.x = 1:nrow(grid.plot), 
                                      .f = function(x) {
                                        segment <- swfscMisc::circle.polygon(x = grid.plot$x[x], 
                                                                             y = grid.plot$y[x], radius = r, poly.type = "cartesian") %>% as.data.frame()
                                        
                                        
                                        nuc <- sp::point.in.polygon(pol.x = segment$x, 
                                                                    pol.y = segment$y, 
                                                                    point.x = sc_dat$x, 
                                                                    point.y = sc_dat$y)
                                        
                                        
                                        intensity <- sc_dat[nuc == 1, ]$var
                                        
                                        out <- data.frame(barcodes = grid.plot$barcodes[x],
                                                          mean=mean(intensity),
                                                          median=median(intensity),
                                                          max=max(intensity),
                                                          min=min(intensity),
                                                          sd=sd(intensity))
                                        
                                        return(out)
                                      }, .progress = T)
  }
  else {
    segments <- map_dfr(.x = 1:nrow(grid.plot), .f = function(x) {
      segment <- swfscMisc::circle.polygon(x = grid.plot$x[x], 
                                           y = grid.plot$y[x], radius = r, poly.type = "cartesian") %>% 
        as.data.frame()
      polygon(segment, border = "red")
      nuc <- sp::point.in.polygon(pol.x = segment$x, pol.y = segment$y, 
                                  point.x = sc_dat$x, point.y = sc_dat$y)
      intensity <- sc_dat[nuc == 1, ]$var
      
      out <- data.frame(barcodes = grid.plot$barcodes[x],
                        mean=mean(intensity),
                        median=median(intensity),
                        max=max(intensity),
                        min=min(intensity),
                        sd=sd(intensity))
      
      return(out)
    })
  }
  
  # Summarize on barcode level
  
  return(segments)
  
}
segments <- runSegmentfromImage(object.list[[1]], Image = getImage(object.list[[1]]))


segments$mean=segments$mean*(myeloid_v+1)
object.list[[1]] <- object.list[[1]] %>% addFeatures(segments, overwrite = T)


plotSurface(object.list[[1]], color_by = "mean", alpha_by = "mean", normalize = T, display_image = F, smooth = T)+
  #scale_color_manual(values=RColorBrewer::brewer.pal(6,"Set3"))+
  scale_color_gradientn(colors=c("#FFFFFF", RColorBrewer::brewer.pal(9, "Reds")),
                        limits=c(0.6,2), oob=scales::squish)+
  scale_alpha_continuous("mean", range = c(0, 1), limits=c(0.6,2), oob=scales::squish)+
  hull+
  SPATA2::ggpLayerScaleBarSI(object.list[[1]])






## Add cellular distributions

single_cell <- readRDS("CellMap_400_T.RDS")
library(scattermore)
ggplot()+
  scattermore::geom_scattermore(data=single_cell, mapping = aes(x,y, color=cell_type_org), pointsize = 2,pixels = c(1024, 1024))+
  scale_color_manual(values=names_col)+
  hull+
  SPATA2::ggpLayerAxesSI(object.list[[1]])+
  Seurat::NoLegend()


library(ggforce)
ggplot(single_cell[ ], aes(x, y)) +
  scattermore::geom_scattermore(data=getCoordsDf(object.list[[2]]), mapping = aes(x,y), color="black", pointsize = 15,pixels = c(1024, 1024))+
  scattermore::geom_scattermore(data=getCoordsDf(object.list[[2]]), mapping = aes(x,y,), color="white", pointsize = 10,pixels = c(1024, 1024))+
  scattermore::geom_scattermore(aes(color = cell_type_org),pointsize = 2,pixels = c(1024, 1024))+
  scale_color_manual(values=names_col)+
  coord_fixed()+
  theme_void()+
  Seurat::NoLegend()


ggplot(single_cell[1:100, ], aes(x, y, group = -1L)) +
  geom_voronoi_tile(aes(fill = cell_type_org), max.radius = 10)+
  scale_fill_manual(values=names_col)+
  coord_fixed()

ggplot(single_cell[1:100, ], aes(x, y, group = -1L)) +
  geom_voronoi_tile(aes(fill = celltypes), max.radius = 10)+
  scale_fill_manual(values=names_col)+
  theme_void()+
  coord_fixed()+
  Seurat::NoLegend()


## 2nd request the Barplots of pattern as ratios


cluster_names=data.frame(C1="Artifact",
                         C2="low_signal_Cell_like",
                         C3="Cell_Diffuse_like",
                         C4="Cell_like",
                         C5="Diffuse_like",
                         C6="Fiber_like") %>% t() %>% as.data.frame()
names(cluster_names)="ALA_cluster_names"
cluster_names$ALA_cluster=as.factor(1:6)



# Revision on the Pattern -------------------------------------------------

object <- object.list[[1]]
org_path <- "/outs/spatial/ALA.png"
file.exists(org_path)

getImagePrediction <- function(object, org_path){
  
  object <- SPATA2::setImageDirHighres(object, org_path)
  object <- SPATA2::loadImageHighres(object)
  #object %>% plotSurface()
  img.list <- 
    SPATA2::getImageSectionsByBarcode(object,
                                      barcodes = SPATA2::getBarcodes(object) )
  
  
  # Scale the images and return the 160 by 160 2 Channels
  img.list_A <- 
    map(.x=1:length(img.list), .f=function(i){
      img <- img.list[[i]]$image[,,1:2]
      img[,,2] <- 0
      img <- img %>% EBImage::resize(., 160,160)
      return(img)
      }, .progress = T)
  
  train_X <- list2array(img.list_A)
  dim(train_X)
  #train_X <- train_X %>% aperm(., c(4,1,2,3))
  train_X <- scales::rescale(train_X, c(0,255))
  
  lastlayer <- runValidation(train_X) %>% as.data.frame()
  names(lastlayer) <- c("AF", "Diff_low", "Diff_bright", "Cellular", "Fiber_like")
  labels <- map(1:nrow(lastlayer), ~names(lastlayer)[which.max(lastlayer[.x,])]) %>% unlist()
  
  lastlayer$barcodes <- SPATA2::getBarcodes(object)
  lastlayer$lables <- labels
  
  object <- addFeatures(object, lastlayer)
  
  #object %>% plotSurfaceInteractive()
  
  return(object)
  
}

object <- getImagePrediction(object, org_path)





