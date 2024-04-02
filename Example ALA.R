reticulate::use_condaenv("torch-gpu")
dyn.load('/opt/R/arm64/gfortran/lib/libgcc_s.2.dylib')
dyn.load('/opt/R/arm64/gfortran/lib/libgfortran.5.dylib')
library(SPATA2)
library(Seurat)
library(EBImage)
library(tidyverse)
library(reticulate)
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
np <- reticulate::import("numpy")
pydicom <- reticulate::import("pydicom")
source_python("/path/R_Project/Script/Run_val_fun.py")

##Examples
#pat_0037_Sample_5B_img9_2_ALA
#pat_0101_Sample_1A_img1_2_ALA
#pat_0142_Sample_1F_img6_2_ALA
#pat_0147_Sample_3A_img5_2_ALA

root <- "/path/Data/All_raw"
all_files <- dir("~/Desktop/All_raw", pattern = "ALA")
img_1 <- "/path/All_raw/pat_0147_Sample_3A_img5_2_ALA.dcm"
img_1 <- "/path/All_raw/pat_0142_Sample_1F_img6_2_ALA.dcm"


image = pydicom$dcmread(img_1)
load_ala_image <-np$asarray(image$pixel_array) %>% EBImage::Image(., colormode="Color")
patches <- getPatches(load_ala_image)
img_array <- patches$image
lastlayer <- runValidation(img_array) %>% as.data.frame()
labels <- map(1:nrow(lastlayer), ~which.max(lastlayer[.x,])) %>% unlist()

coords <- patches$coords
coords$label <- as.character(labels)
coords$x_coord <- coords$start.x+(c(coords$end.x-coords$start.x)/2)
coords$y_coord <- coords$start.y+(c(coords$end.y-coords$start.y)/2)


## Heatmap for lable 4

coords$Lable4 <- lastlayer$V5

dyn.load('/opt/R/arm64/gfortran/lib/libgcc_s.2.dylib')
dyn.load('/opt/R/arm64/gfortran/lib/libgfortran.5.dylib')
library(SPATA2)
smooth <- hlpr_smooth(var_name="Lable_4",variable=coords$Lable4, coords_df =  data.frame(x=coords$start.x,
                                                                               y=coords$start.y,
                                                                               Lable_4 = coords$Lable4),
            subset="Lable_4", smooth_span=0.5)

coords$Lable4 <- smooth

library(scales)
ggplot(coords, mapping = aes(xmin=start.x, xmax=end.x, ymin=start.y, ymax=end.y, fill=Lable4, alpha=Lable4) )+
  scale_fill_gradientn(colors=c("#FFFFFF", RColorBrewer::brewer.pal(9, "Reds")),
                       limits=c(-10,1), oob=scales::squish)+
  geom_rect()+
  theme_void()


#1:AF, 2:diff:bright, 3:Diff_dim, 4:Cellular, 5:Fiber

ggplot(coords, mapping = aes(xmin=start.x, xmax=end.x+2, ymin=start.y, ymax=end.y+2, fill=label) )+
  scale_fill_manual(values=RColorBrewer::brewer.pal(5, "Set2")[c(1,3,4,2,5)] )+
  geom_rect()+
  theme_void()

## plot Image
img_col <- np$asarray(image$pixel_array) %>% scales::rescale(c(0,1)) %>% EBImage::Image(., colormode="Color")
plot(img_col*1.4)
img_grey <- np$asarray(image$pixel_array)[,,1] %>% scales::rescale(c(0,1)) %>% EBImage::Image() 
plot(img_grey+0.3)



## Get the UMAP embedding
source_python("/path/R_Project/Script/RunUMAP_emb.py")
umap_new <- getUMAP(lastlayer) %>% as.data.frame()
umap_new$class <- as.character(labels)

umap_all <- readRDS("~/Desktop/ALA Project_NYU/NYU/umap.RDS")
labels_all <- map(1:nrow(last_layer), ~which.max(last_layer[.x,])) %>% unlist()
umap_all$labels <- as.character(labels_all)

ggplot()+
  scattermore::geom_scattermore(umap_all, mapping=aes(x=V1, y=V2), pointsize = 3, alpha=1, color="black")+
  scattermore::geom_scattermore(umap_all, mapping=aes(x=V1, y=V2), pointsize = 2, alpha=1, color="white")+
  scattermore::geom_scattermore(umap_new, mapping=aes(x=V1, y=V2, color=class), pointsize = 2 )+
  scale_color_manual(values=RColorBrewer::brewer.pal(5, "Set2") )+
  theme_classic()+
  coord_equal()















