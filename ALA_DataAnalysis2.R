

# All samples with ALA an st

folder=paste0("~/Desktop/NIO/GBM/STData/", dir(pattern="+ALA"))
sample=dir(pattern="+ALA")
outs <- paste0(folder, "/outs")

images <- paste0(folder, "/outs/spatial")
SRH <- paste0(images, "/SRH.png")
HE <- paste0(images, "/HE.png")
ALA <- paste0(images, "/ALA.png")
ALA_black <- paste0(images, "/ALA_black.png")

i=6

object <- SPATA2::initiateSpataObject_10X(folder[i], sample_name = sample[i])
img.dim <- SPATA2::getImage(object) %>% dim()

object <- SPATA2::exchangeImage(object, HE[i], resize = img.dim[1:2]*3)
object <- SPATA2::createSegmentation(object)

bc_select <- SPATA2::getFeatureDf(object) %>% dplyr::select(barcodes, histology) %>% filter(histology=="Sample") %>% pull(barcodes)
object_segment <- SPATA2::subsetByBarcodes(object, barcodes = bc_select)

# Read in Nuc Positions

img <- EBImage::readImage(ALA_black[i])
object_segment <- SPATA2::exchangeImage(object_segment, img, resize = dim(img)[1:2]*0.25)


all_pos <- read.csv("~/Desktop/NIO/GBM/Nuc_Classifier/Images/Result_IdentifyPrimaryObjects.csv")
#all_pos <- read.csv("~/Desktop/NIO/GBM/STData/NIO_410+ALA/outs/spatial/Result_IdentifyPrimaryObjects.csv")

out_pos <- 
  all_pos %>% filter(ImageNumber==5) %>% mutate(
  Cell=paste0("Cell_",ObjectNumber),
  x=Location_Center_X,
  y= SPATA2::getImageRange(object_segment)$y[2] - Location_Center_Y + SPATA2::getImageRange(object_segment)$y[1]) %>% 
  dplyr::select(Cell,x,y)


hull <-   
  getCoordsDf(object_segment) %>% 
  ggforce::geom_mark_hull(mapping = ggplot2::aes(x = x, y = y), alpha = 1, color = "black", size = 1, expand = unit(2, "mm"))
ggplot(out_pos, mapping=aes(x,y))+
  geom_point(color="black", size=0.1)+
  theme_classic()+
  hull

object_segment@spatial[[getSampleName(object_segment)]]$Cell_coords <- out_pos


saveRDS(object_segment, paste0(folder[i], "/spataObject_segmented.RDS"))
saveRDS(list(EBImage::readImage(HE[i]),
             EBImage::readImage(SRH[i]),
             EBImage::readImage(ALA[i]),
             EBImage::readImage(ALA_black[i])),
        paste0(folder[i], "/Images.RDS"))




obj <- readRDS(paste0(folder[2], "/spataObject_segmented_410_T.RDS"))



#Images:



color_by = "CHI3L1"
hull <-   
  getCoordsDf(object_segment) %>% 
  ggforce::geom_mark_hull(mapping = ggplot2::aes(x = x, y = y), alpha = 1, color = "white",concavity=1.5, size = 1, expand = unit(2, "mm"))
p1=plotSurface(object_segment, 
            smooth = T,
            color_by = color_by, 
            #alpha_by = color_by, 
            pt_alpha = 0,
            display_image = T, 
            pt_clrsp="Reds",
            limits=c(0.4,1),
            oob=scales::squish)+
  ggdark::dark_theme_bw()+ 
  hull+
  Seurat::NoAxes()+
  Seurat::NoLegend()


color_by = "AIF1"
hull <-   
  getCoordsDf(object_segment) %>% 
  ggforce::geom_mark_hull(mapping = ggplot2::aes(x = x, y = y), alpha = 1, color = "black",concavity=1.5, size = 1, expand = unit(1, "mm"))
p2=plotSurface(object_segment, 
            smooth = T,
            color_by = color_by, 
            alpha_by = color_by,  
            display_image = F, 
            pt_clrsp="Reds",
            limits=c(0.3,0.8),
            oob=scales::squish)+
  hull+
  SPATA2::ggpLayerScaleBarSI(object_segment, sgmt_size=0.5, text_size=4)+
  Seurat::NoAxes()

p1+p2












