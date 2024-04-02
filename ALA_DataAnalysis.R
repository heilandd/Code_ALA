

# Load data
list_data <- readRDS("~/Desktop/NIO/ALA/Metabolomic/list_data.RDS")


# Start by plot ALA 

PpIX <- 722.89 #or# 563.2653
ALA <- 130.0509
Porphobilinogen <- 225.0881



obj <- readRDS(list_data$dir[2])
obj <- obj %>% setActiveExpressionMatrix(mtr_name="MALDI")
#plotSurfaceInteractive(obj)

getGenes(obj)[str_detect(getGenes(obj), pattern = "563")]
PpIX <- c("722.024007907827","722.898908733854")
#PpIX <- c("563.294938237066","563.977501699448")
ALA <- c("130.440238737979","130.756547930667")
Porphobilinogen <- c("225.219672451152","225.492578719235")

#col <- colorRampPalette((RColorBrewer::brewer.pal(9,"Greens")))
col <- colorRampPalette((viridis::plasma(50)))
plotSurface2(joinWith(obj, genes=PpIX, average_genes=T, smooth = T, smooth_span = 0.1), 
             color_by = "mean_genes", 
             pt_alpha = c(joinWith(obj, genes=PpIX, average_genes=T) %>% pull(mean_genes)+0.3),
             image = getImage(obj)
             )+scale_colour_gradientn(colours = col(50), oob = scales::squish, limit=c(0.2, 0.6))


plotSurface2(joinWith(obj, genes=ALA, average_genes=T, smooth = T, smooth_span = 0.1), 
             color_by = "mean_genes", 
             pt_alpha = c(joinWith(obj, genes=ALA, average_genes=T) %>% pull(mean_genes)),
             image = getImage(obj)
)+scale_colour_gradientn(colours = col(50), oob = scales::squish, limit=c(0., 0.8))


plotSurface2(joinWith(obj, genes=Porphobilinogen, average_genes=T, smooth = T, smooth_span = 0.1), 
             color_by = "mean_genes", 
             pt_alpha = c(joinWith(obj, genes=Porphobilinogen, average_genes=T) %>% pull(mean_genes)),
             image = getImage(obj)
)+scale_colour_gradientn(colours = col(50), oob = scales::squish, limit=c(0.2, 0.8))




getFeatureNames(obj)[str_detect(getFeatureNames(obj), pattern = "TAM")]

cell <- "TAM_MG_aging_sig"
col <- colorRampPalette((RColorBrewer::brewer.pal(9,"Greens")))
plotSurface(obj,normalize = F, color_by = cell, display_image = T, alpha_by = cell)+
  scale_colour_gradientn(colours = col(50), oob = scales::squish, limit=c(0, 0.15))


# Co localisation
colors <- readRDS("~/Desktop/ImmunoSpatial/Paper/colors_cell_deconv.RDS")
cor_map <- map(.x=1:6, .f=function(i){
  
  obj <- readRDS(list_data$dir[i])
  
  PpIX <- c("722.024007907827","722.898908733854")
  #PpIX <- c("563.294938237066","563.977501699448")
  ALA <- c("130.440238737979","130.756547930667")
  Porphobilinogen <- c("225.219672451152","225.492578719235")
  modules <- getGeneSets(obj, index="Module")
  
  obj_genes <- SPATA2::setActiveExpressionMatrix(obj, "scaled")
  
  df_list <- list(
  joinWith(obj,genes=PpIX, average_genes =T) %>% select(barcodes, mean_genes) %>% dplyr::rename("PpIX":=mean_genes),
  joinWith(obj,genes=ALA, average_genes =T) %>% select(barcodes, mean_genes) %>% dplyr::rename("5_ALA":=mean_genes),
  joinWith(obj,genes=Porphobilinogen, average_genes =T) %>% select(barcodes, mean_genes) %>% dplyr::rename("Porphobilinogen":=mean_genes),
  joinWith(obj, features = colors$annotation_level_4) %>% select(barcodes, colors$annotation_level_4),
  joinWith(obj_genes, genes = c("AIF1", "CD163", "CD68")) %>% select(-x,-y,-row,-col,-sample),
  joinWith(obj_genes, genes = c("CPOX", "HMBS", "FECH")) %>% select(-x,-y,-row,-col,-sample),
  joinWith(obj_genes, gene_sets = modules ) %>% select(barcodes, modules)
  )
  
  df <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  rownames(df) <- df$barcodes
  df$barcodes <- NULL
  
  cor <- cor(df)
  
  return(cor)
  
})
cor_map <-  Reduce(`+`, cor_map) / 6

meta <- c("CPOX", "HMBS", "FECH")
macro <- c("AIF1", "CD163", "CD68")

modules <- getGeneSets(obj, index="Module")
corrplot::corrplot(cor_map[c("PpIX","Porphobilinogen"),modules], is.corr = F, 
                   col=colorRampPalette(rev(RColorBrewer::brewer.pal(9,"BrBG")))(50))


cor_plot <- cor_map[meta,modules]
cor_plot[cor_plot<0]=0
corrplot::corrplot(cor_plot , is.corr = F, 
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Blues")))(50))


cor_plot <- cor_map[meta,colors$annotation_level_4]
cor_plot[cor_plot<0]=0
corrplot::corrplot(cor_plot , is.corr = F, 
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Reds")))(50))


cor_plot <- cor_map[c("PpIX","Porphobilinogen"),colors$annotation_level_4]
cor_plot[cor_plot<0]=0
corrplot::corrplot(cor_plot , is.corr = F, 
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"Blues")))(50))



corrplot::corrplot(cor_map[meta,modules], is.corr = F, 
                   col=colorRampPalette((RColorBrewer::brewer.pal(9,"BrBG")))(50))


corrplot::corrplot(cor_map[meta,macro] %>% scales::rescale(., c(0,1)), is.corr = F, 
                   col=colorRampPalette(rev(RColorBrewer::brewer.pal(9,"BrBG")))(50))


colors$


obj <- readRDS(list_data$dir[6])
obj <- obj %>% setActiveExpressionMatrix(mtr_name="denoised")

obj <- runAutoencoderDenoising(obj, activation = "relu", bottleneck = 32)


cell <- "CD68"
col <- colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Spectral")))
plotSurface(obj,normalize = T, color_by = cell, display_image = T, alpha_by = cell, smooth=F)+
  scale_colour_gradientn(colours = col(50), oob = scales::squish, limit=c(0.4, 0.9))



obj <- readRDS(list_data$dir[2])
obj <- obj %>% setActiveExpressionMatrix(mtr_name="MALDI")

bc1 <- joinWith(obj, genes=PpIX, average_genes=T) %>% filter(mean_genes>0.9) %>% pull(barcodes)
bc2 <- joinWith(obj, genes=PpIX, average_genes=T) %>% filter(mean_genes==0) %>% pull(barcodes)

MALDI <- getExpressionMatrix(obj)
plot_data1 <- MALDI[,bc1] %>% rowMeans() %>% as.data.frame()
plot_data1 <- plot_data1 %>% rownames_to_column("MHZ") %>% mutate(MHZ=as.numeric(MHZ))
names(plot_data1) <- c("MHZ", "Peak")

plot_data2 <- MALDI[,bc2] %>% rowMeans() %>% as.data.frame()
plot_data2 <- plot_data2 %>% rownames_to_column("MHZ") %>% mutate(MHZ=as.numeric(MHZ))
names(plot_data2) <- c("MHZ", "Peak")

ggplot()+theme_classic()+
  geom_point(data=plot_data1, aes(x=MHZ, y=Peak), color="red")+
  geom_point(data=plot_data2, aes(x=MHZ, y=Peak), color="green")+
  ylim(0,80000)+
  xlim(719,723)+
  #scale_x_continuous(breaks=seq(0,max(plot_data$MHZ), length.out=25))+
  theme(axis.text.x = element_text(angle = 70))







