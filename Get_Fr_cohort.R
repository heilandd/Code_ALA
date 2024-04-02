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
source_python("/Users/path/Desktop/ALA Project_NYU/R_Project/Script/Run_val_fun.py")



root <- "/Users/path/Desktop/Projects/path/Data/All_raw"
all_files <- dir("~/Desktop/Projects/path/Data/All_raw", pattern = "ALA")

patients <- paste0("Patient_", all_files %>% str_split(., "_") %>% map(~.x[[2]]) %>% unlist())
patients %>% unique()

subImage2=paste0(all_files %>% str_split(., "_") %>% map(~.x[[4]]) %>% unlist())
all_files[which(subImage2=="Sample")]
subImage2[which(subImage2=="Sample")] <- c("1A")

meta_data <- data.frame(patient=patients, path=all_files, subfile=subImage2)




# Create a NYU Structure  -------------------------------------------------


create_NIO_folder_root <- ("~/Desktop/All_Institutions/UKF")
library(progress)
pb <- progress_bar$new(
  format = "  Run Analysis [:bar] :percent in :elapsed",
  total = nrow(meta_data), clear = FALSE, width= 60)
for(i in 1:nrow(meta_data)){
  pb$tick()
  setwd(create_NIO_folder_root)
  if(!dir.exists(meta_data$patient[i])){dir.create(meta_data$patient[i])}
  subfolder <- meta_data %>% filter(patient==meta_data$patient[i])
  
  for(subfile in unique(subfolder$subfile)){
    subroot <- paste0(create_NIO_folder_root,"/", meta_data$patient[i])
    setwd(subroot)
    if(!dir.exists(subfile)){dir.create(subfile)}
    
    list.files <- c(subfolder[subfolder$subfile==subfile, ]$path)
    file.copy(paste0(root, "/",list.files), paste0(subroot, "/",subfile,"/",list.files))
    
    file.rename(paste0(subroot, "/",subfile,"/",list.files), paste0(subroot, "/",subfile,"/",paste0("img1_", 1:length(list.files),".dcm")))
    
  }
  
  setwd(create_NIO_folder_root)
}





# run FR Images -----------------------------------------------------------



library(progress)
pb <- progress_bar$new(
  format = "  Run Analysis [:bar] :percent in :elapsed",
  total = length(all_files), clear = FALSE, width= 60)
for( i in 1:length(all_files)){
  pb$tick()
  
  image = pydicom$dcmread(paste0(root,"/",all_files[i]))
  load_ala_image <-np$asarray(image$pixel_array) %>% EBImage::Image(., colormode="Color")
  patches <- getPatches(load_ala_image)
  img_array <- patches$image
  save_file=paste0("/Users/path/Desktop/ALA Project_NYU/Data_all_Patches/Fr_cohort/Image_",i,".npy")
  np$save(save_file, np$asarray(img_array))

  
}



## Run Prediction
# with out save in between
prediction <- map_dfr(.x=1:length(all_files), .f=function(i){
  
  image = pydicom$dcmread(paste0(root,"/",all_files[i]))
  load_ala_image <-np$asarray(image$pixel_array) %>% EBImage::Image(., colormode="Color")
  patches <- getPatches(load_ala_image)
  img_array <- patches$image
  #save_file=paste0("/Users/path/Desktop/ALA Project_NYU/Data_all_Patches/Fr_cohort/Image_",i,".npy")
  #np$save(save_file, np$asarray(img_array))
  
  lastlayer <- runValidation(img_array) %>% as.data.frame()
  
  labels <- map(1:nrow(lastlayer), ~which.max(lastlayer[.x,])) %>% unlist()
  pat <- meta_data$patient[i]
  
  lastlayer$label <- as.character(labels)
  lastlayer$case <- pat
  lastlayer$subImage <- meta_data$subfile[i]
  
  return(lastlayer)
  
})

prediction %>% dim()
saveRDS(prediction, "~/Desktop/ALA Project_NYU/Data_all_Patches/prediction.RDS")


## Annotate regions

clinical_data_Fr <- data.frame(institution="UKF", case=prediction$case, patch=prediction$subImage, label=prediction$label)
clinical_data$X <- NULL
clinical_data <- rbind(clinical_data, clinical_data_Fr)

clinical_data$institution %>% unique()

saveRDS(clinical_data, "~/Desktop/ALA Project_NYU/Data_all_Patches/clinicalData.RDS")



## Add clinical data

NIO_GBM<- read.csv("~/Desktop/Projects/path/Data/NIO_GBM_datasheet_v5_Sample.csv", sep=";")

data <- data.frame(path=NIO_GBM$X5ALA.Image, region=NIO_GBM$X , History=NIO_GBM$History, sex=NIO_GBM$sex) %>% filter(path!="")


meta_data$path=meta_data$path %>% str_remove(., ".dcm")
meta_data <- left_join(meta_data, data)


clinical_data_Fr$case_sub <- paste0(clinical_data_Fr$case,"_", clinical_data_Fr$patch )
meta_data$case_sub <- paste0(meta_data$patient,"_", meta_data$subfile )

clinical_data_Fr <- left_join(clinical_data_Fr, meta_data)

clinical_data <- clinical_data_Fr

dim(clinical_data)
clinical_data$label <- as.character(clinical_data$label)
clinical_data$case <- as.character(clinical_data$case)

clinical_data <- 
clinical_data %>% mutate(region=case_when(
  region %in% c("Subdural_membrane", "Membrane") ~ "Infiltration Meninges",
  region %in% c("Cavitywall") ~ "Cavitywall",
  region %in% c("Tumor_Contrast", "Tumor_Core", "Tumor", "Ctx", "Tumor_Contrast_Periventricular") ~ "Tumor CE-enriched",
  region %in% c("WM") ~ "Outside FLAIR-hyper",
  region %in% c("Infiltration", "Tumor_Border_Ctx") ~ "FLAIR-hyper",
  TRUE~"NA"
  
))

out <- clinical_data %>% group_by(region) %>% count(label) %>%  mutate(per= prop.table(n) * 100) %>% ungroup()
#out$case <- as.character(out$case)


order <- 
  out %>% 
  select(region,label,per) %>% 
  pivot_wider(names_from =  label, values_from = per, values_fill=0) %>%  
  arrange(`1`,`2`) %>% 
  pull(region) 

out$region <- factor(out$region, level = order)

out$region %>% unique()




ggplot(out, aes(y=per, x=region, fill=label))+
  geom_col(position = "fill")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(5, "Set2")[c(1,3,4,2,5)] )+
  theme_classic()


ggplot(out %>% filter(label %in% c("4", "5")), aes(y=per, x=case, fill=label))+
  geom_col()+
  scale_fill_manual(values=RColorBrewer::brewer.pal(5, "Set2")[c(1,4)] )+
  theme_classic()


ggplot(out, aes(y=1, x=case, fill=institution ))+
  geom_col(position = "fill")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Greys") )+
  theme_classic()







