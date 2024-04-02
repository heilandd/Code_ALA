
library(reticulate)
reticulate::use_condaenv("torch-gpu")
library(tidyverse)

np <- reticulate::import("numpy")
last_layer <- np$load("/path_to_file/ALA Project_NYU/Data_all_Patches/last_layer.npy")
labels <- map(1:nrow(last_layer), ~which.max(last_layer[.x,])) %>% unlist()


#get clinical data 

clinical_data <- read.csv("~/Downloads/all_ala_patches.csv", sep=";")
dim(clinical_data)
clinical_data <- readRDS("path_to_file/Data_all_Patches/clinicalData.RDS")
prediction <- readRDS("path_to_file/Data_all_Patches/prediction.RDS")


last_layer_all <- rbind(last_layer, prediction[, 1:5])

dim(last_layer_all)
dim(clinical_data)

## Get the values of each class
clinical_data <- cbind(clinical_data, last_layer_all)
names(clinical_data)[c(5,7,6,8,9)] <- c("AF", "Diffuse_low", "Diffuse_bright", "Cellular", "Fiber_like")

write.csv(clinical_data, "path_to_file/Data_all_Patches/clinical_data_all_new.csv", row.names = F)


## Correlate classes
correlation <- cor(clinical_data[,c(5,7,6,8,9)])
corrplot::corrplot(correlation, type="lower", col=RColorBrewer::brewer.pal(9, "PiYG"))
hclust(as.dist(correlation)) %>% plot()



## Create new lables in two Categories

#AF, Diffuse Low Cellular pattern

mat_val <- clinical_data[,c(5,7,6,8,9)]
names(mat_val) <- c("AF", "Diffuse_low", "Diffuse_bright", "Cellular", "Fiber_like")



## Major classes
vector_x <- apply(mat_val[,1:3], 1, which.max)  
values <- map(.x=1:nrow(mat_val), .f=function(.x){
  mat_val[.x, vector_x[.x]]
} ) %>% unlist()
data_major_classes <- 
  data.frame(class=vector_x, values=values) %>% mutate(class=case_when(class==1~"AF", 
                                                                     class==2~"Diffuse_low", 
                                                                     class==3~"Diffuse_bright"))

data_major_classes <- cbind(clinical_data,data_major_classes)

## Get Subclasses
vector_x <- mat_val[,4:5] %>% mutate(class=case_when(Cellular>c(-5) & Fiber_like<c(-5) ~ "Cellular", 
                                                     Cellular<c(-5) & Fiber_like>c(-10)~"Fiber_like",
                                                     Cellular>c(-5) & Fiber_like>c(-10)~"Mixed",
                                                     Cellular<c(-5) & Fiber_like<c(-10)~"Neither")) %>% 
  mutate(ID=paste0("ID", 1:nrow(.)))

df_values <- do.call("rbind", list(
  vector_x %>% filter(class=="Cellular") %>% select(ID, Cellular) %>% rename("values":=Cellular),
  vector_x %>% filter(class=="Fiber_like") %>% select(ID, Fiber_like) %>% rename("values":=Fiber_like),
  vector_x %>% filter(class=="Mixed") %>% mutate(values=Cellular+Fiber_like) %>% select(ID, values),
  vector_x %>% filter(class=="Neither") %>% mutate(values=Cellular+Fiber_like) %>% select(ID, values)))
vector_x <- left_join(vector_x, df_values)


data_sub_classes <- 
  data.frame(class=vector_x$class, values=vector_x$values)

data_sub_classes <- cbind(clinical_data,data_sub_classes)

dim(data_sub_classes)


# Get the upper class Barplot
data_major_classes
out <- data_major_classes %>% group_by(case) %>% count(class) %>%  mutate(per= prop.table(n) * 100) %>% ungroup()
inst <- data_major_classes[!duplicated(data_major_classes$case), ] %>% select(case, institution)
inst$case <- as.character(inst$case)

order <- 
  out %>% 
  select(case,class,per) %>% 
  pivot_wider(names_from =  class, values_from = per, values_fill=0) %>%  
  arrange(AF,Diffuse_low) %>% 
  pull(case) 
out <- left_join(out, inst)
out$case <- factor(out$case, level = order)

ggplot(out, aes(y=per, x=case, fill=class))+
  geom_col(position = "fill")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(5, "Set2")[c(1,3,4,2,5)] )+
  theme_classic()


out <- data_sub_classes %>% group_by(case) %>% count(class) %>%  mutate(per= prop.table(n) * 100) %>% ungroup()
out$case <- factor(out$case, level = order)
ggplot(out, aes(y=per, x=case, fill=class))+
  geom_col(position = "fill")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(4, "Set3"))+
  theme_classic()

pca <- irlba::prcomp_irlba(mat_val)
summary(pca)
pca_plot <- data.frame(class=names(mat_val), PC1=pca$rotation[,2])
ggplot(pca_plot, aes(x=class, y=PC1))+geom_col()


#### Old barplots .....

dim(clinical_data)
clinical_data$label <- as.character(labels)
clinical_data$case <- as.character(clinical_data$case)



out <- clinical_data %>% group_by(case) %>% count(label) %>%  mutate(per= prop.table(n) * 100) %>% ungroup()
out$case <- as.character(out$case)

inst <- clinical_data[!duplicated(clinical_data$case), ] %>% select(case, institution)
inst$case <- as.character(inst$case)

inst$institution %>% table()

order <- 
  out %>% 
  select(case,label,per) %>% 
  pivot_wider(names_from =  label, values_from = per, values_fill=0) %>%  
  arrange(`1`,`2`) %>% 
  pull(case) 



out$case <- factor(out$case, level = order)

out <- left_join(out, inst)
out$case <- factor(out$case, level = order)

ggplot(out, aes(y=per, x=case, fill=label))+
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





