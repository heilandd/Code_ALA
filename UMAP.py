import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torchvision.transforms as transforms
import torchvision.models as models

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
resnet18 = torch.load('/home/dhh/Desktop/NYU/resnet18_model.pth').to(device)

img = np.load("/home/dhh/Desktop/NYU/patches.npy")
num_images = 40742

last_layer_outputs = np.load("/home/dhh/Desktop/NYU/last_layer.npy")

class1 = r.images_class1_top
class1_img = img[:,:,:,class1]
class1_img = class1_img.transpose(3,0,1,2)
class1_img.shape


# Compute the outputs of the last layer for all images
resnet18.eval()
last_layer_outputs = np.empty((num_images, 5))

for i in range(num_images):
    image = img[:,:,:,i].transpose(2, 0, 1)
    image = torch.from_numpy(image).float().to(device)
    
    output = resnet18(image.unsqueeze(0)).detach().cpu().numpy()
    last_layer_outputs[i] = output
np.save("/home/dhh/Desktop/NYU/last_layer.npy", last_layer_outputs)

import matplotlib.pyplot as plt
import umap
umap_transformer = umap.UMAP(n_neighbors=15, min_dist=0.1)
umap_representation = umap_transformer.fit_transform(last_layer_outputs)


