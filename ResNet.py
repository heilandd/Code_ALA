import torch
from torch.utils.data import Dataset
class NumpyArrayDataset(Dataset):
    def __init__(self, npy_array, transform=None):
        self.npy_array = npy_array
        self.transform = transform

    def __len__(self):
        return len(self.npy_array)

    def __getitem__(self, idx):
        image = self.npy_array[idx]
        if self.transform:
            image = self.transform(image)
        return image

## import the NYU npy


import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torchvision.transforms as transforms
import torchvision.models as models
from torch.utils.data import DataLoader


# Load data:
import numpy as np
img = np.load("/home/dhh/Desktop/NYU/patches.npy")
lables = np.load("/home/dhh/Desktop/NYU/labels.npy")
img = np.transpose(img, (3, 0, 1, 2))

indices = np.random.choice(len(img), 20000, replace=False)
img = img[indices]
lables = lables[indices]


## Split data
import numpy as np
from sklearn.model_selection import train_test_split

# Split the data and labels into training and validation sets with equal label proportions
train_data, val_data, train_labels, val_labels = train_test_split(img, lables, test_size=0.2, stratify=lables)





import matplotlib.pyplot as plt
from PIL import Image

img
plt.imshow(img[280,:,:,:])
plt.axis('off')
plt.show()



# Define hyperparameters
batch_size = 32
learning_rate = 0.001
num_epochs = 2

# Define the ResNet50 model
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
resnet50 = models.resnet50(pretrained=False)
num_classes = 5  
resnet50.fc = nn.Linear(in_features=2048, out_features=num_classes)
resnet50.to(device)

# Define a transformation to preprocess and augment the input images
transform_train = transforms.Compose([
    transforms.RandomResizedCrop(224),
    transforms.RandomHorizontalFlip(),
    transforms.ColorJitter(brightness=0.4, contrast=0.4, saturation=0.4, hue=0.1),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])
transform_val = transforms.Compose([
    transforms.Resize(256),
    transforms.CenterCrop(224),
    transforms.ToTensor(),
    transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225])
])

train_datasettorch.from_numpy(train_labels[i:i+batch_size])

train_dataset = NumpyArrayDataset(train_data, transform=transform_train)
train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

val_dataset = NumpyArrayDataset(val_data, transform=transform_val)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

# Define the loss function and optimizer
criterion = nn.CrossEntropyLoss()
optimizer = optim.Adam(resnet50.parameters(), lr=learning_rate)

# Train the model
for epoch in range(num_epochs):
    # Training loop
    resnet50.train()
    train_loss = 0.0
    train_correct = 0
    for inputs in train_loader:
        inputs = inputs.to(device)
        labels = torch.from_numpy(train_labels[i:i+batch_size]).to(device)  
        optimizer.zero_grad()
        outputs = resnet50(inputs)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()

        train_loss += loss.item() * inputs.size(0)
        _, preds = torch.max(outputs, 1)
        train_correct += torch.sum(preds == labels.data)

    train_loss /= len(train_loader.dataset)
    train_acc = train_correct.double() / len(train_loader.dataset)

    # Validation loop
    resnet50.eval()
    val_loss = 0.0
    val_correct = 0
    with torch.no_grad():
        for inputs in val_loader:
            inputs = inputs.to(device)
            labels = torch.from_numpy(val_labels[i:i+batch_size]).to(device)  # extract label vector from file
            outputs = resnet50(inputs)
            loss = criterion(outputs, labels)

            val_loss




