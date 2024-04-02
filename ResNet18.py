
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import torchvision.transforms as transforms
import torchvision.models as models
from torch.utils.data import DataLoader


train_data = np.load("/home/dhh/Desktop/NYU/train_data.npy")
train_labels = np.load("/home/dhh/Desktop/NYU/train_labels.npy")
val_data = np.load("/home/dhh/Desktop/NYU/val_data.npy")
val_labels = np.load("/home/dhh/Desktop/NYU/val_labels.npy")


import torch
import torch.nn as nn
import torchvision.models as models

# Define the ResNet18 model with 2-channel input
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
resnet18 = models.resnet18(pretrained=False)
num_classes = 5  # replace with the number of classes in your dataset

resnet18.conv1 = nn.Conv2d(in_channels=2, out_channels=64, kernel_size=7, stride=2, padding=3, bias=False)
resnet18.fc = nn.Linear(in_features=512, out_features=num_classes)
resnet18.to(device)

criterion = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(resnet18.parameters())

# Load the data and labels from the NumPy files
data = torch.from_numpy(train_data).float()
data = data.permute(0, 3, 1, 2)  
labels = torch.from_numpy(train_labels).long()

# Train the model
num_epochs = 100
batch_size = 16
num_batches = len(data) // batch_size


for epoch in range(num_epochs):
    resnet18.train()
    epoch_loss = 0.0
    epoch_accuracy = 0.0

    # Shuffle the data and labels at the start of each epoch
    indices = torch.randperm(len(data))
    shuffled_data = data[indices]
    shuffled_labels = labels[indices]

    for i in range(num_batches):
        # Extract a batch of data and labels
        batch_data = shuffled_data[i*batch_size:(i+1)*batch_size].to(device)
        batch_labels = shuffled_labels[i*batch_size:(i+1)*batch_size].to(device)


        # Forward pass
        optimizer.zero_grad()
        outputs = resnet18(batch_data)
        loss = criterion(outputs, batch_labels)
        epoch_loss += loss.item()

        # Backward pass and optimization
        loss.backward()
        optimizer.step()

        # Calculate accuracy
        _, predictions = torch.max(outputs, 1)
        epoch_accuracy += torch.sum(predictions == batch_labels)

    epoch_loss /= num_batches
    epoch_accuracy /= len(data)

    print('Epoch [{}/{}], Loss: {:.4f}, Accuracy: {:.2%}'.format(epoch+1, num_epochs, epoch_loss, epoch_accuracy))


torch.save(resnet18, '/home/dhh/Desktop/NYU/resnet18_model.pth')