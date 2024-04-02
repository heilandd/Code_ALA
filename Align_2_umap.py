import numpy as np
import umap

test_data = np.load("/path/Data_all_Patches/last_layer.npy")

np.isnan(test_data.any())

test_embedding = fitter.transform(test_data)

np.save("/path/Data_all_Patches/umap.npy", test_embedding)
