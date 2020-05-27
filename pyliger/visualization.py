#######################################################################################
#### Dimensionality Reduction
    
# Perform t-SNE dimensionality reduction
def runTSNE(liger_object,
            dims_use, 
            use_raw = False, 
            use_pca = False, 
            perplexity = 30,
            theta = 0.5, 
            method = "Rtsne", 
            fitsne_path = None, 
            rand_seed = 42):
    dims_use = range(1, len(liger_object.H_norm))
    pass

# Perform UMAP dimensionality reduction
def runUMAP(liger_object,
            dims_use, 
            use_raw = False, 
            k = 2, 
            distance = "euclidean",
            n_neighbors = 10,
            min_dist = 0.1, 
            rand_seed = 42):
    dims_use = range(1, len(liger_object.H_norm))
    pass