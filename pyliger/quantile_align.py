import time
import numpy as np
import pandas as pd
from scipy import interpolate
from numba import jit

from ._utilities import refine_clusts
#######################################################################################
#### Quantile Alignment/Normalization
    
def quantile_norm(liger_object, 
                  quantiles = 50, 
                  ref_dataset = None, 
                  min_cells = 20, 
                  dims_use = None, 
                  do_center = False, 
                  max_sample = 1000, 
                  num_trees = None, 
                  refine_knn = True,
                  knn_k = 20,
                  use_ann = False,
                  rand_seed = 1):
    """Quantile align (normalize) factor loadings
    
    This process builds a shared factor neighborhood graph to jointly cluster cells, then quantile
    normalizes corresponding clusters.

    The first step, building the shared factor neighborhood graph, and
    produces a graph representation where edge weights between cells (across all datasets)
    correspond to their similarity in the shared factor neighborhood space. An important parameter
    here is knn_k, the number of neighbors used to build the shared factor space. 
    
    Next we perform quantile alignment for each dataset, factor, and cluster (by
    stretching/compressing datasets' quantiles to better match those of the reference dataset). These
    aligned factor loadings are combined into a single matrix and returned as H_norm.
    
    Parameters
    ----------
    liger_object : liger
        Should run optimizeALS before calling.
    quantiles : int, optional
        Number of quantiles to use for quantile normalization (the default is 50).
    ref_dataset : str, optional
        Name of dataset to use as a "reference" for normalization. By default,
        the dataset with the largest number of cells is used (the default is None).
    min_cells : int, optional
        Minimum number of cells to consider a cluster shared across datasets
        (the default is 20).
    dims_use : list, optional
        Indices of factors to use for shared nearest factor determination
        (the default is list(range(liger_object.adata_list[0].varm['H'].shape[1]))).
    do_center : bool, optional
        Centers the data when scaling factors (useful for less sparse modalities like
        methylation data) (the default is False).
    max_sample : int, optional
        Maximum number of cells used for quantile normalization of each cluster 
        and factor (the default is 1000).
    num_trees : int, optional
        The number of trees used by the approximate nearest neighbor search. 
        Larger number of trees give more precision but consume more memory.
        According to largeVis, they use 10 trees for datasets has less than 100,000
        observations, 20 trees for datasets has less than 1,000,000 observations,
        50 trees for datasets up to 5,000,000 and 100 trees otherwise (the default is None).
    refine_knn : bool, optional
        whether to increase robustness of cluster assignments using KNN graph
        (the default is True).
    knn_k : int, optional
        Number of nearest neighbors for within-dataset knn graph (the default is 20). 
    use_ann : bool, optional
        whether to use approximate nearest neighbor (the default is False)
    rand_seed : int, optional
        Random seed to allow reproducible results (the default is 1).

    Returns
    -------
    liger_object : liger
        liger_object with 'H_norm' and 'clusters' attributes.

    Examples
    --------
    ligerex = quantile_norm(ligerex) # do basic quantile alignment
    ligerex = quantile_norm(ligerex, resolution = 1.2) # higher resolution for more clusters (note that SNF is conserved)
    ligerex = quantile_norm(ligerex, knn_k = 15, resolution = 1.2) # change knn_k for more fine-grained local clustering
    """
    
    np.random.seed(rand_seed)
    
    num_samples = len(liger_object.adata_list)
    
    # set reference dataset
    if ref_dataset is None:
        ns = [adata.shape[1] for adata in liger_object.adata_list]
        ref_dataset_idx = np.argmax(ns)
    else:
        for i in range(num_samples):
            if liger_object.adata_list[i].uns['sample_name'] == ref_dataset:
                ref_dataset_idx = i
                break

    # set indices of factors
    if dims_use is None:
        use_these_factors = list(range(liger_object.adata_list[0].varm['H'].shape[1]))
    else:
        use_these_factors = dims_use
    
    # applied use_these_factors to Hs
    #Hs = [adata.varm['H'][:,use_these_factors] for adata in liger_object.adata_list]
    Hs = [np.loadtxt('/Users/lulu/Desktop/result_H1.txt')[:,use_these_factors], np.loadtxt('/Users/lulu/Desktop/result_H2.txt')[:,use_these_factors]]
    num_clusters = Hs[ref_dataset_idx].shape[1]

    # Max factor assignment
    clusters = []
    col_names = []
    for i in range(num_samples):
        # scale the H matrix by columns, equal to scale_columns_fast in Rcpp
        if do_center:
            scale_H = (Hs[i]-np.mean(Hs[i], axis=0))/np.std(Hs[i], axis=0, ddof=1)
        else:
            scale_H = Hs[i]/(np.sqrt(np.sum(np.square(Hs[i]), axis=0)/(Hs[i].shape[0]-1)))

        # get the index of maximum value for each cell, equal to max_factor in Rcpp
        clusts = np.argmax(scale_H, axis=1)
        
        # increase robustness of cluster assignments using knn graph
        clusts = refine_clusts(Hs[i], clusts, knn_k, use_ann, num_trees)
        
        clusters.append(clusts)
        col_names.append(liger_object.adata_list[i].var['barcodes'])
     
    # all H_matrix used for quantile alignment
    #Hs = [adata.varm['H'] for adata in liger_object.adata_list]
    Hs = [np.loadtxt('/Users/lulu/Desktop/result_H1.txt'), np.loadtxt('/Users/lulu/Desktop/result_H2.txt')]
    # Perform quantile alignment
    for k in range(num_samples):
        for j in range(num_clusters):
            cells2 = clusters[k] == j
            cells1 = clusters[ref_dataset_idx] == j

            for i in range(num_clusters):
                num_cells2 = np.sum(cells2)
                num_cells1 = np.sum(cells1)
                
                # skip clusters having too less cells
                if num_cells1 < min_cells or num_cells2 < min_cells:
                    continue
                
                if num_cells2 == 1:
                    Hs[k][cells2, i] = np.mean(Hs[ref_dataset_idx][cells1,i])
                    continue
                
                # maxiumn number of cells used for quantile normalization
                q2 = np.quantile(np.random.permutation(Hs[k][cells2, i])[0:min(num_cells2, max_sample)], np.linspace(0,1,num=quantiles+1))
                q1 = np.quantile(np.random.permutation(Hs[ref_dataset_idx][cells1, i])[0:min(num_cells1, max_sample)], np.linspace(0,1,num=quantiles+1))

                if np.sum(q1) == 0 or np.sum(q2) == 0 or len(np.unique(q1)) < 2 or len(np.unique(q2)) < 2:
                    new_vals = np.repeat(0, num_cells2) 
                else: 
                    warp_func = interpolate.interp1d(q2, q1)
                    new_vals = warp_func(Hs[k][cells2, i])
                Hs[k][cells2, i] = new_vals
                
        if k == 0:
            H_norm = Hs[k]
        else:
            H_norm = np.vstack((H_norm, Hs[k]))
    
    # combine clusters into one
    clusters = np.array(clusters).flatten()
    col_names = np.array(col_names).flatten()
    
    # assign clusters and H_norm attributes to liger_object
    liger_object.clusters = pd.DataFrame(clusters, index=col_names)
    liger_object.H_norm = pd.DataFrame(H_norm, index=col_names)
    
    return liger_object


def louvain_cluster(liger_object, 
                   resolution = 1.0, 
                   k = 20, 
                   prune = 1 / 15, 
                   eps = 0.1, 
                   nRandomStarts = 10,
                   nIterations = 100, 
                   random_seed = 1):
    """Louvain algorithm for community detection
    
    After quantile normalization, users can additionally run the Louvain algorithm 
    for community detection, which is widely used in single-cell analysis and excels at merging 
    small clusters into broad cell classes.

    Parameters
    ----------
    liger_object : liger object
        Should run quantile_norm before calling.
    resolution : float, optional
        Value of the resolution parameter, use a value above (below) 1.0 if you want 
        to obtain a larger (smaller) number of communities (the default is 1.0).
    k : int, optional
        The maximum number of nearest neighbours to compute (the default is 20).
    prune : float, optional
        Sets the cutoff for acceptable Jaccard index when
        computing the neighborhood overlap for the SNN construction. Any edges with
        values less than or equal to this will be set to 0 and removed from the SNN
        graph. Essentially sets the strigency of pruning (0 --- no pruning, 1 ---
        prune everything) (the default is 1/15).
    eps : float, optional
        The error bound of the nearest neighbor search (the default is 0.1).
    nRandomStarts : int, optional
        Number of random starts (the default is 10).
    nIterations : int, optional
        Maximal number of iterations per random start (the default is 100).
    random_seed : int, optional
        Seed of the random number generator (the default is 1).

    Returns
    -------
    liger_object : liger object
        object with refined 'clusters' attribute.
        
    Examples
    --------
    >>> ligerex = louvainCluster(ligerex, resulotion = 0.3) # liger object, factorization complete
    """
    
    current_time = time.strftime("%Y-%m-%d_%H_%M_%S", time.localtime())
    output_path = 'edge_' + current_time + '.txt'
    
   



def imputeKNN(liger_object, 
              reference, 
              queries, 
              knn_k = 20,
              weight = True, 
              norm = True, 
              scale = False):
    """Impute the query cell expression matrix
    
    Impute query features from a reference dataset using KNN.

    Parameters
    ----------
    liger_object : liger obect
        object.
    reference : str
        Dataset containing values to impute into query dataset(s).
    queries : TYPE
        Dataset to be augmented by imputation. If not specified, will pass in all datasets.
    knn_k : int, optional
        The maximum number of nearest neighbors to search (the default is 20).
    weight : bool, optional
        Whether to use KNN distances as weight matrix (the default is False).
    norm : bool, optional
        Whether normalize the imputed data with default parameters (the default is True).
    scale : bool, optional
        Whether scale but not center the imputed data with default parameters (default True).

    Returns
    -------
    liger_object : liger obect
        object with raw data in raw.data slot replaced by imputed data (genes by cells)
        # TODO: may not be replaced?
    
    Examples
    --------
    >>> 
    """
    pass

# Perform Wilcoxon rank-sum test
def runWilcoxon(liger_object, compare_method, data_use = "all"):
    pass

# Linking genes to putative regulatory elements
def linkGenesAndPeaks(gene_counts, peak_counts, path_to_coords, genes_list = None, dist = "spearman", 
                      alpha = 0.05):
    pass

# Export predicted gene-pair interaction
def makeInteractTrack(corr_mat, genes_list, output_path, path_to_coords):
    pass

# Analyze biological interpretations of metagene
def runGSEA(liger_object, gene_sets = [], mat_w = True, mat_v = 0, custom_gene_sets = []):
    pass
