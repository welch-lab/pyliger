import os
import time
import numpy as np
import pandas as pd
import louvain, leidenalg
from scipy import interpolate
from scipy.stats import ranksums
from scipy.sparse import csr_matrix, vstack


from ._utilities import refine_clusts, compute_snn, run_knn, build_igraph
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
        Whether to use approximate nearest neighbor (the default is False)
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
        use_these_factors = list(range(liger_object.adata_list[0].obsm['H'].shape[1]))
    else:
        use_these_factors = dims_use
    
    # applied use_these_factors to Hs
    Hs = [adata.obsm['H'][:, use_these_factors] for adata in liger_object.adata_list]
    num_clusters = Hs[ref_dataset_idx].shape[1]

    ### Max factor assignment
    clusters = []
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
        
        # assign clusters to corresponding adata
        liger_object.adata_list[i].obs['cluster'] = clusts
    
    # all H_matrix used for quantile alignment
    Hs = [adata.obsm['H'] for adata in liger_object.adata_list]
    
    ### Perform quantile alignment
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
        
        # assign H_norm to corresponding adata
        liger_object.adata_list[k].obsm['H_norm'] = Hs[k]         
    
    return liger_object

#from memory_profiler import profile

#@profile
def louvain_cluster(liger_object, 
                    resolution = 1.0, 
                    k = 20, 
                    prune = 1 / 15, 
                    random_seed = 1,
                    n_starts = 10):
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
    random_seed : int, optional
        Seed of the random number generator (the default is 1).
    n_starts : int, optional
        The number of random starts to be used (the default is 10). 

    Returns
    -------
    liger_object : liger object
        object with refined 'cluster'.
        
    Examples
    --------
    >>> ligerex = louvain_cluster(ligerex, resulotion = 0.3) # liger object, factorization complete
    """
    ### 1. Compute snn
    H_norm = np.vstack([adata.obsm['H_norm'] for adata in liger_object.adata_list])
    knn = run_knn(H_norm, k)
    snn = compute_snn(knn, prune=prune)

    ### 2. Get igraph from snn
    g = build_igraph(snn)
    
    ### 3. Run louvain
    np.random.seed(random_seed)
    max_quality = -1
    for i in range(n_starts): # random starts to improve stability
        seed = np.random.randint(0,1000)
        kwargs = {'weights': g.es['weight'], 'resolution_parameter': resolution, 'seed': seed} # parameters setting
        part = louvain.find_partition(g, louvain.RBConfigurationVertexPartition, **kwargs)
        
        if part.quality() > max_quality:
            cluster = part.membership
            max_quality = part.quality()
    
    ### 4. Assign cluster results
    cluster = np.array(cluster).flatten()
    idx = 0
    for i in range(len(liger_object.adata_list)):        
        liger_object.adata_list[i].obs['cluster'] = cluster[idx:(idx+liger_object.adata_list[i].shape[0])]
        idx = liger_object.adata_list[i].shape[0]
    
    return liger_object

def leiden_cluster(liger_object, 
                   resolution = 1.0, 
                   k = 20, 
                   prune = 1 / 15, 
                   random_seed = 1,
                   n_iterations = -1,
                   n_starts = 10):
    """Leiden algorithm for community detection
    
    After quantile normalization, users can additionally run the Leiden algorithm 
    for community detection.
    
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
    random_seed : int, optional
        Seed of the random number generator (the default is 1).
    n_iterations : TYPE, optional
        The number of iterations to run the Leiden algorithm.
        A negative value indicates the Leiden algorithm runs until one iterration
        that has no improvement. (the default is -1)
    n_starts : int, optional
        The number of random starts to be used (the default is 10).

    Returns
    -------
    liger_object : liger object
        object with refined 'cluster'.
    
    Examples
    --------
    >>> ligerex = leiden_cluster(ligerex, resulotion = 0.25) # liger object, factorization complete
    """
    ### 1. Compute snn
    H_norm = np.vstack([adata.obsm['H_norm'] for adata in liger_object.adata_list])
    knn = run_knn(H_norm, k)
    snn = compute_snn(knn, prune=prune)
    
    ### 2. Get igraph from snn
    g = build_igraph(snn)
    
    ### 3. Run louvain
    np.random.seed(random_seed)
    max_quality = -1
    for i in range(n_starts): # random starts to improve stability
        seed = np.random.randint(0,1000)
        kwargs = {'weights': g.es['weight'], 'resolution_parameter': resolution, 'seed': seed} # parameters setting
        part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, n_iterations=n_iterations, **kwargs)
        
        if part.quality() > max_quality:
            cluster = part.membership
            max_quality = part.quality()
    
    ### 4. Assign cluster results
    cluster = np.array(cluster).flatten()
    idx = 0
    for i in range(len(liger_object.adata_list)):        
        liger_object.adata_list[i].obs['cluster'] = cluster[idx:(idx+liger_object.adata_list[i].shape[0])]
        idx = liger_object.adata_list[i].shape[0]
    
    return liger_object

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
        Whether scale but not center the imputed data with default parameters (the default is True).

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


def run_wilcoxon(liger_object, 
                 compare_method, 
                 data_use = "all"):
    """Perform Wilcoxon rank-sum test
    
    Perform Wilcoxon rank-sum tests on specified dataset using given method.

    Parameters
    ----------
    liger_object : liger object
        object.
    compare_method : str, optional, 'clusters' or 'datasets'
        This indicates the metric of the test.
    data_use : str or list
        This selects which dataset(s) to use (the default is "all").

    Returns
    -------
    results : pd data frame
   TODO     A -columns data.frame with test results.
    
    
    # check parameter inputs
    if not compare_method == 'datasets' and not compare_method == 'clusters':
        raise ValueError('Parameter *compare.method* should be either *clusters* or *datasets*.')
    if compare_method == 'datasets':
        if len(liger_object.adata_list) < 2:
            raise ValueError('Should have at least TWO inputs to compare between datasets.')
        if isinstance(data_use, list) and len(data_use) < 2:
            raise ValueError('Should have at least TWO inputs to compare between datasets.')
    
    ### Create feature x sample matrix
    if data_use == 'all':
        num_samples = len(liger_object.adata_list)
        sample_names = [adata.uns['sample_name'] for adata in liger_object.adata_list]
        print('Performing Wilcoxon test on ALL datasets: {}'.format(', '.join(sample_names)))
    else:
        num_samples = len(data_use)
        sample_names = data_use
        sample_names_all = [adata.uns['sample_name'] for adata in liger_object.adata_list]
        sample_idx = [sample_names_all.index(name) for name in sample_names]
        print('Performing Wilcoxon test on GIVEN datasets: {}'.format(', '.join(sample_names)))
    
    # get all shared genes of every datasets
    genes = []
    for i in range(num_samples):      
        genes = np.intersect1d(genes, liger_object.adata_list[sample_idx[i]].obs['gene_name'])
     
    # TODO: change to barcodes as rows and shared genes as columns
    # get feature matrix, shared genes as rows and all barcodes as columns
    feature_matrix = []
    for i in range(num_samples):
        idx = liger_object.adata_list[sample_idx[i]].obs['gene_name'].isin(genes).to_numpy()
        feature_matrix.append(liger_object.adata_list[i].layers['norm_data'][idx,])
    feature_matrix = vstack(feature_matrix)
    
    # get labels of clusters and datasets
    cell_source = 
        
    else: # for one dataset only
        print("Performing Wilcoxon test on GIVEN dataset: ")
        feature_matrix = 
        clusters = 
            
    ### Perform wilcoxon test
    if compare_method == 'clusters':
        num_rows = feature_matrix.shape[0]
        if num_rows > 100000:
            print('Calculating Large-scale Input...')
            results = 
            
    if compare_method == 'datasets':
        results = 
        
    return results
"""
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


"""
def _save_adata(liger_object, results, dataset_use, annotation):
    helper function to save results into individual adata
    idx = 0
    for i in range(len(liger_object.adata_list)):        
        liger_object.adata_list[i].obs['cluster'] = clusters[idx:(idx+liger_object.adata_list[i].shape[0])]
        idx = liger_object.adata_list[i].shape[0]
        
    # check dimension
    if len(results.shape) == 1:
        idx = 0
        for i in range(len(liger_object.adata_list)):        
            liger_object.adata_list[i].obs['cluster'] = clusters[idx:(idx+liger_object.adata_list[i].shape[0])]
            idx = liger_object.adata_list[i].shape[0]
    elif len(results.shape) == 2:
        
    
    return liger_object
"""