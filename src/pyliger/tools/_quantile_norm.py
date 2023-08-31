import numpy as np
from scipy import interpolate
from scipy.stats.mstats import mquantiles

from pyliger.clustering._utilities import refine_clusts


def quantile_norm(
    liger_object,
    quantiles=50,
    ref_dataset=None,
    min_cells=20,
    dims_use=None,
    do_center=False,
    max_sample=1000,
    num_trees=None,
    refine_knn=True,
    knn_k=20,
    use_ann=False,
    rand_seed=1,
):
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
        ns = [adata.shape[0] for adata in liger_object.adata_list]
        ref_dataset_idx = np.argmax(ns)
    else:
        for i in range(num_samples):
            if liger_object.adata_list[i].uns["sample_name"] == ref_dataset:
                ref_dataset_idx = i
                break

    # set indices of factors
    if dims_use is None:
        use_these_factors = list(range(liger_object.adata_list[0].obsm["H"].shape[1]))
    else:
        use_these_factors = dims_use

    # applied use_these_factors to Hs
    Hs = [adata.obsm["H"][:, use_these_factors] for adata in liger_object.adata_list]
    num_clusters = Hs[ref_dataset_idx].shape[1]

    ### Max factor assignment
    clusters = []
    for i in range(num_samples):
        # scale the H matrix by columns, equal to scale_columns_fast in Rcpp
        if do_center:
            scale_H = (Hs[i] - np.mean(Hs[i], axis=0)) / np.std(Hs[i], axis=0, ddof=1)
        else:
            scale_H = Hs[i] / (
                np.sqrt(np.sum(np.square(Hs[i]), axis=0) / (Hs[i].shape[0] - 1))
            )

        # get the index of maximum value for each cell, equal to max_factor in Rcpp
        clusts = np.argmax(scale_H, axis=1)

        # increase robustness of cluster assignments using knn graph
        if refine_knn:
            clusts = refine_clusts(Hs[i], clusts, knn_k, use_ann, num_trees)

        clusters.append(clusts)

        # assign clusters to corresponding adata
        liger_object.adata_list[i].obs["cluster"] = clusts

    # all H_matrix used for quantile alignment
    Hs = [np.copy(adata.obsm["H"]) for adata in liger_object.adata_list]

    ### Perform quantile alignment
    for k in range(num_samples):
        for j in range(num_clusters):
            cells2 = clusters[k] == j
            cells1 = clusters[ref_dataset_idx] == j

            for i in use_these_factors:
                num_cells2 = np.sum(cells2)
                num_cells1 = np.sum(cells1)

                # skip clusters having too less cells
                if num_cells1 < min_cells or num_cells2 < min_cells:
                    continue

                if num_cells2 == 1:
                    Hs[k][cells2, i] = np.mean(Hs[ref_dataset_idx][cells1, i])
                    continue

                # maximum number of cells used for quantile normalization
                q2 = mquantiles(
                    np.random.permutation(Hs[k][cells2, i])[
                        0 : min(num_cells2, max_sample)
                    ],
                    np.linspace(0, 1, num=quantiles + 1),
                    alphap=1,
                    betap=1,
                )
                max_H = np.max(Hs[k][cells2, i])
                min_H = np.min(Hs[k][cells2, i])
                if q2[-1] < max_H:
                    q2[-1] = max_H
                if q2[0] > min_H:
                    q2[0] = min_H
                q1 = mquantiles(
                    np.random.permutation(Hs[ref_dataset_idx][cells1, i])[
                        0 : min(num_cells1, max_sample)
                    ],
                    np.linspace(0, 1, num=quantiles + 1),
                    alphap=1,
                    betap=1,
                )

                if (
                    np.sum(q1) == 0
                    or np.sum(q2) == 0
                    or len(np.unique(q1)) < 2
                    or len(np.unique(q2)) < 2
                ):
                    new_vals = np.repeat(0, num_cells2)
                else:
                    # handle ties (zeros) in order to get consistent results with LIGER
                    q1 = _mean_ties(q2, q1)
                    warp_func = interpolate.interp1d(q2, q1)
                    new_vals = warp_func(Hs[k][cells2, i])
                Hs[k][cells2, i] = new_vals

        # assign H_norm to corresponding adata
        liger_object.adata_list[k].obsm["H_norm"] = Hs[k]

    return None


def _mean_ties(x, y):
    """helper function to calculate the mean value of y where ties(zeros) occur in x"""
    idx_zeros = x == 0
    if np.sum(idx_zeros) > 0:
        y[idx_zeros] = np.mean(y[idx_zeros])

    return y
