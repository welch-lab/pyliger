from functools import reduce

import numpy as np
import pandas as pd
import statsmodels.stats.multitest as smt
from scipy.sparse import vstack
from scipy.stats import distributions
from sklearn.preprocessing import normalize as sp_normalize


def run_wilcoxon(liger_object, compare_method, data_use="all"):
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
    """
    ### Check parameter inputs
    if compare_method not in ["datasets", "clusters"]:
        raise ValueError(
            "Parameter *compare.method* should be either *clusters* or *datasets*."
        )
    if compare_method == "datasets":
        if len(liger_object.adata_list) < 2:
            raise ValueError(
                "Should have at least TWO inputs to compare between datasets."
            )
        if isinstance(data_use, list) and len(data_use) < 2:
            raise ValueError(
                "Should have at least TWO inputs to compare between datasets."
            )

    ### Create feature x sample matrix
    if data_use == "all":
        num_samples = len(liger_object.adata_list)
        sample_names = [adata.uns["sample_name"] for adata in liger_object.adata_list]
        sample_idx = list(range(num_samples))
        print(
            "Performing Wilcoxon test on ALL datasets: {}".format(
                ", ".join(sample_names)
            )
        )
    else:
        num_samples = len(data_use)
        sample_names = data_use
        sample_names_all = [
            adata.uns["sample_name"] for adata in liger_object.adata_list
        ]
        sample_idx = [sample_names_all.index(name) for name in sample_names]
        print(
            "Performing Wilcoxon test on GIVEN datasets: {}".format(
                ", ".join(sample_names)
            )
        )

    # get all shared genes of every datasets
    # genes_all = [liger_object.adata_list[idx].var.index.to_numpy() for idx in sample_idx]
    genes_all = [
        liger_object.adata_list[idx].raw.var.index.to_numpy() for idx in sample_idx
    ]
    genes_shared = reduce(np.intersect1d, genes_all)

    # get feature matrix, shared genes as columns and all barcodes as rows, as well as labels of clusters and datasets
    feature_matrix = []
    clusters = np.array([])
    cell_source = np.array([])
    for idx in sample_idx:
        # _, gene_idx, _ = np.intersect1d(liger_object.adata_list[idx].var.index, genes_shared, return_indices=True)
        _, gene_idx, _ = np.intersect1d(
            liger_object.adata_list[idx].raw.var.index,
            genes_shared,
            return_indices=True,
        )
        gene_idx = np.sort(
            gene_idx
        )  # make sure data are extracted in the original order
        # gene_name = liger_object.adata_list[idx].var.index[gene_idx]
        gene_name = liger_object.adata_list[idx].raw.var.index[gene_idx]
        raw_data = liger_object.adata_list[idx].raw.X[:, gene_idx]
        norm_data = sp_normalize(raw_data, axis=1, norm="l1")
        feature_matrix.append(norm_data)
        clusters = np.concatenate(
            (clusters, liger_object.adata_list[idx].obs["cluster"]), axis=None
        )
        cell_source = np.concatenate(
            (
                cell_source,
                np.repeat(
                    liger_object.adata_list[idx].uns["sample_name"],
                    liger_object.adata_list[idx].shape[0],
                ),
            ),
            axis=None,
        )
    feature_matrix = vstack(feature_matrix)

    ### Perform wilcoxon test
    if compare_method == "clusters":  # compare between clusters across datasets
        num_rows = feature_matrix.shape[1]
        if num_rows > 100000:
            print("Calculating Large-scale Input...")
        results = _wilcoxon(np.log(feature_matrix.toarray() + 1e-10), clusters)
        gene_label = np.tile(gene_name, len(np.unique(clusters)))
        cluster_label = np.repeat(np.unique(clusters), len(gene_name))
    elif compare_method == "datasets":  # compare between datasets within each cluster
        groups = np.unique(clusters)
        cluster_label = np.array([])
        gene_label = np.array([])
        data = []
        for group in groups:
            sub_idx = clusters == group
            sub_clusters = clusters[sub_idx]
            sub_cell_source = cell_source[sub_idx]
            sub_label = np.array(
                [
                    str(int(sub_clusters[i])) + "-" + sub_cell_source[i]
                    for i in range(len(sub_clusters))
                ]
            )
            sub_matrix = feature_matrix[sub_idx, :]
            if len(np.unique(sub_cell_source)) == 1:
                print(
                    "Note: Skip Cluster {} since it has only ONE data source.".format(
                        group
                    )
                )
            else:
                data.append(_wilcoxon(np.log(sub_matrix.toarray() + 1e-10), sub_label))
                gene_label = np.concatenate(
                    (gene_label, np.tile(gene_name, len(np.unique(sub_label)))),
                    axis=None,
                )
                cluster_label = np.concatenate(
                    (cluster_label, np.repeat(np.unique(sub_label), len(gene_name))),
                    axis=None,
                )
        results = pd.concat(data, ignore_index=True)

    results.insert(0, "group", cluster_label)
    results.insert(0, "feature", gene_label)
    return results


def _wilcoxon(X, y, groups_use=None, verbose=True, lable_use=None):
    ### Check and possibly correct input values
    if X.shape[0] != len(y):
        raise ValueError("Number of columns of X does not match length of y")

    if groups_use is not None:
        idx_use = np.isin(groups_use, y)
        X = X[idx_use, :]
        y = y[idx_use]

    idx_use = pd.isnull(y)
    if np.sum(~idx_use) < len(y):
        X = X[~idx_use, :]
        y = y[~idx_use]
        if verbose:
            print("Removing NA values from labels")

    ### Compute primary statistics
    clusters, cluster_size = np.unique(y, return_counts=True)
    num_genes = X.shape[1]
    num_clusters = len(clusters)
    statistic = []
    pval = []
    padj = []
    sum_table = np.empty((num_clusters, num_genes), dtype=np.float64)
    # ranks = rankdata(X, axis=0)
    temp_z = []
    for i in range(num_genes):
        ranks, T = _rank(X[:, i])
        temp_pval = []
        for j in range(num_clusters):
            idx = y == clusters[j]
            sum_table[j, i] = np.sum(X[idx, i])
            u, z = _mannwhitneyu(ranks[idx], ranks[~idx], T)
            statistic.append(u)
            temp_pval.append(z)
            temp_z.append(z)
        pval.extend(temp_pval)
        padj.extend(smt.multipletests(temp_pval, method="fdr_bh")[1])
    pval = 2 * distributions.norm.sf(np.abs(temp_z))

    ### Auxiliary Statistics (AvgExpr, LFC, etc)
    cs = np.sum(sum_table, axis=0, dtype=np.float64)
    avg_expr = np.divide(
        sum_table,
        np.reshape(cluster_size, (cluster_size.shape[0], 1)),
        dtype=np.float64,
    )
    lfc = np.ravel(
        [
            avg_expr[i, :] - ((cs - sum_table[i, :]) / (len(y) - cluster_size[i]))
            for i in range(num_clusters)
        ]
    )

    summary = {
        "avgExpr": np.ravel(avg_expr),
        "logFC": lfc,
        "statistic": np.asarray(statistic),
        "pval": np.asarray(pval),
        "padj": np.asarray(padj),
    }
    return pd.DataFrame(data=summary)


from numba import jit, njit


@njit
def _mannwhitneyu(X1, X2, T):
    """"""
    n1 = len(X1)
    n2 = len(X2)

    u1 = n1 * n2 + (n1 * (n1 + 1)) / 2.0 - np.sum(X1)
    u2 = n1 * n2 - u1

    sd = np.sqrt(T * n1 * n2 * (n1 + n2 + 1) / 12.0)

    meanrank = n1 * n2 / 2.0 + 0.5
    bigu = max(u1, u2)

    z = (bigu - meanrank) / sd
    # p = 2 * distributions.norm.sf(np.abs(z))
    u = min(u1, u2)
    return u, z


def _rank(arr):
    """
    Parameters
    ----------
    arr : numpy array
        array needs to be ranked

    Returns
    -------
    ranks : numpy array

    """
    sorter = np.argsort(arr)
    inv = np.empty(len(sorter), dtype=np.intp)
    inv[sorter] = np.arange(len(sorter), dtype=np.intp)
    arr = arr[sorter]

    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]
    count = np.r_[np.nonzero(obs)[0], len(obs)]

    ranks = 0.5 * (count[dense] + count[dense - 1] + 1)

    cnt = np.diff(count).astype(np.float64)
    size = np.float64(arr.size)

    T = 1.0 - (cnt**3 - cnt).sum() / (size**3 - size)

    return ranks, T


"""
    for i in range(num_clusters):
        idx = y == clusters[i]
        ranks = _rank(X[])
        sum_table[i, :] = np.sum(X[idx, :], axis=0, dtype=np.float64)
        temp_pval = []
        for j in range(num_genes):
            try:
                #u, p = mannwhitneyu(X[idx, j], X[~idx, j], alternative='two-sided')
                u, p = _mannwhitneyu(ranks[idx, j], ranks[~idx, j], alternative='two-sided')
            except ValueError:
                u, p = -1.0, 1.0
            statistic.append(u)
            temp_pval.append(p)
        pval.extend(temp_pval)
        padj.extend(smt.multipletests(temp_pval, method='fdr_bh')[1])
"""

"""
def wilcoxon(X, y, groups_use=None, verbose=True, lable_use=None):
    #helper function for wilcoxon tests on single-cell datasets (cell by gene matrix)
    ### Check and possibly correct input values
    if X.shape[0] != len(y):
        raise ValueError('Number of columns of X does not match length of y')

    if groups_use is not None:
        idx_use = np.isin(groups_use, y)
        X = X[idx_use, :]
        y = y[idx_use]

    idx_use = pd.isnull(y)
    if np.sum(~idx_use) < len(y):
        X = X[~idx_use, :]
        y = y[~idx_use]
        if verbose:
            print('Removing NA values from labels')

    ### Compute primary statistics
    clusters, cluster_size = np.unique(y, return_counts=True)
    num_genes = X.shape[1]
    num_clusters = len(clusters)
    statistic = []
    pval = []
    padj = []
    sum_table = np.zeros((num_clusters, num_genes), dtype=np.float64)

    for i in range(num_clusters):
        idx = y == clusters[i]
        sum_table[i, :] = np.sum(X[idx, :], axis=0, dtype=np.float64)
        temp_pval = []
        for j in range(num_genes):
            try:
                u, p = mannwhitneyu(X[idx, j], X[~idx, j], alternative='two-sided')
            except ValueError:
                u, p = -1.0, 1.0
            statistic.append(u)
            temp_pval.append(p)
        pval.extend(temp_pval)
        padj.extend(smt.multipletests(temp_pval, method='fdr_bh')[1])

    ### Auxiliary Statistics (AvgExpr, LFC, etc)
    cs = np.sum(sum_table, axis=0, dtype=np.float64)
    avg_expr = np.divide(sum_table, np.reshape(cluster_size, (cluster_size.shape[0], 1)), dtype=np.float64)
    lfc = np.ravel([avg_expr[i, :] - ((cs - sum_table[i, :]) / (len(y) - cluster_size[i]))
                      for i in range(num_clusters)])

    summary = {'avgExpr': np.ravel(avg_expr),
               'logFC': lfc,
               'statistic': np.asarray(statistic),
               'pval': np.asarray(pval),
               'padj': np.asarray(padj)}
    return pd.DataFrame(data=summary)
"""
