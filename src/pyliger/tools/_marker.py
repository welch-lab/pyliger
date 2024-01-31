import numpy as np
import pandas as pd
from scipy.sparse import vstack
from sklearn.preprocessing import scale

from pyliger.tools import _wilcoxon
from pyliger.tools._metrics import calc_dataset_specificity

# Find shared and dataset-specific markers
def get_factor_markers(
    liger_object,
    dataset1=0,
    dataset2=1,
    factor_share_thresh=10,
    dataset_specificity=None,
    log_fc_thresh=1,
    umi_thresh=30,
    frac_thresh=0,
    pval_thresh=0.05,
    num_genes=30,
    print_genes=False,
):
    """Find shared and dataset-specific markers

    Applies various filters to genes on the shared (W) and dataset-specific (V)  of the
    factorization, before selecting those which load most significantly on each factor (in a shared
    or dataset-specific way).

    Parameters
    ----------
    liger_object : liger
        Should call optimizeALS before calling.
    dataset1 : str
        Name of first dataset (default first dataset by order)
    dataset2 : str
        Name of second dataset (default second dataset by order)
    factor_share_thresh : list
        Use only factors with a dataset specificity less than or equal to threshold (default 10).
    dataset_specificity : TYPE
        Pre-calculated dataset specificity if available. Will calculate if not available.
    log_fc_thresh : float
        Lower log-fold change threshold for differential expression in markers (default 1).
    umi_thresh : int
        Lower UMI threshold for markers (default 30).
    frac_thresh : float
        Lower threshold for fraction of cells expressing marker (default 0).
    pval_thresh : float
        Upper p-value threshold for Wilcoxon rank test for gene expression (default 0.05).
    num_genes : int
        Max number of genes to report for each dataset (default 30).
    print_genes : boolean
        Print ordered markers passing logfc, umi and frac thresholds (default FALSE).
    """
    ### Parameter setting
    if isinstance(dataset1, str):
        dataset1 = liger_object.find_dataset_idx(dataset1)
    if isinstance(dataset2, str):
        dataset2 = liger_object.find_dataset_idx(dataset2)

    if num_genes is None:
        num_genes = len(liger_object.var_genes)

    if dataset_specificity is None:
        dataset_specificity = calc_dataset_specificity(
            liger_object, dataset1=dataset1, dataset2=dataset2, do_plot=False
        )

    factors_use = np.abs(dataset_specificity[2]) <= factor_share_thresh

    if len(factors_use) < 2:
        print(
            "Warning: only {} factors passed the dataset specificity threshold.".format(
                len(factors_use)
            )
        )

    ### Extract values
    # scale H matrices
    Hs_scaled = [scale(matrix) for matrix in liger_object.H]

    # select index of max factor for each cell within the range of factors use
    labels = [np.argmax(H_scaled[:, factors_use], axis=1) for H_scaled in Hs_scaled]

    V1_matrices = []
    V2_matrices = []
    W_matrices = []
    for j in range(len(factors_use)):
        factor = factors_use[j]
        W = liger_object.W
        V1 = liger_object.V[dataset1]
        V2 = liger_object.V[dataset2]

        # if not max factor for any cell in either dataset
        if (
            np.sum(labels[dataset1] == factor) <= 1
            or np.sum(labels[dataset2] == factor) <= 1
        ):
            print(
                "Warning: factor {} did not appear as max in any cell in either dataset.".format(
                    factor
                )
            )
            continue

        # Filter genes by gene_count and cell_frac thresholds
        expr_mat = []
        for dset in [dataset1, dataset2]:
            adata = liger_object.adata_list[dset]
            idx = labels[dset] == factor
            expr_mat.append(adata.layers["norm_data"][idx, :])
        expr_mat = vstack(expr_mat)
        cell_label = np.concatenate(
            np.repeat(dataset1, np.sum(labels[dataset1] == factor)),
            np.repeat(dataset1, np.sum(labels[dataset2] == factor)),
        )
        wilcoxon_result = _wilcoxon(np.log(expr_mat.toarray() + 1e-10), cell_label)

        # filter based on log-fold change
        # log2fc = wilcoxon_result[wilcoxon_result['group'] == dataset1]['logFC'].to_numpy()
        filtered_genes_V1 = wilcoxon_result[
            (wilcoxon_result["logFC"] > log_fc_thresh)
            & (wilcoxon_result["pval"] < pval_thresh)
        ]["feature"].to_numpy()
        filtered_genes_V2 = wilcoxon_result[
            (-wilcoxon_result["logFC"] > log_fc_thresh)
            & (wilcoxon_result["pval"] < pval_thresh)
        ]["feature"].to_numpy()

        W = np.minimum((W + V1), (W + V2))
        V1 = V1[filtered_genes_V1, :]
        V2 = V2[filtered_genes_V2, :]

        if np.sum(filtered_genes_V1) == 0:
            top_genes_V1 = None
        else:
            top_genes_V1 = liger_object.adata_list[dataset1].var.index[
                np.argsort(V1[:, factor])[0:num_genes]
            ]

        if len(filtered_genes_V2) == 0:
            top_genes_V2 = None
        else:
            top_genes_V2 = liger_object.adata_list[dataset2].var.index[
                np.argsort(V2[:, factor])[0:num_genes]
            ]

        top_genes_W = liger_object.adata_list[dataset1].var.index[
            np.argsort(W[:, factor])[0:num_genes]
        ]

        if print_genes:
            print("Factor: ", factor)
            print("Dataset 1")
            print(top_genes_V1)
            print("Shared")
            print(top_genes_W)
            print("Dataset 2")
            print(top_genes_V2)

        pvals = []  # order is V1, V2, W
        top_genes = [top_genes_V1, top_genes_V2, top_genes_W]

        for idx, gene_list in enumerate(top_genes):
            pvals.append(
                wilcoxon_result[
                    (wilcoxon_result["feature"].isin(gene_list))
                    & (wilcoxon_result["group"] == dataset1)
                ]["pval"].to_numpy()
            )

        V1_matrices.append(
            pd.DataFrame(
                {
                    "factor_num": np.repeat(factor, len(top_genes_V1)),
                    "gene": top_genes_V1,
                    "log2fc": wilcoxon_result[
                        (wilcoxon_result["group"] == dataset1)
                        & (wilcoxon_result["feature"].isin(top_genes_V1))
                    ]["logFC"].to_numpy(),
                    "p_value": pvals[0],
                }
            )
        )
        V2_matrices.append(
            pd.DataFrame(
                {
                    "factor_num": np.repeat(factor, len(top_genes_V2)),
                    "gene": top_genes_V1,
                    "log2fc": wilcoxon_result[
                        (wilcoxon_result["group"] == dataset1)
                        & (wilcoxon_result["feature"].isin(top_genes_V2))
                    ]["logFC"].to_numpy(),
                    "p_value": pvals[1],
                }
            )
        )
        W_matrices.append(
            pd.DataFrame(
                {
                    "factor_num": np.repeat(factor, len(top_genes_W)),
                    "gene": top_genes_V1,
                    "log2fc": wilcoxon_result[
                        (wilcoxon_result["group"] == dataset1)
                        & (wilcoxon_result["feature"].isin(top_genes_W))
                    ]["logFC"].to_numpy(),
                    "p_value": pvals[2],
                }
            )
        )

    V1_genes = pd.concat(V1_matrices)
    V2_genes = pd.concat(V2_matrices)
    W_genes = pd.concat(W_matrices)

    output_list = [V1_genes, W_genes, V2_genes]
    # cutoff only applies to dataset-specific dfs
    for idx, df in enumerate(output_list):
        if idx != 1:
            output_list[idx] = df[df["p_value"] < pval_thresh]

    return output_list


"""
def get_factor_markers(liger_object,
                       dataset1 = 0,
                       dataset2 = 1,
                       factor_share_thresh = 10,
                       dataset_specificity = None,
                       log_fc_thresh = 1,
                       umi_thresh = 30,
                       frac_thresh = 0,
                       pval_thresh = 0.05,
                       num_genes = 30,
                       print_genes = False):
    ### Parameter setting
    if isinstance(dataset1, str):
        dataset1 = liger_object.find_dataset_idx(dataset1)
    if isinstance(dataset2, str):
        dataset2 = liger_object.find_dataset_idx(dataset2)

    if num_genes is None:
        num_genes = len(liger_object.var_genes)

    if dataset_specificity is None:
        dataset_specificity = calc_dataset_specificity(liger_object, dataset1=dataset1, dataset2=dataset2, do_plot=False)

    factors_use = np.abs(dataset_specificity[2]) <= factor_share_thresh

    if len(factors_use) < 2:
        print('Warning: only {} factors passed the dataset specificity threshold.'.format(len(factors_use)))

    ### Extract values
    # scale H matrices
    Hs_scaled = [scale(matrix) for matrix in liger_object.H]

    # select index of max factor for each cell within the range of factors use
    labels = [np.argmax(H_scaled[:, factors_use], axis=1) for H_scaled in Hs_scaled]

    V1_matrices = []
    V2_matrices = []
    W_matrices = []
    for j in range(len(factors_use)):
        factor = factors_use[j]
        W = liger_object.W
        V1 = liger_object.V[dataset1]
        V2 = liger_object.V[dataset2]

        # if not max factor for any cell in either dataset
        if np.sum(labels[dataset1] == factor) <= 1 or np.sum(labels[dataset2] == factor) <= 1:
            print('Warning: factor {} did not appear as max in any cell in either dataset.'.format(factor))
            continue

        # Filter genes by gene_count and cell_frac thresholds
        gene_info = []
        for dset in [dataset1, dataset2]:
            adata = liger_object.adata_list[dset]
            idx = labels[dset] == factor
            gene_counts = np.ravel(np.sum(adata.X[idx, adata.uns['var_gene_idx']], axis=0))
            cell_fracs = np.ravel(np.sum(adata.X[idx, adata.uns['var_gene_idx']] > 0, axis=0)) / np.sum(idx)
            norm_counts = adata.layers['norm_data'][idx, :]
            mean_value = np.mean(norm_counts[:, adata.uns['var_gene_idx']], axis=0)
            gene_info.append(pd.DataFrame({'gene_counts': gene_counts, 'cell_fracs': cell_fracs, 'norm_counts': norm_counts, 'mean_value': mean_value}))

        # filter based on umi_thresh
        idx_gene_counts = np.logical_or(gene_info[0]['gene_counts'] > umi_thresh, gene_info[1]['gene_counts'] > umi_thresh)
        # filter based on fraction of cells
        idx_cell_fracs = np.logical_or(gene_info[0]['cell_fracs'] > frac_thresh, gene_info[1]['cell_fracs'] > frac_thresh)
        idx = np.logical_and(idx_gene_counts, idx_cell_fracs)  # bool index
        initial_filtered = liger_object.var_genes[idx]  # gene string name

        # filter based on log-fold change
        log2fc = np.log2(gene_info[0]['mean_value'] / gene_info[1]['mean_value'])
        idx_log2fc_v1 = log2fc > log_fc_thresh
        idx_log2fc_v2 = log2fc * (-1) > log_fc_thresh

        filtered_genes_V1 = np.logical_and(idx_log2fc_v1, idx) # bool index
        filtered_genes_V2 = np.logical_and(idx_log2fc_v2, idx) # bool index

        W = np.minimum((W+V1), (W+V2))
        V1 = V1[filtered_genes_V1, :]
        V2 = V2[filtered_genes_V2, :]

        if np.sum(filtered_genes_V1) == 0:
            top_genes_V1 = None
        else:
            top_genes_V1 = initial_filtered[np.argsort(V1[:, factor])[0:num_genes]] # gene name string

        if len(filtered_genes_V2) == 0:
            top_genes_V2 = None
        else:
            top_genes_V2 = initial_filtered[np.argsort(V2[:, factor])[0:num_genes]] # gene name string

        top_genes_W = initial_filtered[np.argsort(W[:, factor])[0:num_genes]]

        if print_genes:
            print('Factor: ', factor)
            print('Dataset 1')
            print(top_genes_V1)
            print('Shared')
            print(top_genes_W)
            print('Dataset 2')
            print(top_genes_V2)

        pvals = [] # order is V1, V2, W
        top_genes = [top_genes_V1, top_genes_V2, top_genes_W]

        for idx, gene_list in enumerate(top_genes):
            if gene_list is None:
                pvals.append(None)
                continue

            for gene in gene_list:
                idx_1 = liger_object.adata_list[dataset1].var.index.get_loc(gene)
                idx_2 = liger_object.adata_list[dataset2].var.index.get_loc(gene)
                pvals.append(mannwhitneyu(gene_info[dataset1]['norm_counts'][idx_1, :], gene_info[dataset2]['norm_counts'][idx_2, :]).pvalue)

        V1_matrices.append(pd.DataFrame([np.repeat(factor, len(top_genes_V1)),
                                         top_genes_V1,
                                         gene_info[dataset1]['gene_counts'][filtered_genes_V1],
                                         gene_info[dataset2]['gene_counts'][filtered_genes_V1],
                                         gene_info[dataset1]['cell_fracs'][filtered_genes_V1],
                                         gene_info[dataset2]['cell_fracs'][filtered_genes_V1],
                                         log2fc[filtered_genes_V1],
                                         pvals[0]]))
        V2_matrices.append(pd.DataFrame([np.repeat(factor, len(top_genes_V2)),
                                         top_genes_V1,
                                         gene_info[dataset1]['gene_counts'][filtered_genes_V2],
                                         gene_info[dataset2]['gene_counts'][filtered_genes_V2],
                                         gene_info[dataset1]['cell_fracs'][filtered_genes_V2],
                                         gene_info[dataset2]['cell_fracs'][filtered_genes_V2],
                                         log2fc[filtered_genes_V1],
                                         pvals[1]]))
        W_matrices.append(pd.DataFrame([np.repeat(factor, len(top_genes_W)),
                                        top_genes_W,
                                        gene_info[dataset1]['gene_counts'][filtered_genes_V1],
                                        gene_info[dataset2]['gene_counts'][filtered_genes_V1],
                                        gene_info[dataset1]['cell_fracs'][filtered_genes_V1],
                                        gene_info[dataset2]['cell_fracs'][filtered_genes_V1],
                                        log2fc[filtered_genes_V1],
                                        pvals[2]]))
    V1_genes =
    V2_genes =
    W_genes =

    df_cols = ['factor_num', 'gene', 'counts1', 'counts2', 'fracs1', 'fracs2', 'log2fc', 'p_value']
    output_list = [V1_genes, W_genes, V2_genes]

    # cutoff only applies to dataset-specific dfs


    return output_list

"""
