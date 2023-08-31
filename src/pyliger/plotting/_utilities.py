import numpy as np


def get_gene_values(
    liger_object, gene, methylation_indices=None, log2scale=False, scale_factor=10000
):
    """"""

    if methylation_indices is None:
        methylation_indices = []

    gene_vals_total = []
    for idx, adata in enumerate(liger_object.adata_list):
        if adata.isbacked:
            gene_vals = _get_gene_values_disk(adata, gene)  # To be finished
        else:
            gene_vals = _get_gene_values_memory(adata, gene)

        if log2scale and idx not in methylation_indices:
            gene_vals = np.log2(gene_vals * scale_factor + 1)

        gene_vals_total.append(gene_vals)

    return np.concatenate(gene_vals_total)


def _get_gene_values_memory(adata, gene):
    # if gene is highly variable gene
    if gene in adata.var.index:
        gene_vals = np.ravel(adata[:, gene].layers["norm_data"].toarray())

    # recover gene value from raw backup
    elif gene in adata.raw.var.index:
        from sklearn.preprocessing import normalize as sp_normalize

        normalized_data = sp_normalize(adata.raw.X)
        idx = adata.raw.var.index.get_loc(gene)
        gene_vals = np.ravel(normalized_data[:, idx].toarray())

    # use 0 for not existing gene
    else:
        gene_vals = np.zeros(adata.shape[0], dtype=np.int32)
    return gene_vals


def _get_gene_values_disk():
    return None
