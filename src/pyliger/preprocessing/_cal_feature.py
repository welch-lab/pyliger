from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix


def make_feature_matrix(
    file_dir, gene_file, promoter_file, filter_threshold, return_separate
):
    """

    :param file_dir:
    :param gene_file:
    :param promoter_file:
    :param filter_threshold:
    :return:
    """
    ### 0. build file path
    p_gene = Path(file_dir, gene_file)
    p_promoter = Path(file_dir, promoter_file)

    ### 1. extract barcodes (equivalent to step 3 from R tutorial)
    genes_bc = pd.read_csv(p_gene, sep="\t", comment="t", header=None, usecols=[3, 6])
    promoters_bc = pd.read_csv(
        p_promoter, sep="\t", comment="t", header=None, usecols=[3, 6]
    )
    genes_bc = genes_bc.sort_values(by=[3])
    promoters_bc = promoters_bc.sort_values(by=[3])
    num_genes = genes_bc.shape[0]

    ### 2. make matrices for gene counts and promoter counts (equivalent to step 4 from R tutorial)
    genes_bc_total, genes_bc_counts, genes_indptr, _ = _extract_barcodes(genes_bc)
    promoters_bc_total, promoters_bc_counts, _, promoters_num_bc = _extract_barcodes(
        promoters_bc
    )

    genes_bc_unique, genes_bc_inverse = np.unique(genes_bc_total, return_inverse=True)

    # select barcodes
    genes_mtx = csr_matrix(
        (genes_bc_counts, genes_bc_inverse, genes_indptr),
        shape=(num_genes, len(genes_bc_unique)),
    )
    filter_idx = np.sum(genes_mtx, axis=0) > filter_threshold
    genes_mtx = genes_mtx[:, filter_idx.A1]

    row_idx = np.repeat(np.arange(num_genes), promoters_num_bc)
    mask = np.isin(promoters_bc_total, genes_bc_unique[filter_idx.A1])

    promoters_bc_total = promoters_bc_total[mask]
    promoters_bc_counts = promoters_bc_counts[mask]
    row_idx = row_idx[mask]
    promoters_bc_unique, promoters_bc_inverse = np.unique(
        promoters_bc_total, return_inverse=True
    )
    promoters_mtx = csr_matrix(
        (promoters_bc_counts, (row_idx, promoters_bc_inverse)),
        shape=(num_genes, len(promoters_bc_unique)),
    )

    if return_separate:
        return None
    else:
        sum_mtx = genes_mtx + promoters_mtx
        barcodes = pd.DataFrame(genes_bc_unique[filter_idx.A1])
        barcodes.columns = ["barcodes"]
        barcodes = barcodes.set_index("barcodes")

        gene_name = pd.DataFrame(genes_bc[3])
        gene_name.columns = ["gene_name"]
        gene_name = gene_name.set_index("gene_name")

        atac_anndata = AnnData(
            csr_matrix(sum_mtx.transpose()),
            obs=barcodes,
            var=gene_name,
            uns={"sample_name": "test_1", "data_type": "ATAC"},
        )
        return atac_anndata


def _extract_barcodes(bedmat):
    """"""
    barcodes_total = []
    barcodes_value = []
    indptr = [0]
    num_barcodes = []
    for i in range(bedmat.shape[0]):
        row = bedmat.iloc[i]

        if type(row[6]) == str:
            barcodes = row[6].split(";")
        else:
            barcodes = []

        barcodes_count = Counter(barcodes)
        barcodes_total.extend(list(barcodes_count))
        barcodes_value.extend(list(barcodes_count.values()))
        indptr.append(indptr[i] + len(barcodes_count))
        num_barcodes.append(len(barcodes_count))

    return np.asarray(barcodes_total), np.asarray(barcodes_value), indptr, num_barcodes
