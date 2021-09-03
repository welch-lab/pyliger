import numpy as np
from anndata import AnnData
from scipy.sparse import csr_matrix

from ..pyliger import Liger
from .._utilities import _remove_missing_obs, _h5_idx_generator, _merge_sparse_data_all

from typing import List, Optional


def create_liger(adata_list: List,
                 make_sparse: bool = True,
                 take_gene_union: bool = False,
                 remove_missing: bool = True,
                 chunk_size: Optional[int] = 1000
                 ) -> Liger:
    """Create a liger object.

    This function initializes a liger object with the raw data passed in. It requires a list of
    expression (or another single-cell modality) matrices (cell by gene) for at least two datasets.
    By default, it converts all passed data into Compressed Sparse Row matrix (CSR matrix) to reduce
    object size. It initializes cell_data with nUMI and nGene calculated for every cell.

    Parameters
    ----------
    adata_list : list
        List of AnnData objects which store expression matrices (cell by gene).
        Should be named by dataset.
    make_sparse : bool, optional
        Whether to convert raw_data into sparse matrices (the default is True).
    take_gene_union : bool, optional
        Whether to fill out raw_data matrices with union of genes across all
        datasets (filling in 0 for missing data) (requires make_sparse=True)
        (the default is False).
    remove_missing : bool, optional
        Whether to remove cells not expressing any measured genes, and genes not
        expressed in any cells (if take_gene_union=True, removes only genes not
        expressed in any dataset) (the default is True).
    chunk_size : int, optional


    Returns
    -------
    liger_object : liger object
        object with raw_data slot set.
    TODO: update the docstring for returns

    Examples
    --------
    >>> adata1 = AnnData(np.arange(12).reshape((4, 3)))
    >>> adata2 = AnnData(np.arange(12).reshape((4, 3)))
    >>> ligerex = pyliger.create_liger([adata1, adata2])

    """
    # On-disk mode (set for online learning approach)
    if adata_list[0].isbacked:
        processed_list = []
        for adata in adata_list:
            processed_list.append(_initialization_online(adata, chunk_size, remove_missing))

        liger_object = Liger(processed_list)

    # In-memory mode
    else:
        liger_object = _create_liger_matrix(adata_list, make_sparse, take_gene_union, remove_missing)

    return liger_object


def _initialization_online(adata, chunk_size, remove_missing):
    """"""

    # calculate row sum and sum of squares using raw data
    gene_sum = np.zeros(adata.shape[1])
    gene_sum_sq = np.zeros(adata.shape[1])
    nUMI = np.zeros(adata.shape[0])
    nGene = np.zeros(adata.shape[0])
    for left, right in _h5_idx_generator(chunk_size, adata.shape[0]):
        raw_data = adata.X[left:right, :]
        gene_sum += np.ravel(np.sum(raw_data, axis=0))
        gene_sum_sq += np.ravel(np.sum(raw_data.power(2), axis=0))
        nUMI[left:right] = np.ravel(np.sum(raw_data, axis=1))
        nGene[left:right] = raw_data.getnnz(axis=1)

    import os
    # create results folder if not exists
    if not os.path.exists('results'):
        os.makedirs('results')

    # save AnnData
    file_path = './results/' + adata.uns['sample_name'] + '.h5ad'
    if remove_missing:
        idx_missing = gene_sum == 0
    else:
        idx_missing = np.repeat(False, adata.shape[1])
    processed_adata = adata[:, ~idx_missing].copy(file_path)
    processed_adata.var['gene_sum'] = gene_sum[~idx_missing]
    processed_adata.var['gene_sum_sq'] = gene_sum_sq[~idx_missing]
    processed_adata.obs['nUMI'] = nUMI
    processed_adata.obs['nGene'] = nGene

    return processed_adata


def _create_liger_matrix(adata_list, make_sparse, take_gene_union, remove_missing):
    """"""
    from scipy.sparse import csr_matrix

    num_samples = len(adata_list)

    # Make matrix sparse
    if make_sparse:
        for idx, adata in enumerate(adata_list):
            # force raw data to be csr matrix
            adata_list[idx].X = csr_matrix(adata_list[idx].X, dtype=np.int)
            # check if dimnames exist
            if not adata.obs.index.name or not adata.var.index.name:
                raise ValueError('Raw data must have both row (cell) and column (gene) names.')
            # check whether cell name is unique or not
            if not adata.obs.index.is_unique and adata.X.shape[0] > 1:
                raise ValueError('At least one cell name is repeated across datasets; '
                                 'please make sure all cell names are unique.')

    # Take gene union (requires make_sparse=True)
    if take_gene_union and make_sparse:
        merged_data = _merge_sparse_data_all(adata_list)
        if remove_missing:
            missing_genes = np.ravel(np.sum(merged_data.X, axis=0)) == 0
            if np.sum(missing_genes) > 0:
                print('Removing {} genes not expressed in any cells across merged datasets.'.format(
                    np.sum(missing_genes)))
                # show gene name when the total of missing genes is less than 25
                if np.sum(missing_genes) < 25:
                    print(merged_data.var.index[missing_genes])
                # save data after removing missing genes
                merged_data = merged_data[:, ~missing_genes].copy()
        # fill out raw_data matrices with union of genes across all datasets
        for i in range(num_samples):
            adata_list[i] = merged_data[merged_data.obs.index == adata_list[i].index, :].copy()

    # Remove missing cells
    for idx, adata in enumerate(adata_list):
        if remove_missing:
            adata_list[idx] = _remove_missing_obs(adata, use_rows=True)
            # remove missing genes if not already merged
            if not take_gene_union:
                adata_list[idx] = _remove_missing_obs(adata, use_rows=False)

    # Create liger object based on raw data list
    liger_object = Liger(adata_list)

    # Initialize cell_data for liger_object with nUMI, nGene, and dataset
    for idx, adata in enumerate(adata_list):
        liger_object.adata_list[idx].obs['nUMI'] = np.ravel(np.sum(adata.X, axis=1))
        liger_object.adata_list[idx].obs['nGene'] = adata.X.getnnz(axis=1)
            #np.count_nonzero(adata.X.toarray(), axis=1)
        liger_object.adata_list[idx].obs['dataset'] = np.repeat(adata.uns['sample_name'], adata.obs.index.shape[0])

    return liger_object


