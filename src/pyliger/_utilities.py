import os

import h5sparse
import numpy as np
from anndata import AnnData


def _merge_sparse_data_all(adata_list, library_names=None):
    """Function to merge all sparse data into a single one

    Function takes in a list of DGEs, with gene row names and cell column names,
    and merges them into a single DGE.
    Also adds library_names to cell_names if expected to be overlap (common with 10X barcodes)

    Parameters
    ----------
    adata_list : AnnData
        List of AnnData objects which store expression matrices (gene by cell).
    library_names : list, optional
        (the default is None)

    Returns
    -------
    merged_adata : AnnData
        AnnData object stores expression matrix across datasets
    """

    merged_adata = AnnData()

    for adata in adata_list:
        merged_adata = merged_adata.concatenate(adata, join="inner")

    return merged_adata


def _remove_missing_obs(adata, slot_use="raw_data", use_rows=True):
    """Remove cells/genes with no expression across any genes/cells

    Removes cells/genes from chosen slot with no expression in any genes or cells respectively.

    Parameters
    ----------
    TODO: change to adata
    adata : AnnData object
        object (scale_data or norm_data must be set).
    slot_use : str, optional, 'raw_data' or 'scale_data'
        The data slot to filter (the default is 'raw_data').
    use_rows : bool, optional
        Treat each row as a cell (the default is True).

    Returns
    -------
    liger_object : liger object
        object with modified raw_data (or chosen slot) (dataset names preserved).

    Examples
    --------
    >>> adata = _remove_missing_obs(adata)
    """
    removed = str(
        np.where(
            slot_use in ["raw_data", "norm_data"] and use_rows == True, "cells", "genes"
        )
    )
    expressed = str(np.where(removed == "cells", " any genes", ""))

    data_type = adata.uns["sample_name"]
    if slot_use == "raw_data":
        filter_data = adata.X
    elif slot_use == "scale_data":
        filter_data = adata.layers["scale_data"]

    if use_rows:
        missing = np.ravel(np.sum(filter_data, axis=1)) == 0
    else:
        missing = np.ravel(np.sum(filter_data, axis=0)) == 0
    if np.sum(missing) > 0:
        # logging.info('Removing {} {} not expressing{} in {}.'.format(np.sum(missing), removed, expressed, data_type))
        print(
            "Removing {} {} not expressing{} in {}.".format(
                np.sum(missing), removed, expressed, data_type
            )
        )
        if use_rows:
            # show gene name when the total of missing is less than 25
            if np.sum(missing) < 25:
                print(adata.obs.index[missing])
            adata = adata[~missing, :].copy()
        else:
            # show cell name when the total of missing is less than 25
            if np.sum(missing) < 25:
                print(adata.var.index[missing])
            adata = adata[:, ~missing].copy()

    return adata


################################## For Use of hdf5 ################################
def _h5_idx_generator(chunk_size, matrix_size):
    """ """
    previous_idx = 0
    if matrix_size < chunk_size:
        current_idx = matrix_size
    else:
        current_idx = chunk_size
    num_chunk = np.ceil(matrix_size / chunk_size).astype(int)
    iter = 0
    while current_idx <= matrix_size and iter < num_chunk:
        yield int(previous_idx), int(current_idx)
        previous_idx += chunk_size
        current_idx += chunk_size
        if current_idx > matrix_size:
            current_idx = matrix_size
        iter += 1
    return None


def merge_H5(
    file_list,
    library_names,
    new_filename,
    format_type="10X",
    data_name=None,
    genes_name=None,
    barcodes_name=None,
):
    return None


def nonneg(x, eps=1e-16):
    """Given a input matrix, set all negative values to be zero"""
    x[x < eps] = eps
    return x


def _create_h5_using_adata(adata, chunk_size):
    # create h5 file.
    if not os.path.isdir("./results"):
        os.mkdir("./results")

    file_name = "./results/" + adata.uns["sample_name"] + ".hdf5"
    with h5sparse.File(file_name, "w") as f:
        for left, right in _h5_idx_generator(chunk_size, adata.shape[0]):
            if "raw_data" not in f.keys():
                f.create_dataset(
                    "raw_data",
                    data=adata[left:right, :].X,
                    chunks=(chunk_size,),
                    maxshape=(None,),
                )
            else:
                f["raw_data"].append(adata[left:right, :].X)
    return None
