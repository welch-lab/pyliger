from pathlib import Path
from typing import Optional

import h5sparse
import numpy as np
from sklearn.preprocessing import normalize as sp_normalize

from pyliger._utilities import _h5_idx_generator, _remove_missing_obs
from pyliger.pyliger import Liger

PARRENT_PATH = Path(__file__).parent


def normalize(
    liger_object: Liger, remove_missing: bool = True, chunk_size: Optional[int] = 1000
) -> None:
    """Normalize raw datasets to row sums

    This function normalizes data to account for total gene expression across a cell.

    Parameters
    ----------
    liger_object : liger object
        liger object with raw_data.
    chunk_size :

    remove_missing :

    Returns
    -------
    liger_object : liger object
        liger object with norm_data.

    Examples
    --------
    >>> adata1 = AnnData(np.arange(12).reshape((4, 3)))
    >>> adata2 = AnnData(np.arange(12).reshape((4, 3)))
    >>> ligerex = create_liger([adata1, adata2])
    >>> ligerex = normalize(ligerex)
    """
    # Iterate through each sample
    for idx, adata in enumerate(liger_object.adata_list):
        # On-disk mode (set for online learning approach)
        if adata.isbacked:
            norm_sum, norm_sum_sq = _normalize_online(adata, chunk_size)

        # In-memory mode
        else:
            norm_data, norm_sum, norm_sum_sq = _normalize_matrix(adata, remove_missing)
            liger_object.adata_list[idx].layers["norm_data"] = norm_data

        # save row sum and sum of squares for further use
        liger_object.adata_list[idx].var["norm_sum"] = norm_sum
        liger_object.adata_list[idx].var["norm_sum_sq"] = norm_sum_sq
        liger_object.adata_list[idx].var["norm_mean"] = norm_sum / adata.shape[0]

    return None


def _normalize_online(adata, chunk_size):
    """"""
    norm_sum = np.zeros(adata.shape[1])
    norm_sum_sq = np.zeros(adata.shape[1])

    # create h5 file for each individual sample.
    file_path = "./results/" + adata.uns["sample_name"] + ".hdf5"
    with h5sparse.File(file_path, "r+") as f:
        for left, right in _h5_idx_generator(chunk_size, adata.shape[0]):
            # normalize data and store normalized data as sparse matrix in h5 file
            # norm_data = sp_normalize(adata.X[left:right, :], axis=1, norm='l1')
            norm_data = sp_normalize(
                f["raw_data"][left:right][:, ~adata.uns["idx_missing"]],
                axis=1,
                norm="l1",
            )
            if "norm_data" not in f.keys():
                f.create_dataset(
                    "norm_data", data=norm_data, chunks=(chunk_size,), maxshape=(None,)
                )
            else:
                f["norm_data"].append(norm_data)

            # calculate row sum and sum of squares using normalized data
            norm_sum = norm_sum + np.ravel(np.sum(norm_data, axis=0))
            norm_sum_sq = norm_sum_sq + np.ravel(np.sum(norm_data.power(2), axis=0))

    return norm_sum, norm_sum_sq


def _normalize_matrix(adata, remove_missing):
    """"""
    if remove_missing:
        adata = _remove_missing_obs(adata, slot_use="raw_data", use_rows=True)
    norm_data = sp_normalize(adata.X, axis=1, norm="l1")
    norm_sum = np.ravel(np.sum(norm_data, axis=0))
    norm_sum_sq = np.ravel(np.sum(norm_data.power(2), axis=0))

    return norm_data, norm_sum, norm_sum_sq
