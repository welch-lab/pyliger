from typing import Optional

import h5sparse
import numpy as np
from scipy.sparse import csr_matrix
from sklearn.utils.sparsefuncs import inplace_column_scale

from pyliger._utilities import _h5_idx_generator


def scale_not_center(
    liger_object, remove_missing=True, chunk_size: Optional[int] = 1000
) -> None:
    """Scale genes by root-mean-square across cells

    This function scales normalized gene expression data after variable genes have been selected.
    Note that the data is not mean-centered before scaling because expression values must remain
    positive (NMF only accepts positive values). It also removes cells which do not have any
    expression across the genes selected, by default.

    Parameters
    ----------
    liger_object : liger object
        Should call normalize and selectGenes before calling.
    remove_missing : bool, optional
        Whether to remove cells from scale_data with no gene expression
        (the default is True).
    chunk_size :

    Returns
    -------
    liger_object : liger object
        Object with scale_data layer.

    Examples
    --------
    >>> adata1 = AnnData(np.arange(12).reshape((4, 3)))
    >>> adata2 = AnnData(np.arange(12).reshape((4, 3)))
    >>> ligerex = create_liger([adata1, adata2])
    >>> ligerex = normalize(ligerex)
    >>> ligerex = select_genes(ligerex) # select genes
    >>> ligerex = scale_not_center(ligerex)
    """

    for idx, adata in enumerate(liger_object.adata_list):
        var_gene_idx = adata.var.index.isin(liger_object.var_genes)

        # On-disk mode (set for online learning approach)
        if adata.isbacked:
            liger_object.adata_list[idx] = _scale_online(
                adata, var_gene_idx, chunk_size
            )

        # In-memory mode
        else:
            liger_object.adata_list[idx] = _scale_matrix(adata, var_gene_idx)

    return None


def _scale_online(adata, var_gene_idx, chunk_size):
    file_path = "./results/" + adata.uns["sample_name"] + ".hdf5"
    with h5sparse.File(file_path, "r+") as f:
        for left, right in _h5_idx_generator(chunk_size, adata.shape[0]):
            scale_data = csr_matrix(
                f["norm_data"][left:right][:, var_gene_idx]
                / np.sqrt(
                    adata.var["norm_sum_sq"][var_gene_idx].to_numpy()
                    / (adata.shape[0] - 1)
                ),
                dtype=np.float64,
            )
            if "scale_data" not in f.keys():
                f.create_dataset(
                    "scale_data",
                    data=scale_data,
                    chunks=(chunk_size,),
                    maxshape=(None,),
                )
            else:
                f["scale_data"].append(scale_data)

    # slice adata after keeping a raw version
    # adata.raw = adata
    file_name = "./results/" + adata.uns["sample_name"] + ".h5ad"
    adata = adata[:, var_gene_idx].copy(filename=file_name)

    return adata


def _scale_matrix(adata, var_gene_idx):
    """"""
    # slice adata after keeping a raw version
    adata.raw = adata
    adata = adata[:, var_gene_idx].copy()

    # calculate scale data
    scale_data = adata.copy().layers["norm_data"]
    scaler = 1 / np.sqrt(adata.var["norm_sum_sq"].to_numpy() / (adata.shape[0] - 1))
    inplace_column_scale(scale_data, scaler)
    adata.layers["scale_data"] = scale_data

    # adata.layers['scale_data'] = csr_matrix(
    #    adata.layers['norm_data'] / np.sqrt(adata.var['norm_sum_sq'].to_numpy() / (adata.shape[0] - 1)),
    #    dtype=np.float64)
    # scale_data = csr_matrix(
    #    selected_data / np.sqrt(np.sum(np.square(selected_data.toarray()), axis=0) / (selected_data.shape[0] - 1)),
    #    dtype=np.float64)

    # numerical_idx = np.nonzero(var_gene_idx)[0]
    # row_idx, col_idx = scale_data.nonzero()
    # col_idx = np.take(numerical_idx, col_idx)
    # adata.layers['scale_data'] = csr_matrix((scale_data.data, (row_idx, col_idx)),
    #                                        shape=sample_shape, dtype=np.float64)

    # liger_object.adata_list[i] = liger_object.adata_list[i][:, idx].copy()
    # temp_norm = liger_object.adata_list[i].layers['norm_data']
    # liger_object.adata_list[i].layers['scale_data'] = csr_matrix(temp_norm / np.sqrt(np.sum(np.square(temp_norm.toarray()), axis=0) / (temp_norm.shape[0] - 1)))

    # if remove_missing:
    #    liger_object = _remove_missing_obs(liger_object, slot_use='scale_data', use_rows=False)
    return adata
