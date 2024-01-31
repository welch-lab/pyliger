import h5sparse
import numpy as np
from scipy.sparse import vstack
from tqdm import tqdm

from pyliger._utilities import _h5_idx_generator
from pyliger.preprocessing._initialization import _initialization_online
from pyliger.preprocessing._normalization import _normalize_online
from pyliger.preprocessing._scale import _scale_online
from pyliger.factorization._utilities import (
    _init_V_online,
    _init_W,
    _update_V_HALS,
    _update_W_HALS,
    nnlsm_blockpivot,
    nonneg,
)


# from memory_profiler import profile
# @profile
def online_iNMF(
    liger_object,
    X_new=None,
    projection=False,
    W_init=None,
    V_init=None,
    H_init=None,
    A_init=None,
    B_init=None,
    k=20,
    value_lambda=5.0,
    max_epochs=5,
    miniBatch_max_iters=1,
    miniBatch_size=5000,
    h5_chunk_size=1000,
    rand_seed=1,
    verbose=True,
) -> None:
    """Perform online iNMF on scaled datasets

    Perform online integrative non-negative matrix factorization to represent multiple single-cell datasets
    in terms of H, W, and V matrices. It optimizes the iNMF objective function using online learning (non-negative
    least squares for H matrix, hierarchical alternating least squares for W and V matrices), where the
    number of factors is set by k. The function allows online learning in 3 scenarios: (1) fully observed datasets;
    (2) iterative refinement using continually arriving datasets; and (3) projection of new datasets without updating
    the existing factorization. All three scenarios require fixed memory independent of the number of cells.

    For each dataset, this factorization produces an H matrix (cells by k), a V matrix (k by genes),
    and a shared W matrix (k by genes). The H matrices represent the cell factor loadings.
    W is identical among all datasets, as it represents the shared components of the metagenes
    across datasets. The V matrices represent the dataset-specific components of the metagenes.

    Parameters
    ----------
        liger_object:
        X_new: AnnData object
        projection:
        W_init:
        V_init:
        H_init:
        A_init:
        B_init:
        k:
        value_lambda:
        max_epochs:
        miniBatch_max_iters:
        miniBatch_size:
        h5_chunk_size:
        rand_seed:
        verbose:

    Returns
    -------

    """
    matrices_init_dict = {
        "W_init": W_init,
        "V_init": V_init,
        "H_init": H_init,
        "A_init": A_init,
        "B_init": B_init,
    }

    factorization_params = {
        "k": k,
        "value_lambda": value_lambda,
        "max_epochs": max_epochs,
        "miniBatch_max_iters": miniBatch_max_iters,
        "miniBatch_size": miniBatch_size,
        "h5_chunk_size": h5_chunk_size,
        "rand_seed": rand_seed,
    }

    # Online Learning Scenario 1
    if X_new is None:
        if verbose:  # print start information
            print("Starting Online iNMF...")
        _online_iNMF_from_scratch(
            liger_object, verbose, matrices_init_dict, factorization_params
        )
    else:
        # check whether X_new needs to be processed
        for idx, adata in enumerate(X_new):
            if "scale_data" not in adata.layers.keys():
                X_new[idx] = _preprocessing(
                    adata, h5_chunk_size, liger_object.var_genes
                )

        # Online Learning Scenario 3
        if projection:
            if verbose:
                print("Metagene projection")
            _projection(liger_object, X_new, W_init, k, miniBatch_size)
        # Online Learning Scenario 2
        else:
            _online_iNMF_refine(
                liger_object, X_new, verbose, matrices_init_dict, factorization_params
            )

    return None


def _online_iNMF_from_scratch(
    liger_object, verbose, matrices_init_dict, factorization_params
):
    """ """
    ### 0. Extract required information
    # prepare basic dataset profiles
    num_files = liger_object.num_samples  # number of total input hdf5 files
    num_genes = len(liger_object.var_genes)  # number of variable genes
    num_cells = [
        adata.shape[0] for adata in liger_object.adata_list
    ]  # number of cells in each hdf5 files
    prev_info = None
    Xs = [_get_scale_data(adata) for adata in liger_object.adata_list]

    # import scipy.io
    # from scipy.sparse import csr_matrix
    # Xs = [csr_matrix(scipy.io.mmread('/Users/lulu/Documents/GitHub/pyliger/results/stim.mtx').transpose()), csr_matrix(scipy.io.mmread('/Users/lulu/Documents/GitHub/pyliger/results/ctrl.mtx').transpose())]
    # Xs = [csr_matrix(scipy.io.mmread('/Users/lulu/Desktop/cells.mtx').transpose())]

    ### 1. Run Online_iNMF
    W, Vs, A, B = _online_iNMF_cal_W_V(
        Xs,
        num_genes,
        num_cells,
        verbose,
        prev_info,
        **matrices_init_dict,
        **factorization_params,
    )

    Hs = _online_iNMF_cal_H(
        Xs,
        W,
        Vs,
        num_genes,
        num_cells,
        verbose,
        factorization_params["value_lambda"],
        factorization_params["miniBatch_size"],
    )

    ### 2. Sava results and close hdf5 files
    for file in Xs:  # close all files in the end
        if isinstance(file, h5sparse.h5sparse.File):
            file.close()

    for i in range(num_files):
        liger_object.adata_list[i].obsm["H"] = Hs[i].transpose()
        liger_object.adata_list[i].varm["W"] = W
        liger_object.adata_list[i].varm["V"] = Vs[i]
        liger_object.adata_list[i].varm["B"] = B[i]
        liger_object.adata_list[i].uns["A"] = A[i]

    return None


def _online_iNMF_refine(
    liger_object, X_new, verbose, matrices_init_dict, factorization_params
):
    """ """
    if verbose:
        print("{} new datasets detected.".format(len(X_new)))

    ### 0. Extract required information
    # prepare basic dataset profiles
    num_genes = len(liger_object.var_genes)  # number of variable genes

    if matrices_init_dict["W_init"] is None:
        matrices_init_dict["W_init"] = liger_object.W
    prev_info = {
        "A_prev": [adata.uns["A"] for adata in liger_object.adata_list],
        "B_prev": [adata.varm["B"] for adata in liger_object.adata_list],
        "V_prev": [adata.varm["V"] for adata in liger_object.adata_list],
    }

    for i, adata in enumerate(X_new):
        # open hdf5 files
        X = [_get_scale_data(adata)]

        # import scipy.io
        # from scipy.sparse import csr_matrix
        # X = [csr_matrix(scipy.io.mmread('/Users/lulu/Desktop/nuclei.mtx').transpose())]

        num_cells = [adata.shape[0]]

        # calculate W, V, H matrices for the new datasets
        W, Vs, A, B = _online_iNMF_cal_W_V(
            X,
            num_genes,
            num_cells,
            verbose,
            prev_info,
            **matrices_init_dict,
            **factorization_params,
        )

        Hs = _online_iNMF_cal_H(
            X,
            W,
            Vs,
            num_genes,
            num_cells,
            verbose,
            factorization_params["value_lambda"],
            factorization_params["miniBatch_size"],
        )

        # Sava results and close hdf5 files
        if isinstance(X[0], h5sparse.h5sparse.File):
            X[0].close()
        X_new[i].obsm["H"] = Hs[0].transpose()
        X_new[i].varm["W"] = W
        X_new[i].varm["V"] = Vs[0]
        X_new[i].varm["B"] = B[0]
        X_new[i].uns["A"] = A[0]

        matrices_init_dict["W_init"] = W

    # update H matrics for existing datasets
    for i, adata in enumerate(liger_object.adata_list):
        X = [_get_scale_data(adata)]
        # X = [csr_matrix(scipy.io.mmread('/Users/lulu/Desktop/cells.mtx').transpose())]
        Vs = [adata.varm["V"]]
        num_cells = [adata.shape[0]]
        Hs = _online_iNMF_cal_H(
            X,
            W,
            Vs,
            num_genes,
            num_cells,
            verbose,
            factorization_params["value_lambda"],
            factorization_params["miniBatch_size"],
        )

        if isinstance(X[0], h5sparse.h5sparse.File):
            X[0].close()
        liger_object.adata_list[i].obsm["H"] = Hs[0].transpose()
        liger_object.adata_list[i].varm["W"] = W

    for adata in X_new:
        liger_object.adata_list.append(adata)

    return None


def _projection(liger_object, X_new, W_init, k, miniBatch_size):
    """"""
    ### 0. Extract required information
    # prepare basic dataset profiles
    num_files = len(X_new)  # number of total input hdf5 files
    num_genes = len(liger_object.var_genes)  # number of variable genes
    num_cells = [
        adata.shape[0] for adata in X_new
    ]  # number of cells in each hdf5 files

    ### 1. Initialization (W, V_i, A(HiHi^t), B(XiHit))
    W = liger_object.W if W_init is None else W_init

    # open hdf5 files
    Xs = [_get_scale_data(adata) for adata in X_new]

    # import scipy.io
    # from scipy.sparse import csr_matrix
    # Xs = [csr_matrix(scipy.io.mmread('/Users/lulu/Desktop/nuclei.mtx').transpose())]

    Hs = []
    Vs = []
    for i in range(num_files):
        num_batch = np.ceil(num_cells[i] / miniBatch_size).astype(int)
        H_miniBatch = []
        for batch in range(num_batch):
            if batch != num_batch:
                X_miniBatch = (
                    Xs[i]["scale_data"][
                        batch * miniBatch_size : (batch + 1) * miniBatch_size
                    ]
                    .transpose()
                    .toarray()
                )
                # X_miniBatch = Xs[i][
                #              batch * miniBatch_size:(batch + 1) * miniBatch_size].transpose().toarray()
            else:
                X_miniBatch = (
                    Xs[i]["scale_data"][batch * miniBatch_size : num_cells[i]]
                    .transpose()
                    .toarray()
                )
                # X_miniBatch = Xs[i][batch * miniBatch_size:num_cells[i]].transpose().toarray()

            H_miniBatch.append(nnlsm_blockpivot(A=W, B=X_miniBatch)[0])

        Hs.append(np.hstack(H_miniBatch))
        Vs.append(np.zeros((num_genes, k)))

    # Sava results and close hdf5 files
    for i in range(num_files):
        if isinstance(Xs[i], h5sparse.h5sparse.File):
            Xs[i].close()
        X_new[i].obsm["H"] = Hs[i].transpose()
        X_new[i].varm["W"] = W
        X_new[i].varm["V"] = Vs[i]

        # add object into liger object TODO: think of another way of appending adata into liger object
        liger_object.adata_list.append(X_new[i])

    return None


def _online_iNMF_cal_W_V(
    Xs,
    num_genes,
    num_cells,
    verbose,
    prev_info,
    W_init,
    V_init,
    H_init,
    A_init,
    B_init,
    k,
    value_lambda,
    max_epochs,
    miniBatch_max_iters,
    miniBatch_size,
    h5_chunk_size,
    rand_seed,
):
    """ """
    ### 1. Initialization (W, V_i, A(HiHi^t), B(XiHit))
    num_files = len(Xs)

    np.random.seed(seed=rand_seed)
    W = _init_W(num_genes, k, rand_seed) if W_init is None else W_init
    Vs = (
        [
            _init_V_online(num_cells[i], k, Xs[i], h5_chunk_size, rand_seed)
            for i in range(num_files)
        ]
        if V_init is None
        else V_init
    )
    A = (
        [np.zeros((k, k)) for i in range(num_files)] if A_init is None else A_init
    )  # A = HiHi^t
    B = (
        [np.zeros((num_genes, k)) for i in range(num_files)]
        if B_init is None
        else B_init
    )  # B = XiHit

    # save information older than 2 epochs
    A_old = [np.zeros((k, k)) for i in range(num_files)]
    B_old = [np.zeros((num_genes, k)) for i in range(num_files)]

    ### 2. Iteration starts here
    # determine batch size for each dataset
    miniBatch_sizes = np.asarray(
        [
            np.round((num_cells[i] / np.sum(num_cells)) * miniBatch_size).astype(int)
            for i in range(num_files)
        ]
    )

    # determine total number of iterations
    total_iters = np.floor(num_cells[0] * max_epochs / miniBatch_sizes[0]).astype(int)

    # generate all slicing index for each iteration
    all_idx = [
        _generate_idx(
            total_iters, miniBatch_sizes[i], max_epochs, h5_chunk_size, num_cells[i]
        )
        for i in range(num_files)
    ]

    epoch = 0
    for iter in tqdm(range(total_iters)):
        # track epochs
        epoch_prev = epoch
        epoch = ((iter + 1) * miniBatch_sizes[0]) // num_cells[0]

        # cells from each dataset that are used in this iteration
        Xs_miniBatch = []
        for i in range(num_files):
            miniBatch_idx = all_idx[i][iter]
            X_miniBatch = vstack(
                [Xs[i]["scale_data"][left:right] for left, right in miniBatch_idx]
            )
            # X_miniBatch = vstack([Xs[i][left:right] for left, right in miniBatch_idx]) ###########
            Xs_miniBatch.append(X_miniBatch.transpose().toarray())

        ## 1) update H_i_minibatch by ANLS
        H_miniBatch = []
        for i in range(num_files):
            H_miniBatch.append(
                nnlsm_blockpivot(
                    A=np.vstack(((W + Vs[i]), np.sqrt(value_lambda) * Vs[i])),
                    B=np.vstack(
                        (Xs_miniBatch[i], np.zeros((num_genes, miniBatch_sizes[i])))
                    ),
                )[0]
            )

        ## 2) updata A and B matrices
        for i in range(num_files):
            A[i], B[i], A_old[i], B_old[i] = _update_A_B(
                A[i],
                B[i],
                A_old[i],
                B_old[i],
                H_miniBatch[i],
                Xs_miniBatch[i],
                miniBatch_sizes[i],
                iter,
                epoch,
                epoch_prev,
            )

        ## 3) update W, V_i by HALS
        iter_miniBatch = 1
        while iter_miniBatch <= miniBatch_max_iters:
            # update W
            if prev_info is None:
                W = _update_W_HALS(A, B, W, Vs)
            else:
                W = _update_W_HALS(
                    prev_info["A_prev"] + A,
                    prev_info["B_prev"] + B,
                    W,
                    prev_info["V_prev"] + Vs,
                )

            # update V_i
            Vs = _update_V_HALS(A, B, W, Vs, value_lambda)

            iter_miniBatch += 1

    return W, Vs, A, B


def _online_iNMF_cal_H(
    Xs, W, Vs, num_genes, num_cells, verbose, value_lambda, miniBatch_size
):
    """"""
    # Calculate metagene loadings (H metrics)
    num_files = len(Xs)
    if verbose:
        print("Calculate metagene loadings...")

    Hs = []
    for i in range(num_files):
        num_batch = np.ceil(num_cells[i] / miniBatch_size).astype(int)
        H_miniBatch = []
        for batch in range(num_batch):
            if batch != num_batch:
                X_miniBatch = (
                    Xs[i]["scale_data"][
                        batch * miniBatch_size : (batch + 1) * miniBatch_size
                    ]
                    .transpose()
                    .toarray()
                )
                # X_miniBatch = Xs[i][batch * miniBatch_size:(batch + 1) * miniBatch_size].transpose().toarray()

            else:
                X_miniBatch = (
                    Xs[i]["scale_data"][batch * miniBatch_size : num_cells[i]]
                    .transpose()
                    .toarray()
                )
                # X_miniBatch = Xs[i][batch * miniBatch_size:num_cells[i]].transpose().toarray()

            H_miniBatch.append(
                nnlsm_blockpivot(
                    A=np.vstack(((W + Vs[i]), np.sqrt(value_lambda) * Vs[i])),
                    B=np.vstack(
                        (X_miniBatch, np.zeros((num_genes, X_miniBatch.shape[1])))
                    ),
                )[0]
            )

        Hs.append(np.hstack(H_miniBatch))

    return Hs


def _generate_idx(num_iter, miniBatch_size, max_epochs, h5_chunk_size, num_cell):
    """helper function to generate miniBatch index for each iteration"""
    idx_dict = {}

    # permutate chunks for all+1 epochs (one more set for extreme case)
    all_idx = np.concatenate(
        [_chunk_permutation(num_cell, h5_chunk_size) for j in range(max_epochs + 1)]
    )

    # assign chunks for each iteration
    for i in range(num_iter):
        total_sum = 0
        temp_list = []
        while total_sum < miniBatch_size:
            left, right = all_idx[0]
            diff = right - left
            missing = int(miniBatch_size - total_sum)
            if diff < missing:
                temp_list.append((left, right))
                all_idx = all_idx[1:]
                total_sum += diff
            elif diff == missing:
                temp_list.append((left, right))
                all_idx = all_idx[1:]
                total_sum += diff
            else:
                temp_list.append((left, left + missing))
                all_idx[0] = [left + missing, right]
                total_sum += missing
        idx_dict[i] = temp_list

    return idx_dict


def _chunk_permutation(num_cell, h5_chunk_size):
    """helper function to permutate chunks"""
    # all_chunks = np.random.permutation([(left, right) for left, right in _h5_idx_generator(h5_chunk_size, num_cell)])
    all_chunks = np.asarray(
        [(left, right) for left, right in _h5_idx_generator(h5_chunk_size, num_cell)]
    )

    return all_chunks


def _get_scale_data(adata):
    """helper function to open hdf5 files if AnnData is on-disk, otherwise get scale data from AnnData object"""
    if adata.isbacked:
        return h5sparse.File("./results/" + adata.uns["sample_name"] + ".hdf5", "r")
    else:
        return adata.layers


def _update_A_B(
    A,
    B,
    A_old,
    B_old,
    H_miniBatch,
    X_miniBatch,
    miniBatch_size,
    iter,
    epoch,
    epoch_prev,
):
    """helper function to update A B matrices"""
    if iter == 0:
        scale_param = 0
    elif iter == 1:
        scale_param = 1 / miniBatch_size
    else:
        scale_param = (iter - 1) / iter

    if epoch > 0 and epoch - epoch_prev == 1:  # remove information older than 2 epochs
        A = A - A_old
        A_old = scale_param * A
        B = B - B_old
        B_old = scale_param * B
    else:  # otherwise scale the old information
        A_old = scale_param * A_old
        B_old = scale_param * B_old

    t_H_miniBatch = H_miniBatch.transpose()
    A = scale_param * A + (H_miniBatch @ t_H_miniBatch) / miniBatch_size
    B = scale_param * B + (X_miniBatch @ t_H_miniBatch) / miniBatch_size

    # replace the zero values in the diagonal of A with 1e-15
    A_diag = np.diagonal(A).copy()
    A_diag[A_diag == 0] = 1e-15
    np.fill_diagonal(A, A_diag)

    return A, B, A_old, B_old


def _preprocessing(adata, h5_chunk_size, var_genes):
    """helper function to preprocess AnnData object
    TODO: in the future, consider have preprocessing support for AnnData objects"""
    # initialization
    adata = _initialization_online(adata, h5_chunk_size, remove_missing=False)
    # normalization
    norm_sum, norm_sum_sq = _normalize_online(adata, h5_chunk_size)
    norm_sum = nonneg(norm_sum)  # to avoid zero-division
    norm_sum_sq = nonneg(norm_sum_sq)  # to avoid zero-division
    adata.var["norm_sum"] = norm_sum
    adata.var["norm_sum_sq"] = norm_sum_sq
    adata.var["norm_mean"] = norm_sum / adata.shape[0]
    # scale
    var_gene_idx = adata.var.index.isin(var_genes)
    adata = _scale_online(adata, var_gene_idx, h5_chunk_size)

    return adata


"""
    ### 1. Initialization (W, V_i, A(HiHi^t), B(XiHit))
    np.random.seed(seed=rand_seed)
    W = _init_W(num_genes, k, rand_seed) if W_init is None else W_init
    Vs = [_init_V_online(num_cells[i], k, Xs[i], h5_chunk_size, rand_seed) for i in range(num_files)] if V_init is None else V_init
    A = [np.zeros((k, k)) for i in range(num_files)] if A_init is None else A_init # A = HiHi^t
    B = [np.zeros((num_genes, k)) for i in range(num_files)] if B_init is None else B_init # B = XiHit

    # save information older than 2 epochs
    A_old = [np.zeros((k, k)) for i in range(num_files)]
    B_old = [np.zeros((num_genes, k)) for i in range(num_files)]

    ### 2. Iteration starts here
    # determine batch size for each dataset
    miniBatch_sizes = np.asarray([np.round((num_cells[i] / np.sum(num_cells)) * miniBatch_size).astype(int) for i in range(num_files)])

    # determine total number of iterations
    total_iters = np.floor(num_cells[0] * max_epochs / miniBatch_sizes[0]).astype(int)

    # generate all slicing index for each iteration
    all_idx = [_generate_idx(total_iters, miniBatch_sizes[i], max_epochs, h5_chunk_size, num_cells[i]) for i in range(num_files)]

    epoch = 0
    for iter in tqdm(range(total_iters)):
        # track epochs
        epoch_prev = epoch
        epoch = ((iter+1) * miniBatch_sizes[0]) // num_cells[0]

        # cells from each dataset that are used in this iteration
        Xs_miniBatch = []
        for i in range(num_files):
            miniBatch_idx = all_idx[i][iter]
            X_miniBatch = vstack([Xs[i]['scale_data'][left:right] for left, right in miniBatch_idx])
            #X_miniBatch = vstack([Xs[i][left:right] for left, right in miniBatch_idx]) ###########
            Xs_miniBatch.append(X_miniBatch.transpose().toarray())

        ## 1) update H_i_minibatch by ANLS
        H_miniBatch = []
        for i in range(num_files):
            H_miniBatch.append(nnlsm_blockpivot(A=np.vstack(((W + Vs[i]), np.sqrt(value_lambda) * Vs[i])),
                                                B=np.vstack((Xs_miniBatch[i], np.zeros((num_genes, miniBatch_sizes[i])))))[0])

        ## 2) updata A and B matrices
        for i in range(num_files):
            A[i], B[i], A_old[i], B_old[i] = _update_A_B(A[i], B[i], A_old[i], B_old[i], H_miniBatch[i], Xs_miniBatch[i], miniBatch_sizes[i], iter, epoch, epoch_prev)

        ## 3) update W, V_i by HALS
        iter_miniBatch = 1
        while iter_miniBatch <= miniBatch_max_iters:
            W = _update_W_HALS(A, B, W, Vs)  # update W
            Vs = _update_V_HALS(A, B, W, Vs, value_lambda)  # update V_i

            iter_miniBatch += 1

    # Calculate metagene loadings (H metrics)
    if verbose:
        print('Calculate metagene loadings...')

    Hs = []
    for i in range(num_files):
        num_batch = np.ceil(num_cells[i] / miniBatch_size).astype(int)
        H_miniBatch = []
        for batch in range(num_batch):
            if batch != num_batch:
                X_miniBatch = Xs[i]['scale_data'][batch*miniBatch_size:(batch+1)*miniBatch_size].transpose().toarray()
                #X_miniBatch = Xs[i][batch * miniBatch_size:(batch + 1) * miniBatch_size].transpose().toarray()
            else:
                X_miniBatch = Xs[i]['scale_data'][batch*miniBatch_size:num_cells[i]].transpose().toarray()
                #X_miniBatch = Xs[i][batch * miniBatch_size:num_cells[i]].transpose().toarray()

            H_miniBatch.append(nnlsm_blockpivot(A=np.vstack(((W + Vs[i]), np.sqrt(value_lambda) * Vs[i])),
                                                B=np.vstack((X_miniBatch, np.zeros((num_genes, X_miniBatch.shape[1])))))[0])

        Hs.append(np.hstack(H_miniBatch))


        ### 1. Initialization (W, V_i, A(HiHi^t), B(XiHit))
        np.random.seed(seed=rand_seed)
        W = _init_W(num_genes, k, rand_seed) if W_init is None else W_init
        V = _init_V_online(num_cells, k, X, h5_chunk_size, rand_seed) if V_init is None else V_init
        A = np.zeros((k, k)) if A_init is None else A_init  # A = HiHi^t
        B = np.zeros((num_genes, k)) if B_init is None else B_init  # B = XiHit

        # save information older than 2 epochs
        A_old = np.zeros((k, k))
        B_old = np.zeros((num_genes, k))

        # determine total number of iterations
        total_iters = np.floor(num_cells[0] * max_epochs / miniBatch_size).astype(int)

        # generate all slicing index for each iteration
        all_idx = _generate_idx(total_iters, miniBatch_size, max_epochs, h5_chunk_size, num_cells[i])

        ### 2. Iteration starts here
        epoch = 0
        for iter in tqdm(range(total_iters)):
            # track epochs
            epoch_prev = epoch
            epoch = ((iter+1) * miniBatch_size) // num_cells[i]

            # cells that are used in this iteration
            miniBatch_idx = all_idx[iter]
            X_miniBatch = vstack([X['scale_data'][left:right] for left, right in miniBatch_idx])
            # X_miniBatch = vstack([X[left:right] for left, right in miniBatch_idx]) ###########
            X_miniBatch = X_miniBatch.transpose().toarray()

            ## 1) update H_i_minibatch by ANLS
            H_miniBatch = nnlsm_blockpivot(A=np.vstack(((W + V), np.sqrt(value_lambda) * V)),
                                           B=np.vstack((X_miniBatch, np.zeros((num_genes, miniBatch_size)))))[0]

            ## 2) updata A and B matrices
            A, B, A_old, B_old = _update_A_B(A, B, A_old, B_old, H_miniBatch, X_miniBatch, miniBatch_size, iter, epoch, epoch_prev)

            ## 3) update W, V_i by HALS
            iter_miniBatch = 1
            while iter_miniBatch <= miniBatch_max_iters:
                W = _update_W_HALS([A], [B], W, [V])  # update W
                V = _update_V_HALS([A], [B], W, [V], value_lambda)[0]  # update V

                iter_miniBatch += 1


        # Calculate metagene loadings (H metrics)
        if verbose:
            print('Calculate metagene loadings...')

        Hs = []
        for i in range(num_files):
            num_batch = np.ceil(num_cells[i] / miniBatch_size).astype(int)
            H_miniBatch = []
            for batch in range(num_batch):
                if batch != num_batch:
                    X_miniBatch = Xs[i]['scale_data'][batch*miniBatch_size:(batch+1)*miniBatch_size].transpose().toarray()
                    #X_miniBatch = Xs[i][
                    #              batch * miniBatch_size:(batch + 1) * miniBatch_size].transpose().toarray()
                else:
                    X_miniBatch = Xs[i]['scale_data'][batch*miniBatch_size:num_cells[i]].transpose().toarray()
                    #X_miniBatch = Xs[i][batch * miniBatch_size:num_cells[i]].transpose().toarray()

                H_miniBatch.append(nnlsm_blockpivot(A=np.vstack(((W + Vs[i]), np.sqrt(value_lambda) * Vs[i])),
                                                    B=np.vstack((X_miniBatch, np.zeros((num_genes, X_miniBatch.shape[1])))))[0])

            Hs.append(np.hstack(H_miniBatch))

                    # Close all files in the end
        X.close()

        # Sava results and close hdf5 files
        X_new[i].obsm['H'] = Hs[i].transpose()
        X_new[i].varm['W'] = W
        X_new[i].varm['V'] = V

"""
