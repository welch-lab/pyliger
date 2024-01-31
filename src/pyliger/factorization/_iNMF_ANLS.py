import numpy as np
from tqdm import tqdm

from pyliger.factorization._utilities import nnlsm_blockpivot


def optimize_ALS(
    liger_object,
    k,
    value_lambda=5.0,
    thresh=1e-6,
    max_iters=30,
    nrep=1,
    H_init=None,
    W_init=None,
    V_init=None,
    rand_seed=1,
    print_obj=False,
):
    """Perform iNMF on scaled datasets

    Perform integrative non-negative matrix factorization to return factorized H, W, and V matrices.
    It optimizes the iNMF objective function using block coordinate descent (alternating non-negative
    least squares), where the number of factors is set by k. TODO: include objective function equation here in documentation (using deqn)

    For each dataset, this factorization produces an H matrix (cells by k), a V matrix (k by genes),
    and a shared W matrix (k by genes). The H matrices represent the cell factor loadings.
    W is held consistent among all datasets, as it represents the shared components of the metagenes
    across datasets. The V matrices represent the dataset-specific components of the metagenes.

    Args:
        liger_object(liger):
            Should normalize, select genes, and scale before calling.
        k(int):
            Inner dimension of factorization (number of factors). Run suggestK to determine
            appropriate value; a general rule of thumb is that a higher k will be needed for datasets with
            more sub-structure.
        value_lambda(float): optional, (default 5.0)
            Regularization parameter. Larger values penalize dataset-specific effects more
            strongly (ie. alignment should increase as lambda increases). Run suggestLambda to determine
            most appropriate value for balancing dataset alignment and agreement.
        thresh(float): optional, (default 1e-6)
            Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh.
        max_iters(int): optional, (default 30)
            Maximum number of block coordinate descent iterations to perform
        nrep(int): optional, (default 1)
            Number of restarts to perform (iNMF objective function is non-convex, so taking the
            best objective from multiple successive initializations is recommended). For easier
            reproducibility, this increments the random seed by 1 for each consecutive restart, so future
            factorizations of the same dataset can be run with one rep if necessary.
        H_init(): optional, (default None)
            Initial values to use for H matrices.
        W_init(): optional, (default None)
            Initial values to use for W matrix.

        V_init(): optional, (default None)
            Initial values to use for V matrices.
        rand_seed(seed): optional, (default 1)
            Random seed to allow reproducible results
        print_obj(bool): optional, (default False)
            Print objective function values after convergence.

    Return:
        liger_object(liger):
            liger object with H, W, and V attributes.

    Usage:
        >>> adata1 = AnnData(np.arange(12).reshape((4, 3)))
        >>> adata2 = AnnData(np.arange(12).reshape((4, 3)))
        >>> ligerex = pyliger.create_liger([adata1, adata2])
        >>> ligerex = pyliger.normalize(ligerex)
        >>> ligerex = pyliger.select_genes(ligerex) # select genes
        >>> ligerex = pyliger.scale_not_center(ligerex)
        >>> ligerex = pyliger.optimize_ALS(ligerex, k = 20, value_lambda = 5, nrep = 3) # get factorization using three restarts and 20 factors
    """
    ### 0. Extract required information
    # prepare basic dataset profiles
    N = liger_object.num_samples  # number of total input hdf5 files
    ns = [
        adata.shape[0] for adata in liger_object.adata_list
    ]  # number of cells in each hdf5 files
    num_genes = len(liger_object.var_genes)  # number of variable genes
    X = [adata.layers["scale_data"].toarray() for adata in liger_object.adata_list]

    if k >= np.min(ns):
        raise ValueError(
            "Select k lower than the number of cells in smallest dataset: {}".format(
                np.min(ns)
            )
        )

    if k >= num_genes:
        raise ValueError(
            "Select k lower than the number of variable genes: {}".format(
                len(liger_object.var_genes)
            )
        )

    best_obj = np.Inf

    for j in range(nrep):
        np.random.seed(seed=rand_seed + j - 1)

        ### 1. Initialization (W, V_i, H_i)
        W = np.abs(np.random.uniform(0, 2, (k, num_genes)))
        V = [np.abs(np.random.uniform(0, 2, (k, num_genes))) for i in range(N)]
        H = [np.abs(np.random.uniform(0, 2, (ns[i], k))) for i in range(N)]

        if W_init is not None:
            W = W_init

        if V_init is not None:
            V = V_init

        if H_init is not None:
            H = H_init

        delta = 1
        sqrt_lambda = np.sqrt(value_lambda)

        # Initial training obj
        obj_train_approximation = 0
        obj_train_penalty = 0
        for i in range(N):
            obj_train_approximation += np.linalg.norm(X[i] - H[i] @ (W + V[i])) ** 2
            obj_train_penalty += np.linalg.norm(H[i] @ V[i]) ** 2

        obj0 = obj_train_approximation + value_lambda * obj_train_penalty

        ### 2. Iteration starts here
        for iter in tqdm(range(max_iters)):
            if delta > thresh:
                ## 1) update H matrix
                for i in range(N):
                    H[i] = nnlsm_blockpivot(
                        A=np.hstack(((W + V[i]), sqrt_lambda * V[i])).transpose(),
                        B=np.hstack((X[i], np.zeros((ns[i], num_genes)))).transpose(),
                    )[0].transpose()

                ## 2) update V matrix
                for i in range(N):
                    V[i] = nnlsm_blockpivot(
                        A=np.vstack((H[i], sqrt_lambda * H[i])),
                        B=np.vstack(((X[i] - H[i] @ W), np.zeros((ns[i], num_genes)))),
                    )[0]

                ## 3) update W matrix
                W = nnlsm_blockpivot(
                    A=np.vstack(H),
                    B=np.vstack([(X[i] - H[i] @ V[i]) for i in range(N)]),
                )[0]

                obj_train_prev = obj0
                obj_train_approximation = 0
                obj_train_penalty = 0
                for i in range(N):
                    obj_train_approximation += (
                        np.linalg.norm(X[i] - H[i] @ (W + V[i])) ** 2
                    )
                    obj_train_penalty += np.linalg.norm(H[i] @ V[i]) ** 2
                obj0 = obj_train_approximation + value_lambda * obj_train_penalty
                delta = np.absolute(obj_train_prev - obj0) / (
                    (obj_train_prev + obj0) / 2
                )
            else:
                continue

        if obj0 < best_obj:
            final_W = W
            final_H = H
            final_V = V
            best_obj = obj0
            best_seed = rand_seed + i - 1

        if print_obj:
            print("Objective: {}".format(best_obj))

    # liger_object.W = final_W.transpose()

    ### 3. Save results into the liger_object
    for i in range(N):
        liger_object.adata_list[i].obsm["H"] = final_H[i]
        liger_object.adata_list[i].varm["W"] = final_W.transpose()
        liger_object.adata_list[i].varm["V"] = final_V[i].transpose()
        # idx = liger_object.adata_list[i].uns['var_gene_idx']
        # shape = liger_object.adata_list[i].shape
        # save_W = np.zeros((shape[1], k))
        # save_W[idx, :] = final_W.transpose()
        # save_V = np.zeros((shape[1], k))
        # save_V[idx, :] = final_V[i].transpose()
        # liger_object.adata_list[i].obsm['H'] = final_H[i]
        # liger_object.adata_list[i].varm['W'] = save_W
        # liger_object.adata_list[i].varm['V'] = save_V

    return None
