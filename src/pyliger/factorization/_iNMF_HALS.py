import time

import numpy as np

from pyliger.factorization._utilities import (
    _init_H,
    _init_V,
    _init_W,
    _update_H_HALS,
    _update_V_HALS,
    _update_W_HALS,
)


def iNMF_HALS(
    liger_object, k=20, value_lambda=5.0, thresh=1e-4, max_iters=25, nrep=1, rand_seed=1
):
    """Perform iNMF on scaled datasets using HALS method

    Parameters
    ----------
    liger_object : liger object
        Should normalize, select genes, and scale before calling.
    k : int, optional
        Inner dimension of factorization (number of factors). Run suggestK to determine
        appropriate value; a general rule of thumb is that a higher k will be needed for datasets with
        more sub-structure (the default is 20).
    value_lambda : float, optional
        Regularization parameter. Larger values penalize dataset-specific effects more
        strongly (ie. alignment should increase as lambda increases). Run suggestLambda to determine
        most appropriate value for balancing dataset alignment and agreement (the default is 5.0).
    thresh : float, optional
        Convergence threshold. Convergence occurs when |obj0-obj|/(mean(obj0,obj)) < thresh
        (the default is 1e-4).
    max_iters : int, optional
        Maximum number of block coordinate descent iterations to perform (the default is 25).
    nrep : int, optional
        Number of restarts to perform (iNMF objective function is non-convex, so taking the
        best objective from multiple successive initializations is recommended). For easier
        reproducibility, this increments the random seed by 1 for each consecutive restart, so future
        factorizations of the same dataset can be run with one rep if necessary (the default is 1).
    rand_seed : int, optional
        Random seed to allow reproducible results (the default is 1).

    Returns
    -------
    liger_object : liger object
        liger object with H, W, and V annotations.

    Examples
    --------
    >>> adata1 = AnnData(np.arange(12).reshape((4, 3)))
    >>> adata2 = AnnData(np.arange(12).reshape((4, 3)))
    >>> ligerex = pyliger.create_liger([adata1, adata2])
    >>> pyliger.normalize(ligerex)
    >>> pyliger.select_genes(ligerex) # select genes
    >>> pyliger.scale_not_center(ligerex)
    >>> pyliger.iNMF_HALS(ligerex, k = 20, value_lambda = 5, nrep = 3)
    """

    num_samples = len(liger_object.adata_list)
    X = [adata.layers["scale_data"].transpose() for adata in liger_object.adata_list]

    num_cells = [X[i].shape[1] for i in range(num_samples)]
    num_genes = X[0].shape[0]

    # set seed
    np.random.seed(seed=rand_seed)

    # allow restarts
    num_rep = 0
    while num_rep < nrep:
        # Initialization
        V = _init_V(num_cells, num_samples, k, X)
        W = _init_W(num_genes, k, rand_seed)
        H = _init_H(num_cells, num_samples, k)

        # Initial training obj
        obj_train_approximation = 0
        obj_train_penalty = 0
        for i in range(num_samples):
            obj_train_approximation += (
                np.linalg.norm(X[i].toarray() - (W + V[i]) @ H[i]) ** 2
            )
            obj_train_penalty += np.linalg.norm(V[i] @ H[i]) ** 2

        obj_train = obj_train_approximation + value_lambda * obj_train_penalty

        # print start information
        print("Initial Training Obj: {}".format(obj_train))

        iteration = 1
        total_time = 0  # track the total amount of time used for learning
        delta = np.inf

        # Perform block coordinate descent
        while delta > thresh and iteration <= max_iters:
            iter_start_time = time.time()

            # update H matrix
            H = _update_H_HALS(H, V, W, X, value_lambda)
            # VitVi = [V[i].transpose() @ V[i] for i in range(num_samples)]
            # W_Vi = [W + V[i] for i in range(num_samples)]
            # W_Vi_sq = [W_Vi[i].transpose() @ W_Vi[i] for i in range(num_samples)]
            # for i in range(num_samples):
            #    for j in range(k):
            #        H[i][j, :] = nonneg(H[i][j, :] + (
            #                    W_Vi[i][:, j].transpose() @ X[i] - W_Vi[i][:, j].transpose() @ W_Vi[i] @ H[
            #                i] - value_lambda * VitVi[i][j, :] @ H[i]) / (
            #                                        W_Vi_sq[i][j, j] + value_lambda * VitVi[i][j, j]))

            # update W matrix
            XHt = [X[i] @ H[i].transpose() for i in range(num_samples)]
            HHt = [H[i] @ H[i].transpose() for i in range(num_samples)]
            # A = List()
            # [A.append(x) for x in HHt]
            # B = List()
            # [B.append(x) for x in XHt]
            # temp_V = List()
            # [temp_V.append(x) for x in V]
            W = _update_W_HALS(HHt, XHt, W, V)

            # update V matrix
            V = _update_V_HALS(HHt, XHt, W, V, value_lambda)

            # training obj
            obj_train_prev = obj_train
            obj_train_approximation = 0
            obj_train_penalty = 0
            for i in range(num_samples):
                obj_train_approximation += np.linalg.norm(X[i] - (W + V[i]) @ H[i]) ** 2
                obj_train_penalty += np.linalg.norm(V[i] @ H[i]) ** 2
            obj_train = obj_train_approximation + value_lambda * obj_train_penalty
            delta = np.absolute(obj_train_prev - obj_train) / (
                (obj_train_prev + obj_train) / 2
            )

            # print result information
            total_time += time.time() - iter_start_time
            print(
                "Iter: {}, Total time: {}, Training Obj: {}".format(
                    iteration, total_time, obj_train
                )
            )

            if obj_train <= obj_train_prev:
                final_H = H
                final_W = W
                final_V = V

            iteration += 1

        num_rep += 1

    liger_object.W = final_W
    # Save results into the liger_object
    for i in range(num_samples):
        liger_object.adata_list[i].obsm["H"] = final_H[i].transpose()
        liger_object.adata_list[i].varm["W"] = final_W
        liger_object.adata_list[i].varm["V"] = final_V[i]

    return None
