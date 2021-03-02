import h5sparse
import numpy as np

from ..preprocessing import *
from ._utilities import nonneg, nnlsm_blockpivot, _init_W, _init_V, _update_W_HALS, _update_V_HALS
from .._utilities import _h5_idx_generator

def online_iNMF(liger_object,
                X_new = None,
                projection = False,
                W_init = None,
                V_init = None,
                H_init = None,
                A_init = None,
                B_init = None,
                k = 20,
                value_lambda = 5.0,
                max_epochs = 5,
                miniBatch_max_iters = 1,
                miniBatch_size = 5000,
                h5_chunk_size = 1000,
                rand_seed = 1):
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
        miniBatch_max_iters:
        miniBatch_size:
        h5_chunk_size:
        rand_seed:

    Returns
    -------

    """
    pass
    """
    # if there is new dataset
    if X_new:
        # check whether X_new needs to be processed
        #for adata in X_new:
            #if 'scale_data' not in adata.layers.keys():

        # assuming only one new dataset arrives at a time
        #liger_object.adata_list = liger_object.adata_list + X_new
        #TODO: deal with new datasets
        pass

    ### Extract required information and initialize algorithm
    num_files = liger_object.num_samples # number of total input hdf5 files
    num_new_files = num_files if X_new else 1  # number of new input hdf5 files since last step
    num_prev_files = 0 if X_new else num_files - num_new_files  # number of input hdf5 files processed in last step

    file_idx = set(range(num_files))
    file_idx_new = set(range(num_prev_files+1, num_files))
    file_idx_prev = file_idx-file_idx_new

    file_names = liger_object.sample_names
    num_genes = liger_object.num_var_genes # number of the selected genes
    gene_names = liger_object.var_genes # genes selected for analysis

    cell_barcodes = [adata.obs.index.to_numpy() for adata in liger_object.adata_list] # cell barcodes for each dataset
    num_cells = [adata.shape[0] for adata in liger_object.adata_list] # number of cells in each dataset
    num_cells_new = num_cells[num_prev_files:num_files]

    miniBatch_sizes = np.round((num_cells[i]/np.sum(num_cells[file_idx_new])) * miniBatch_size for i in file_idx_new:)
    miniBatch_sizes_orig = miniBatch_sizes

    Xs = []
    for adata in liger_object.adata_list:
        file_path = './results/' + adata.uns['sample_name'] + '.hdf5'
        Xs.append(h5sparse.File(file_path,'r'))

    if not projection:
        # set seed
        np.random.seed(seed=rand_seed)

        # Initialization
        if X_new is None:
            W = _init_W(num_genes, k, rand_seed)
            V = _init_V(num_cells, num_files, k, Xs)
        else:
            W = W_init if W_init else liger_object.W
            V = V_init if V_init else liger_object.V

        # H_i matrices initialization
        if X_new is None:
            H = []
            H_miniBatch = []
        else:
            H_miniBatch = []

        # A = HiHi^t, B = XiHit
        A_old = []
        B_old = []

        if X_new is None:
            A = 0
            B = 0
            A_old = 0
            B_old = 0
        else:
            pass

        iter = 1
        epoch = np.repeat(0, num_files) # intialize the number of epoch for each dataset
        epoch_prev = np.repeat(0, num_files) # intialize the previous number of epoch for each dataset
        epoch_next = np.repeat(False, num_files)
        sqrt_lambda = np.sqrt(value_lambda)
        total_time = 0 # track the total amount of time used for the online learning

        # chunk permutation
        num_chunks, chunk_idx = _chunk_permutation(num_files, num_cells, file_idx_new, h5_chunk_size)

        # print start information
        print('Starting Online iNMF...')

        total_iters = np.floor(np.sum(num_cells_new) * max_epochs / miniBatch_size)
        while epoch[file_idx_new[0]] < max_epochs:
            # track epochs

            miniBatch_idx = _track_epoch(max_epochs, num_cells_new, epoch, epoch_next, epoch_prev)
            miniBatch_idx = np.zeros(np.nan, num_files) # indices of samples in each dataest used for this iteration
            # indices of samples in each dataest used for this iteration
            if (max_epochs * num_cells_new[0] - (iter-1) * miniBatch_sizes[file_idx_new[0]]) >= miniBatch_sizes[file_idx_new[0]]:
                for i in file_idx_new:
                    epoch[i] = (iter * miniBatch_sizes[i]) // num_cells[i] # caculate the current epoch
                    # if current iter cycles through the data and start a new cycle
                    if epoch_prev[i] != epoch[i] and ((iter * miniBatch_sizes[i]) % num_cells[i]) != 0:
                        epoch_next[i] = True
                        epoch_prev[i] = epoch[i]

                        # shuffle dataset before the next epoch
                        _ = _dataset_shuffle()

                        all_idx[i] = all_idx[i].pop(0) # remove the first element 0
                        miniBatch_idx[i] =
                    # if current iter finishes this cycle without start a a new cycle
                    elif epoch_prev[i]:
                        epoch_next[i] = True
                    # if current iter stays within a single cycle
                    else:
                        miniBatch_idx[i]
            else:
                for i in file_idx_new:
                    miniBatch_sizes[i]
                    miniBatch_idx[i]
                epoch[file_idx_new[0]] = max_epochs # last epoch

            if len(miniBatch_idx[file_idx_new[0]] == miniBatch_sizes_orig[file_idx_new[0]]):
                X_miniBatch = [Xs[miniBatch_idx[]]]

                # update H_i by ANLS Hi_minibatch[[i]]
                for i in file_idx_new:
                    H_miniBatch = [nnlsm_blockpivot(A=np.hstack(((W + V[i]), sqrt_lambda * V[i])),
                                B=np.hstack((X_miniBatch[i], np.zeros((miniBatch_sizes[i], num_genes)))))[
                        0].transpose()]

                # updata A and B matrices
                if iter == 1:
                    scale_


                # update W, V_i by HALS
                iter_miniBatch = 1
                delta_miniBatch = np.inf
                max_iters_miniBatch = miniBatch_max_iters

                while iter_miniBatch <= max_iters_miniBatch:
                    # update W matrix
                    W = _update_W_HALS(A, B, W, V)
                    #for j in range(k):
                    #    W_update_numerator = np.repeat(0, num_genes)
                    #    W_update_denominator = 0
                    #    for i in file_idx:
                    #        W_update_numerator = W_update_numerator + B[i][:, j] - ((W + V[i]) @ A[i])[:, j]
                    #        W_update_denominator += A[i][j, j]
                    #    W[:, j] = nonneg(W[:, j] + W_update_numerator / W_update_denominator)

                    # update V matrix
                    V[file_idx_new] = _update_V_HALS(A, B, W, V[file_idx_new], value_lambda)
                    #for j in range(k):
                    #    for i in file_idx_new:
                    #        V[i][:, j] = nonneg(
                    #            V[i][:, j] + (B[i][:, j] - (W + (1 + value_lambda) * V[i]) @ A[i][:, j])
                    #            / ((1 + value_lambda) * A[i][j, j]))

                    iter_miniBatch += 1

                epoch_next = np.repeat(False, num_files) # reset epoch change indicator
                iter += 1

        print('Calculate metagene loadings...')
        H

    else:
        print('Metagene projection')

    # Save results into the liger_object
    for i in range(num_files):
        liger_object.adata_list[i].obsm['H'] = H[i]
        liger_object.adata_list[i].varm['W'] = W.transpose()
        liger_object.adata_list[i].varm['V'] = V[i].transpose()

    # Close all files in the end
    for file in Xs:
        file.close()

    return liger_object


def _online_iNMF_from_scratch():


    return

def _online_iNMF_refine():
    return

def _online_iNMF(W, H, V, ):

    return None


def _update_cell_loading():
    return None

def _chunk_permutation(num_files, num_cells, file_idx_new, h5_chunk_size):
    num_chunks = [np.ceil(num_cell/h5_chunk_size) for num_cell in num_cells]
    chunk_idx = [np.random.permutation(num_chunk) for num_chunk in num_chunks]
    #all_idx = [list(range(left, right)) for left, right in _h5_idx_generator(h5_chunk_size, num_cell) for num_cell in num_cells]

    return num_chunks, chunk_idx


def _track_epoch():

    return


def _dataset_shuffle():
    miniBatch_idx[i] = all_idx[i][]
    chunk_idx[i] =
    all_idx[i] = 0
    for j in chunk_idx[i]:
        if j != num_chunks[i]:
            all_idx[] =
        else:
            all_idx[] =


def _projection():

    return None
"""