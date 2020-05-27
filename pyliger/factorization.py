import time
import numpy as np

from .utilities import nonneg

#######################################################################################
#### Factorization

def optimizeALS(liger_object,
                k,
                value_lambda = 5.0,
                thresh = 1e-6,
                max_iters = 30,
                nrep = 1,
                H_init = None,
                W_init = None,
                V_init = None,
                rand_seed = 1,
                print_obj = False):
    """ Perform iNMF on scaled datasets
    
    Perform integrative non-negative matrix factorization to return factorized H, W, and V matrices.
    It optimizes the iNMF objective function using block coordinate descent (alternating non-negative
    least squares), where the number of factors is set by k. TODO: include objective function
    equation here in documentation (using deqn)
    
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
        >>> ligerex = createLiger([adata1, adata2])
        >>> ligerex = normalize(ligerex)
        >>> ligerex = selectGenes(ligerex) # select genes
        >>> ligerex = scaleNotCenter(ligerex)
        >>> ligerex = optimizeALS(ligerex, k = 20, lambda = 5, nrep = 3) # get factorization using three restarts and 20 factors
    """
    """
    for :
        if :
            raise ValueError('All values in "object" must be a matrix')
    
    E = liger_object
    N = 
    ns = 
    if k >= :
        raise ValueError('Select k lower than the number of cells in smallest dataset: {}'. format())
        
    tmp = 
    g = 
    if k >= g:
        raise ValueError('Select k lower than the number of variable genes: {}'. format(g))
        
    W_m = csr_matrix()
    V_m = []
    for i in range(N):
        V_m.append(csr_matrix())
        
    tmp = 
    best_obj = np.Inf
    run_stats = csr_matrix()
    for i in range(1, nrep+1):
        np.random.seed(seed = rand_seed + i - 1)
        start_time = time.time()
        W = 
        V = []
        for i in range(N):
            V.append()
        
        H = []
        for i in ns:
            H.append()
        
        tmp = 
        if W_init is not None:
            W = W_init
            
        if V_init is not None:
            V = V_init
        
        if H_init is not None:
            H = H_init
            
        delta = 1
        iters = 0
        pb = 
        sqrt_lambda = np.sqrt(value_lambda)
        obj0 = 
        
        tmp = 
        
        while delta > thresh and iters < max_iters:
            for i in range():
                H =
            
            tmp = 
            for i in range():
                V = 
                
            tmp = 
            W = 
            
            tmp =
            obj = 
            
            
        if print_obj:
            print('Objective: {}'.format())
        
    print('Best results with seed {}.'.format(best_seed))
    
    liger_object.H = H_m
    liger_object.V = V_m
    liger_object.W = W_m
    
    return liger_object
    """
    pass

def iNMF_HALS(liger_object,
              k = 20,
              value_lambda = 5.0,
              thresh = 1e-4,
              max_iters = 25,
              nrep = 1,
              rand_seed = 1): 
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
    >>> ligerex = createLiger([adata1, adata2])
    >>> ligerex = normalize(ligerex)
    >>> ligerex = selectGenes(ligerex) # select genes
    >>> ligerex = scaleNotCenter(ligerex)
    >>> ligerex = iNMF_HALS(ligerex, k = 20, lambda = 5, nrep = 3)
    """

    adata_list = liger_object.adata_list
    num_samples = len(adata_list)
    num_cells = [adata_list[i].shape[1] for i in range(num_samples)]
    num_genes = adata_list[0].shape[0]
    
    # set seed
    np.random.seed(seed = rand_seed)

    # allow restarts
    num_rep = 0
    while num_rep < nrep:
        # initialization
        V = [adata_list[i].layers['scale_data'][:,np.random.choice(list(range(num_cells[i])), k)].toarray() for i in range(num_samples)]
        W = np.abs(np.random.uniform(0, 2, (num_genes, k)))
        H = [np.random.uniform(0, 2, (k, num_cells[i])) for i in range(num_samples)]
        
        # normalize columns of dictionaries
        W = W/np.sqrt(np.sum(np.square(W), axis=0))
        V = [V[i]/np.sqrt(np.sum(np.square(V[i]), axis=0)) for i in range(num_samples)]
        
        # initial training obj
        obj_train_approximation = 0
        obj_train_penalty = 0
        for i in range(num_samples):
            obj_train_approximation += np.linalg.norm(np.float64(adata_list[i].layers['scale_data'].toarray())-np.float64((W+V[i]))@np.float64((H[i])))**2
            obj_train_penalty += np.linalg.norm(V[i]@H[i])**2

        obj_train = obj_train_approximation + value_lambda*obj_train_penalty

        # print start information
        print('Initial Training Obj: {}'.format(obj_train))
            
        iteration = 1
        total_time = 0 # track the total amount of time used for learning
        delta = np.inf
        
        # perform block coordinate descent
        while delta > thresh and iteration <= max_iters:
            
            iter_start_time = time.time()

            VitVi = [V[i].transpose()@V[i] for i in range(num_samples)]
            W_Vi = [W+V[i] for i in range(num_samples)]
            W_Vi_sq = [W_Vi[i].transpose()@W_Vi[i] for i in range(num_samples)]
            W_VitXi = [W_Vi[i].transpose()@adata_list[i].layers['scale_data'] for i in range(num_samples)]

            for i in range(num_samples):
                for j in range(k):
                    H[i][j,:] = nonneg((H[i][j,:]*W_Vi_sq[i][j,j] + W_VitXi[i][j,:] - (W_Vi_sq[i]@H[i])[j,:])/(W_Vi_sq[i][j,j] + value_lambda * VitVi[i][j,j]))
                    
            XHt = [adata_list[i].layers['scale_data']@H[i].transpose() for i in range(num_samples)]
            HHt = [H[i]@H[i].transpose() for i in range(num_samples)]

            for j in range(k):
                W_update_numerator = np.repeat(0, num_genes)
                W_update_denominator = 0
                for i in range(num_samples):
                    W_update_numerator[:] = W_update_numerator[:] + XHt[i][:,j] - ((W + V[i])@HHt[i])[:,j]
                    W_update_denominator += HHt[i][j,j]
                W[:,j] = nonneg(W[:,j] + W_update_numerator / W_update_denominator)
                
            # update (Di).k
            for j in range(k):
                for i in range(num_samples):
                    V[i][:,j] = nonneg(V[i][:,j] / (1 + value_lambda) + (XHt[i][:,j]-((W+V[i])@HHt[i])[:,j])/((1+value_lambda)*HHt[i][j,j]))
                    
            total_time += (time.time()-iter_start_time)
            
            # training obj
            obj_train_prev = obj_train
            obj_train_approximation = 0
            obj_train_penalty = 0
            for i in range(num_samples):
                obj_train_approximation += np.linalg.norm(adata_list[i].layers['scale_data'] - (W + V[i])@H[i])**2
                obj_train_penalty += np.linalg.norm(V[i]@H[i])**2
            obj_train = obj_train_approximation + value_lambda*obj_train_penalty
            delta = np.absolute(obj_train_prev - obj_train)/((obj_train_prev + obj_train)/2)
            
            # result
            print('Iter: {}, Total time: {}, Obj Delta: {}'.format(iteration, total_time, delta))
            
            iteration += 1
            
        num_rep += 1
    
    # Save results to the liger_object
    for i in range(num_samples):
        liger_object.adata_list[i].varm['H'] = H[i].transpose()
        liger_object.adata_list[i].obsm['W'] = W
        liger_object.adata_list[i].obsm['V'] = V[i]
        
    return liger_object

# Perform factorization for new value of k
def optimizeNewK(liger_object, k_new, value_lambda = None, thresh = 1e-4, max_iters = 100,
                 rand_seed = 1):
    pass

# Perform factorization for new data
def optimizeNewData(liger_object, new_data, which_datasets, add_to_existing = True,
                     value_lambda = None, thresh = 1e-4, max_iters = 100):
    pass

# Perform factorization for subset of data
def optimizeSubset(liger_object, cell_subset = None, cluster_subset = None, 
                    value_lambda = None, thresh = 1e-4, max_iters = 100, datasets_scale = None):
    pass

# Perform factorization for new lambda value
def optimizeNewLambda(liger_object, new_lambda, thresh = 1e-4, max_iters = 100, rand_seed = 1):
    pass

# Visually suggest appropriate lambda value
def suggestLambda(liger_object, k, lambda_test = None, rand_seed = 1, num_cores = 1,
                  thresh = 1e-4, max_iters = 100, knn_k = 20, k2 = 500, ref_dataset = None,
                  resolution = 1, gen_new = False, nrep = 1, return_data = False, return_raw = False):
    pass

# Visually suggest appropiate k value
# k_test range is set from 5 to 55
def suggestK(liger_object, k_test = range(5, 55, 5), value_lambda = 5, thresh = 1e-4, max_iters = 100,
             num_cores = 1, rand_seed = 1, gen_new = False, nrep = 1, plot_log2 = True,
             return_data = False, return_raw = False):
    
    pass
