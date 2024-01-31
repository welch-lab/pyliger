
#######################################################################################
#### Factorization TODO


# Perform factorization for new value of k
def optimizeNewK(
    liger_object, k_new, value_lambda=None, thresh=1e-4, max_iters=100, rand_seed=1
):
    pass


# Perform factorization for new data
def optimizeNewData(
    liger_object,
    new_data,
    which_datasets,
    add_to_existing=True,
    value_lambda=None,
    thresh=1e-4,
    max_iters=100,
):
    pass


# Perform factorization for subset of data
def optimizeSubset(
    liger_object,
    cell_subset=None,
    cluster_subset=None,
    value_lambda=None,
    thresh=1e-4,
    max_iters=100,
    datasets_scale=None,
):
    pass


# Perform factorization for new lambda value
def optimizeNewLambda(
    liger_object, new_lambda, thresh=1e-4, max_iters=100, rand_seed=1
):
    pass


# Visually suggest appropriate lambda value
def suggestLambda(
    liger_object,
    k,
    lambda_test=None,
    rand_seed=1,
    num_cores=1,
    thresh=1e-4,
    max_iters=100,
    knn_k=20,
    k2=500,
    ref_dataset=None,
    resolution=1,
    gen_new=False,
    nrep=1,
    return_data=False,
    return_raw=False,
):
    pass


# Visually suggest appropiate k value
# k_test range is set from 5 to 55
def suggestK(
    liger_object,
    k_test=range(5, 55, 5),
    value_lambda=5,
    thresh=1e-4,
    max_iters=100,
    num_cores=1,
    rand_seed=1,
    gen_new=False,
    nrep=1,
    plot_log2=True,
    return_data=False,
    return_raw=False,
):
    pass
