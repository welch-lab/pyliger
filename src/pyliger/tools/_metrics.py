import numpy as np


def calc_dataset_specificity(liger_object, dataset1=None, dataset2=None, do_plot=True):
    """Calculate a dataset-specificity score for each factor

    This score represents the relative magnitude of the dataset-specific components of each factor's
    gene loadings compared to the shared components for two datasets. First, for each dataset we
    calculate the norm of the sum of each factor's shared loadings (W) and dataset-specific loadings
    (V). We then determine the ratio of these two values and subtract from 1... TODO: finish
    description.

    Parameters
    ----------
        liger_object:
        dataset1:
        dataset2:
        do_plot:

    Returns
    -------

    """
    ### Extract values for use
    W = liger_object.W
    V1 = liger_object.V[dataset1]
    V2 = liger_object.V[dataset2]

    ### Calculation
    pct1 = np.linalg.norm((V1 + W), axis=0)
    pct2 = np.linalg.norm((V2 + W), axis=0)

    #    if do_plot:

    return [pct1, pct2, 100 * (1 - (pct1 / pct2))]


# Calculate agreement metric
def calcAgreement(
    liger_object,
    dr_method="NMF",
    ndims=40,
    k=15,
    use_aligned=True,
    rand_seed=42,
    by_dataset=False,
):
    pass


# Calculate alignment metric
def calcAlignment(
    liger_object,
    k=None,
    rand_seed=1,
    cells_use=None,
    cells_comp=None,
    clusters_use=None,
    by_cell=False,
    by_dataset=False,
):
    pass


# Calculate alignment for each cluster
def calcAlignmentPerCluster(liger_object, rand_seed=1, k=None, by_dataset=False):
    pass


# Calculate adjusted Rand index
def calcARI(liger_object, clusters_compare):
    pass


# Calculate purity
def calcPurity(liger_object, classes_compare):
    pass


# Calculate proportion mitochondrial contribution
def getProportionMito(liger_object, use_norm=False):
    pass
