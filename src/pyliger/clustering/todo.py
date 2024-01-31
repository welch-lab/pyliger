#####TODO
def imputeKNN(
    liger_object, reference, queries, knn_k=20, weight=True, norm=True, scale=False
):
    """Impute the query cell expression matrix

    Impute query features from a reference dataset using KNN.

    Parameters
    ----------
    liger_object : liger obect
        object.
    reference : str
        Dataset containing values to impute into query dataset(s).
    queries : TYPE
        Dataset to be augmented by imputation. If not specified, will pass in all datasets.
    knn_k : int, optional
        The maximum number of nearest neighbors to search (the default is 20).
    weight : bool, optional
        Whether to use KNN distances as weight matrix (the default is False).
    norm : bool, optional
        Whether normalize the imputed data with default parameters (the default is True).
    scale : bool, optional
        Whether scale but not center the imputed data with default parameters (the default is True).

    Returns
    -------
    liger_object : liger obect
        object with raw data in raw.data slot replaced by imputed data (genes by cells)
        # TODO: may not be replaced?

    Examples
    --------
    >>>
    """
    pass


# Linking genes to putative regulatory elements
def linkGenesAndPeaks(
    gene_counts,
    peak_counts,
    path_to_coords,
    genes_list=None,
    dist="spearman",
    alpha=0.05,
):
    pass


# Export predicted gene-pair interaction
def makeInteractTrack(corr_mat, genes_list, output_path, path_to_coords):
    pass


# Analyze biological interpretations of metagene
def runGSEA(liger_object, gene_sets=[], mat_w=True, mat_v=0, custom_gene_sets=[]):
    pass
