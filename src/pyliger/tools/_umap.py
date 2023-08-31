import numpy as np
import pandas as pd
import umap
import umap.plot


def run_umap(
    liger_object,
    use_raw=False,
    dims_use=None,
    k=2,
    distance="euclidean",
    n_neighbors=10,
    min_dist=0.1,
    rand_seed=42,
):
    """Perform UMAP dimensionality reduction

    Run UMAP on the normalized cell factors (or raw cell factors) to generate a 2D embedding for
    visualization (or general dimensionality reduction). Has option to run on subset of factors.

    Note that running multiple times will overwrite tsne_coords values. It is generally
    recommended to use this method for dimensionality reduction with extremely large datasets.

    Parameters
    ----------
    liger_object : liger
        Should run quantile_norm before calling with defaults.
    use_raw : bool, optional
        Whether to use un-aligned cell factor loadings (H matrices) (the default is False).
    dims_use : list
        Factors to use for computing tSNE embedding (the default is using all cells).
    k : int, optional
        Number of dimensions to reduce to (the default is 2).
    distance : str, optional
        Mtric used to measure distance in the input space. A wide variety of metrics are
        already coded, and a user defined function can be passed as long as it has been JITd by numba.
        (the default is "euclidean", alternatives: "cosine", "manhattan", "hamming").
    n_neighbors : int, optional
        Number of neighboring points used in local approximations of manifold
        structure. Larger values will result in more global structure being preserved at the loss of
        detailed local structure. In general this parameter should often be in the range 5 to 50, with
        a choice of 10 to 15 being a sensible default (the default is 10).
    min_dist : float, optional
        Controls how tightly the embedding is allowed compress points together. Larger
        values ensure embedded points are more evenly distributed, while smaller values allow the
        algorithm to optimise more accurately with regard to local structure. Sensible values are in
        the range 0.001 to 0.5, with 0.1 being a reasonable default (the default is 0.1).
    rand_seed : int, optional
        Random seed for reproducibility (the default is 42).

    Returns
    -------
    liger_object : liger
        object with tsne_coords attribute.

    Examples
    --------
    >>> ligerex = quantile_norm(ligerex) # generate H_norm by quantile aligning factor loadings
    >>> ligerex = run_umap(ligerex) # get tsne_coords for normalized data
    >>> ligerex = run_umap(ligerex, use_raw = T) # get tsne.coords for raw factor loadings
    """
    np.random.seed(rand_seed)

    ### Run umap
    if use_raw:
        raw_data = np.vstack([adata.obsm["H"] for adata in liger_object.adata_list])

        # if H_norm not set yet
        if "H_norm" not in liger_object.adata_list[0].obs_keys():
            dims_use = list(range(raw_data.shape[1]))
        else:
            H_norm = np.vstack(
                [adata.obsm["H_norm"] for adata in liger_object.adata_list]
            )
            dims_use = list(range(H_norm.shape[1]))

        umap_coords = umap.UMAP(
            n_components=k,
            metric=distance,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=rand_seed,
        ).fit_transform(raw_data[:, dims_use])
    else:
        # if dims_use is None:
        #    H_norm = np.vstack([adata.obsm['H_norm'] for adata in liger_object.adata_list])
        #    dims_use = list(range(H_norm.shape[1]))

        # tsne_coords = umap.UMAP(n_components=k, metric=distance,
        #                        n_neighbors=n_neighbors, min_dist=min_dist,
        #                        random_state=rand_seed).fit(H_norm[dims_use, :])
        H_norm = np.vstack([adata.obsm["H_norm"] for adata in liger_object.adata_list])
        umap_coords = umap.UMAP(
            n_components=k,
            metric=distance,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=rand_seed,
        ).fit_transform(H_norm)

    ### Save umap results
    liger_object.save_obsm(umap_coords, "umap_coords")

    liger_object.tsne_coords = pd.DataFrame(umap_coords, columns=["tsne1", "tsne2"])
    return None


# Perform t-SNE dimensionality reduction
def runTSNE(
    liger_object,
    dims_use,
    use_raw=False,
    use_pca=False,
    perplexity=30,
    theta=0.5,
    method="Rtsne",
    fitsne_path=None,
    rand_seed=42,
):
    # dims_use = range(1, len(liger_object.H_norm))
    pass
