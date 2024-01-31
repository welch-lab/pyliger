import warnings

import h5sparse
import matplotlib.pyplot as plt
import numexpr as ne
import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from sklearn.utils.sparsefuncs import mean_variance_axis

from pyliger._utilities import _h5_idx_generator


def select_genes(
    liger_object,
    var_thresh: float = 0.1,
    alpha_thresh: float = 0.99,
    num_genes=None,
    tol: float = 0.0001,
    datasets_use=None,
    combine="union",
    keep_unique=False,
    capitalize=False,
    do_plot=False,
    cex_use=0.3,
    chunk_size=1000,
):
    """Select a subset of informative genes

    This function identifies highly variable genes from each dataset and combines these gene sets
    (either by union or intersection) for use in downstream analysis. Assuming that gene
    expression approximately follows a Poisson distribution, this function identifies genes with
    gene expression variance above a given variance threshold (relative to mean gene expression).
    It also provides a log plot of gene variance vs gene expression (with a line indicating expected
    expression across genes and cells). Selected genes are plotted in green.

    Parameters
    ----------
    liger_object : liger object
        Should have already called normalize.
    var_thresh : float, optional
        Variance threshold. Main threshold used to identify variable genes. Genes with
        expression variance greater than threshold (relative to mean) are selected.
        (higher threshold -> fewer selected genes). Accepts single value or vector with separate
        var_thresh for each dataset (the default is 0.1).
    alpha_thresh : float, optional
        Alpha threshold. Controls upper bound for expected mean gene expression
        (lower threshold -> higher upper bound) (the default is 0.99).
    num_genes : int, optional
        Number of genes to find for each dataset. Optimises the value of var_thresh
        for each dataset to get this number of genes. Accepts single value or vector with same length
        as number of datasets (the default is None).
    tol : float, optional
        Tolerance to use for optimization if num.genes values passed in (the default is 0.0001).
    datasets_use : list, optional
        List of datasets to include for discovery of highly variable genes
        (the default is using all datatsets).
    combine : str, optional, 'union' or 'intersect'
        How to combine variable genes across experiments (the default is 'union').
    keep_unique : bool, optional
        Keep genes that occur (i.e., there is a corresponding column in raw_data) only
        in one dataset (the default is False).
    capitalize : bool, optional
        Capitalize gene names to match homologous genes (ie. across species)
        (the default is False).
    do_plot : bool, optional
        Display log plot of gene variance vs. gene expression for each dataset.
        Selected genes are plotted in green (the default is False).
    cex_use : float, optional
        Point size for plot (the default is 0.3).

    Returns
    -------
    liger_object : liger object
        Object with var_genes attribute.

    Examples
    --------
    >>> adata1 = AnnData(np.arange(12).reshape((4, 3)))
    >>> adata2 = AnnData(np.arange(12).reshape((4, 3)))
    >>> ligerex = pyliger.create_liger([adata1, adata2])
    >>> pyliger.normalize(ligerex)
    >>> pyliger.select_genes(ligerex) # use default selectGenes settings
    >>> pyliger.select_genes(ligerex, var_thresh=0.8) # select a smaller subset of genes
    """
    num_samples = liger_object.num_samples

    if datasets_use is None:
        datasets_use = list(range(num_samples))

    # Expand if only single var_thresh passed
    if isinstance(var_thresh, int) or isinstance(var_thresh, float):
        var_thresh = np.repeat(var_thresh, num_samples)
    if num_genes is not None:
        num_genes = np.repeat(num_genes, num_samples)

    if not np.array_equal(
        np.intersect1d(datasets_use, list(range(num_samples))), datasets_use
    ):
        datasets_use = np.intersect1d(datasets_use, list(range(num_samples)))

    genes_use = np.array([])
    # Iterate through datasets
    for i in datasets_use:
        adata = liger_object.adata_list[i]
        if capitalize:
            liger_object.adata_list[i].var.index = liger_object.adata_list[
                i
            ].var.index.str.upper()

        trx_per_cell = liger_object.adata_list[i].obs["nUMI"].to_numpy()
        gene_expr_mean = liger_object.adata_list[i].var["norm_mean"].to_numpy()

        # On-disk mode (set for online learning approach)
        if adata.isbacked:
            gene_expr_var = _calc_var_online(adata, gene_expr_mean, chunk_size)

        # In-memory mode
        else:
            # each gene's expression variance (across all cells)
            gene_expr_var = _calc_var_matrix(adata)

        nolan_constant = np.mean(ne.evaluate("1 / trx_per_cell"))
        alpha_thresh_corrected = alpha_thresh / adata.shape[1]
        gene_mean_upper = gene_expr_mean + norm.ppf(
            1 - alpha_thresh_corrected / 2
        ) * np.sqrt(ne.evaluate("gene_expr_mean * nolan_constant") / adata.shape[0])
        # base_gene_lower = np.log10(ne.evaluate('gene_expr_mean * nolan_constant'))
        base_gene_lower = ne.evaluate("log10(gene_expr_mean * nolan_constant)")

        def num_var_genes(x, num_genes_des):
            # This function returns the difference between the desired number of genes and
            # the number actually obtained when thresholded on x
            y = np.sum(
                (gene_expr_var / nolan_constant)
                > gene_mean_upper & np.log10(gene_expr_var, dtype=np.float64)
                > (base_gene_lower + x)
            )
            return np.abs(num_genes_des - y)

        if num_genes is not None:
            # Optimize to find value of x which gives the desired number of genes for this dataset
            # if very small number of genes requested, var.thresh may need to exceed 1
            optimized = minimize(
                fun=num_var_genes, x0=[0], agrs=num_genes[i], tol=tol, bounds=[(0, 1.5)]
            )
            var_thresh[i] = optimized.x
            if var_thresh[i].shape[0] > 1:
                warnings.warn(
                    "Returned number of genes for dataset {} differs from requested by {}. Lower tol or alpha_thresh for better results.".format(
                        i, optimized.x.shape[0]
                    )
                )

        temp = var_thresh[i]
        select_gene = ne.evaluate(
            "((gene_expr_var / nolan_constant) > gene_mean_upper) & (log10(gene_expr_var) > (base_gene_lower + temp))"
        )
        # select_gene = ((gene_expr_var / nolan_constant) > gene_mean_upper) & (
        #            np.log10(gene_expr_var) > (base_gene_lower + var_thresh[i]))
        genes_new = adata.var.index.to_numpy()[select_gene]

        # TODO: graph needs to be improved
        if do_plot:
            plt.plot(np.log10(gene_expr_mean), np.log10(gene_expr_var))
            plt.title(liger_object.adata_list[i].uns["sample_name"])
            plt.xlabel("Gene Expression Mean (log10)")
            plt.ylabel("Gene Expression Variance (log10)")

            plt.plot(
                np.log10(gene_expr_mean),
                np.log10(gene_expr_mean) + nolan_constant,
                c="p",
            )
            plt.scatter(
                np.log10(gene_expr_mean[select_gene]),
                np.log10(gene_expr_var[select_gene]),
                c="g",
            )
            plt.legend("Selected genes: " + str(len(genes_new)), loc="best")

            plt.show()

        if combine == "union":
            genes_use = np.union1d(genes_use, genes_new)

        if combine == "intersect":
            if genes_use.shape[0] == 0:
                genes_use = genes_new
            genes_use = np.intersect1d(genes_use, genes_new)

    if not keep_unique:
        for i in range(num_samples):
            genes_use = np.intersect1d(genes_use, liger_object.adata_list[i].var.index)

    if genes_use.shape[0] == 0:
        warnings.warn(
            'No genes were selected; lower var_thresh values or choose "union" for combine parameter'
        )

    liger_object.var_genes = genes_use
    for idx, adata in enumerate(liger_object.adata_list):
        var_gene_idx = (
            liger_object.adata_list[idx].var.index.isin(genes_use).nonzero()[0]
        )
        liger_object.adata_list[idx].uns["var_gene_idx"] = var_gene_idx
    return None


def _calc_var_online(adata, gene_expr_mean, chunk_size):
    """"""
    sum_sq_dev = np.zeros(adata.shape[1])

    file_path = "./results/" + adata.uns["sample_name"] + ".hdf5"
    with h5sparse.File(file_path, "r") as f:
        for left, right in _h5_idx_generator(chunk_size, adata.shape[0]):
            sum_sq_dev += np.sum(
                np.square(
                    np.absolute(
                        np.subtract(
                            f["norm_data"][left:right].toarray(), gene_expr_mean
                        )
                    )
                ),
                axis=0,
            )
    gene_expr_var = sum_sq_dev / (adata.shape[0] - 1)
    adata.var["norm_var"] = gene_expr_var

    return gene_expr_var


def _calc_var_matrix(adata):
    """"""
    # gene_expr_var = np.ravel(
    #    np.var(adata.layers['norm_data'].toarray(), axis=0, dtype=np.float64, ddof=1))

    num_sampels = adata.shape[0]
    _, gene_expr_var = mean_variance_axis(adata.layers["norm_data"], axis=0)
    gene_expr_var = ne.evaluate(
        "gene_expr_var * num_sampels / (num_sampels-1)"
    )  # change ddof
    adata.var["norm_var"] = gene_expr_var

    return gene_expr_var
