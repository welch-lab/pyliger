import umap
import umap.plot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#######################################################################################
#### Dimensionality Reduction
    
# Perform t-SNE dimensionality reduction
def runTSNE(liger_object,
            dims_use, 
            use_raw = False, 
            use_pca = False, 
            perplexity = 30,
            theta = 0.5, 
            method = "Rtsne", 
            fitsne_path = None, 
            rand_seed = 42):
    dims_use = range(1, len(liger_object.H_norm))
    pass


def run_umap(liger_object,
            use_raw = False, 
            dims_use = None, 
            k = 2, 
            distance = "euclidean",
            n_neighbors = 10,
            min_dist = 0.1, 
            rand_seed = 42):
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
        Factors to use for computing tSNE embedding (the default is 1:ncol(H.norm)).
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
    ligerex <- quantile_norm(ligerex) # generate H_norm by quantile aligning factor loadings
    ligerex <- run_umap(ligerex) # get tsne_coords for normalized data
    ligerex <- run_umap(ligerex, use_raw = T) # get tsne.coords for raw factor loadings
    """

    raw_data = None
    if use_raw:
        for i in range(len(liger_object.adata_list)):
            if raw_data is None:
                raw_data = liger_object.adata_list[i].varm['H']
            else:
                raw_data = np.vstack((raw_data, liger_object.adata_list[i].varm['H']))
        
        # if H_norm not set yet
        if not hasattr(liger_object, 'H_norm'):
            dims_use = list(range(raw_data.shape[1]))
        else:
            dims_use = list(range(liger_object.H_norm.shape[1]))
            
        tsne_coords = umap.UMAP(n_components=k, metric=distance,
                                n_neighbors=n_neighbors, min_dist=min_dist, 
                                random_state=rand_seed).fit_transform(raw_data[:,dims_use])
    else:
        if dims_use is None:
            dims_use = list(range(liger_object.H_norm.shape[1]))
        
        H_norm = liger_object.H_norm.to_numpy()  
        
        tsne_coords = umap.UMAP(n_components=k, metric=distance,
                                n_neighbors=n_neighbors, min_dist=min_dist, 
                                random_state=rand_seed).fit(H_norm[:,dims_use])
    
    #liger_object.tsne_coords = pd.DataFrame(tsne_coords, index=liger_object.H_norm.index,
    #                                        columns = ['tsne1', 'tsne2'])
    liger_object.tsne_coords = tsne_coords
    return liger_object

#######################################################################################
#### Visualization

def plot_by_dataset_and_cluster(liger_object, 
                                clusters = None, 
                                title = None, 
                                pt_size = 0.3,
                                text_size = 3, 
                                do_shuffle = True, 
                                rand_seed = 1,
                                axis_labels = None, 
                                do_legend = True, 
                                legend_size = 5,
                                return_plots = False):
    """Plot t-SNE coordinates of cells across datasets
    
    Generates two plots of all cells across datasets, one colored by dataset and one colored by
    cluster. These are useful for visually examining the alignment and cluster distributions,
    respectively. If clusters have not been set yet (quantileAlignSNF not called), will plot by
    single color for second plot. It is also possible to pass in another clustering (as long as
    names match those of cells).

    Parameters
    ----------
    liger_object : TYPE
        DESCRIPTION.
    clusters : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is None.
    pt_size : TYPE, optional
        DESCRIPTION. The default is 0.3.
    text_size : TYPE, optional
        DESCRIPTION. The default is 3.
    do_shuffle : TYPE, optional
        DESCRIPTION. The default is True.
    rand_seed : TYPE, optional
        DESCRIPTION. The default is 1.
    axis_labels : TYPE, optional
        DESCRIPTION. The default is None.
    do_legend : TYPE, optional
        DESCRIPTION. The default is True.
    legend_size : TYPE, optional
        DESCRIPTION. The default is 5.
    return_plots : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    
    
    return None
    

# Plot specific feature on t-SNE coordinates
def plotFeature(liger_object, feature, by_dataset = True, discrete = None, title = None, 
                pt_size = 0.3, text_size = 3, do_shuffle = True, rand_seed = 1, do_labels = False,
                axis_labels = None, do_legend = True, legend_size = 5, option = 'plasma', 
                cols_use = None, zero_color = '#F5F5F5', return_plots = False):
    pass

# Plot scatter plots of unaligned and aligned factor loadings
def plotFactors(liger_object, num_genes = 10, cells_highlight = None, plot_tsne = False):
    pass

# Generate word clouds and t-SNE plots
def plotWordClouds(liger_object, dataset1 = None, dataset2 = None, num_genes = 30, min_size = 1,
                   max_size = 4, factor_share_thresh = 10, log_fc_thresh = 1,
                   umi_thresh = 30, frac_thresh = 0, pval_thresh = 0.05,
                   do_spec_plot = True, return_plots = False):
    pass

# Generate t-SNE plots and gene loading plots
def plotGeneLoadings(liger_object, dataset1 = None, dataset2 = None, num_genes_show = 12,
                     num_genes = 30, mark_top_genes = True, factor_share_thresh = 10,
                     log_fc_thresh = 1, umi_thresh = 30, frac_thresh = 0,
                     pval_thresh = 0.05, do_spec_plot = True, max_val = 0.1, pt_size = 0.1,
                     option = "plasma", zero_color = "#F5F5F5", return_plots = False,
                     axis_labels = None, do_title = False):
    pass

# Plot violin plots for gene expression
def plotGeneViolin(liger_object, gene, methylation_indices = None,
                   by_dataset = True, return_plots = False):
    pass

# Plot gene expression on dimensional reduction (t-SNE) coordinates
def plotGene(liger_object, gene, use_raw = False, use_scaled = False, scale_by = 'dataset', 
             log2scale = None, methylation_indices = None, plot_by = 'dataset', 
             set_dr_lims = False, pt_size = 0.1, min_clip = None, max_clip = None, 
             clip_absolute = False, points_only = False, option = 'plasma', cols_use = None, 
             zero_color = '#F5F5F5', axis_labels = None, do_legend = True, return_plots = False):
    pass

# Plot expression of multiple genes
def plotGenes(liger_object, genes):
    pass

# Generate a river (Sankey) plot
def makeRiverplot(liger_object, cluster1, cluster2, cluster_consensus = None, min_frac = 0.05,
                  min_cells = 10, river_yscale = 1, river_lty = 0, river_node_margin = 0.1,
                  label_cex = 1, label_col = "black", lab_srt = 0, river_usr = None,
                  node_order = "auto"):
    pass

#' Plot cluster proportions by dataset
def plotClusterProportions(liger_object, return_plot = False):
    pass

# Plot heatmap of cluster/factor correspondence
def plotClusterFactors(liger_object, use_aligned = False, Rowv = np.nan, Colv = "Rowv", col = None,
                       return_data = False):
    pass
