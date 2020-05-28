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

# Perform UMAP dimensionality reduction
def runUMAP(liger_object,
            dims_use, 
            use_raw = False, 
            k = 2, 
            distance = "euclidean",
            n_neighbors = 10,
            min_dist = 0.1, 
            rand_seed = 42):
    dims_use = range(1, len(liger_object.H_norm))
    pass

#######################################################################################
#### Visualization

# Plot t-SNE coordinates of cells across datasets
def plotByDatasetAndCluster(liger_object, clusters = None, title = None, pt_size = 0.3,
                            text_size = 3, do_shuffle = True, rand_seed = 1,
                            axis_labels = None, do_legend = True, legend_size = 5,
                            return_plots = False):
    pass

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
