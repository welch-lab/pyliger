# from ._qc import plot_qc


###############TODO
# Plot specific feature on t-SNE coordinates
def plotFeature(
    liger_object,
    feature,
    by_dataset=True,
    discrete=None,
    title=None,
    pt_size=0.3,
    text_size=3,
    do_shuffle=True,
    rand_seed=1,
    do_labels=False,
    axis_labels=None,
    do_legend=True,
    legend_size=5,
    option="plasma",
    cols_use=None,
    zero_color="#F5F5F5",
    return_plots=False,
):
    pass


# Plot scatter plots of unaligned and aligned factor loadings
def plotFactors(liger_object, num_genes=10, cells_highlight=None, plot_tsne=False):
    pass


# Generate word clouds and t-SNE plots
def plotWordClouds(
    liger_object,
    dataset1=None,
    dataset2=None,
    num_genes=30,
    min_size=1,
    max_size=4,
    factor_share_thresh=10,
    log_fc_thresh=1,
    umi_thresh=30,
    frac_thresh=0,
    pval_thresh=0.05,
    do_spec_plot=True,
    return_plots=False,
):
    pass


# Plot violin plots for gene expression
def plotGeneViolin(
    liger_object, gene, methylation_indices=None, by_dataset=True, return_plots=False
):
    pass


# Plot expression of multiple genes
def plotGenes(liger_object, genes):
    pass


# Generate a river (Sankey) plot
def makeRiverplot(
    liger_object,
    cluster1,
    cluster2,
    cluster_consensus=None,
    min_frac=0.05,
    min_cells=10,
    river_yscale=1,
    river_lty=0,
    river_node_margin=0.1,
    label_cex=1,
    label_col="black",
    lab_srt=0,
    river_usr=None,
    node_order="auto",
):
    pass


# ' Plot cluster proportions by dataset
def plotClusterProportions(liger_object, return_plot=False):
    pass


# Plot heatmap of cluster/factor correspondence
def plotClusterFactors(
    liger_object,
    use_aligned=False,
    Rowv=np.nan,
    Colv="Rowv",
    col=None,
    return_data=False,
):
    pass
