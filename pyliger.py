# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

## The LIGER Class
class Liger(object):
    """ Main LIGER class
    
    The liger object is created from two or more single cell datasets. To construct a
    liger object, the user needs to provide at least two expression (or another
    single-cell modality) matrices. The class provides functions for data
    preprocessing, integrative analysis, and visualization.
    
    Args:
        
    Attributes:
        raw_data(list): List of raw data matrices, one per experiment/dataset (genes by cells)
        norm_data(list): List of normalized matrices (genes by cells)
        scale_data(list): List of scaled matrices (cells by genes)
#        cell_data(): Dataframe of cell attributes across all datasets (nrows equal to total number
        cells across all datasets)
#        var.genes(): Subset of informative genes shared across datasets to be used in matrix
        factorization
        H(list): Cell loading factors (one matrix per dataset, dimensions cells by k)
        H_norm(): Normalized cell loading factors (cells across all datasets combined into single
        matrix)
        W(): Shared gene loading factors (k by genes)
        V(list): Dataset-specific gene loading factors (one matrix per dataset, dimensions k by genes)
        tsne_coords(): Matrix of 2D coordinates obtained from running t-SNE on H.norm or H matrices
        alignment_clusters(): Initial joint cluster assignments from shared factor alignment
        clusters(): Joint cluster assignments for cells
        snf(list): List of values associated with shared nearest factor matrix for use in clustering and
        alignment (out_summary contains edge weight information between cell combinations)
        agg_data(list): Data aggregated within clusters
        parameters(list): List of parameters used throughout analysis
        version(): Version of package used to create object
    """
    
    __slots__ = ('raw_data', 'norm_data', 'scale_data', 'cell_data', 'var_genes',
                 'H', 'H_norm', 'W', 'V', 'tsne_coords', 'alignment_clusters',
                 'clusters', 'agg_data', 'parameters', 'snf', 'version')
    
    def __init__(self, raw_data):
        self.raw_data = raw_data
        
                
    """
    @property
    def norm_data(self):
        return self._norm_data
    @norm_data.setter
    def norm_data(self, value):
        self._norm_data = value
        
    @property
    def scale_data(self):
        return self._scale_data
    @scale_data.setter
    def scale_data(self, value):
        self.scale_data = value
   """
       
    @property
    def cell_data(self):
        return self._cell_data
    @cell_data.setter
    def cell_data(self, value):
        self.cell_data = value
        
    def show(self):
        print("An object of class liger with {} datasets and {} total cells.".format(len(self.raw_data),len(self.cell_data)))

   
#######################################################################################
#### Data Preprocessing

# Read 10X alignment data (including V3)        
def read10X(sample_dirs, sample_names, merge = True, num_cells = None, min_umis = 0,
            use_filtered = False, reference = None, data_type = "rna"):
    pass


def createLiger(raw_data, 
                make_sparse = True, 
                take_gene_union = False,
                remove_missing = True):
    """ Create a liger object. 
    
    This function initializes a liger object with the raw data passed in. It requires a list of
    expression (or another single-cell modality) matrices (gene by cell) for at least two datasets.
    By default, it converts all passed data into sparse matrices (dgCMatrix) to reduce object size.
    It initializes cell.data with nUMI and nGene calculated for every cell.
    
    Args:
        raw_data(list): List of expression matrices (gene by cell). Should be named by dataset.
        make_sparse(bool): Whether to convert raw data into sparse matrices (default TRUE).
        take.gene.union(bool): Whether to fill out raw.data matrices with union of genes across all
        datasets (filling in 0 for missing data) (requires make.sparse=T) (default FALSE).
        remove_missing(bool): Whether to remove cells not expressing any measured genes, and genes not
        expressed in any cells (if take.gene.union = T, removes only genes not expressed in any
        dataset) (default TRUE).
        
    Return:
        liger_object(liger): object with raw.data slot set.
    
    Example:
        >>> Y = matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4, byrow = T)
        >>> Z = matrix(c(1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2), nrow = 4, byrow = T)
        >>> ligerex = createLiger(list(y_set = Y, z_set = Z))
        
    """
    if make_sparse:
        if :
            
            # Check if dimnames exist
            if :
                raise ValueError('Raw data must have both row (gene) and column (cell) names.')
            
    if 
        raise ValueError('At least one cell name is repeated across datasets; please make sure all cell names are unique.')

    if take_gene_union:
        merged_data = 
        if remove_missing:
            missing_genes = 
            if len(missing_genes) > 0:
                print()
    
    liger_object = Liger(raw_data)
    
    # remove missing cells
    if remove_missing:
        liger_object = removeMissingObs(liger_object, use_cols = True)
        
        #remove missing genes if not already merged
        if not take_gene_union:
            liger_object = removeMissingObs(liger_object, use_cols = False)
    
    # Initialize cell_data for liger_object with nUMI, nGene, and dataset
    nUMI = 
    nGene = 
    dataset = 
    liger_object.cell_data = 
    
    return liger_object


# Normalize raw datasets to column sums
def normalize(liger_object):
    return liger_object

# Select a subset of informative genes
def selectGenes(liger_object, var_thresh = 0.1, alpha_thresh = 0.99, num_genes = None,
                tol = 0.0001, combine = 'union', keep_unique = False, capitalize = False, 
                do_plot = False, cex_use = 0.3):
    datasets_use = range(1,len(liger_object.raw_data))
    pass

# Scale genes by root-mean-square across cells
def scaleNotCenter(liger_object, remove_missing = True):
    pass


def removeMissingObs(liger_object, 
                     slot_use = "raw_data", 
                     use_cols = True):
    """Remove cells/genes with no expression across any genes/cells
    
    Removes cells/genes from chosen slot with no expression in any genes or cells respectively.
    
    Args:
        liger_object(liger): object (scale.data or norm.data must be set).
        slot_use(str): The data slot to filter (takes "raw_data" and "scale_data") (default "raw_data").
        use_cols(bool): Treat each column as a cell (default TRUE).
        
    Return:
        liger_object(liger): object with modified raw_data (or chosen slot) (dataset names preserved).
        
    Example:
        >>>ligerex = removeMissingObs(ligerex)
    """
    
    return liger_object


#######################################################################################
#### Factorization

# Perform iNMF on scaled datasets
def optimizeALS_list(
  liger_object,
  k,
  value_lambda = 5.0,
  thresh = 1e-6,
  max_iters = 30,
  nrep = 1,
  H_init = None,
  W_init = None,
  V_init = None,
  rand_seed = 1,
  print_obj = False,
):
    pass


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


#######################################################################################
#### Quantile Alignment/Normalization
    
# Quantile align (normalize) factor loadings
def quantile_norm(liger_object, quantiles = 50, ref_dataset = None, min_cells = 20, knn_k = 20, 
                  dims_use = None, do_center = False, max_sample = 1000, eps = 0.9, refine_knn = True):
    pass

# Louvain algorithm for community detection
def louvainCluster(liger_object, resolution = 1.0, k = 20, prune = 1 / 15, eps = 0.1, nRandomStarts = 10,
                   nIterations = 100, random_seed = 1):
    pass

# Impute the query cell expression matrix
def imputeKNN(liger_object, reference, queries, knn_k = 20, weight = True, norm = True, scale = False):
    pass

# Perform Wilcoxon rank-sum test
def runWilcoxon(liger_object, compare_method, data_use = "all"):
    pass

# Linking genes to putative regulatory elements
def linkGenesAndPeaks(gene_counts, peak_counts, path_to_coords, genes_list = None, dist = "spearman", 
                      alpha = 0.05):
    pass

# Export predicted gene-pair interaction
def makeInteractTrack(corr_mat, genes_list, output_path, path_to_coords):
    pass

# Analyze biological interpretations of metagene
def runGSEA(liger_object, gene_sets = [], mat_w = True, mat_v = 0, custom_gene_sets = []):
    pass


#######################################################################################
#### Dimensionality Reduction
    
# Perform t-SNE dimensionality reduction
def runTSNE(liger_object, dims_use, use_raw = False, use_pca = False, perplexity = 30,
            theta = 0.5, method = "Rtsne", fitsne_path = None, rand_seed = 42):
    dims_use = range(1, len(liger_object.H_norm))
    pass

# Perform UMAP dimensionality reduction
def runUMAP(liger_object, dims_use, use_raw = False, k = 2, distance = "euclidean",
            n_neighbors = 10, min_dist = 0.1, rand_seed = 42):
    dims_use = range(1, len(liger_object.H_norm))
    pass


#######################################################################################
#### Metrics

# Calculate a dataset-specificity score for each factor
def calcDatasetSpecificity(liger_object, dataset1 = None, dataset2 = None, do_plot = True):
    pass

# Calculate agreement metric
def calcAgreement(liger_object, dr_method = "NMF", ndims = 40, k = 15, use_aligned = True,
                  rand_seed = 42, by_dataset = False):
    pass

# Calculate alignment metric
def calcAlignment(liger_object, k = None, rand_seed = 1, cells_use = None, cells_comp = None,
                  clusters_use = None, by_cell = False, by_dataset = False):
    pass

# Calculate alignment for each cluster
def calcAlignmentPerCluster(liger_object, rand_seed = 1, k = None, by_dataset = False):
    pass

# Calculate adjusted Rand index
def calcARI(liger_object, clusters_compare):
    pass

# Calculate purity
def calcPurity(liger_object, classes_compare):
    pass

# Calculate proportion mitochondrial contribution
def getProportionMito(liger_object, use_norm = False):
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

#######################################################################################
#### Marker/Cell Analysis

# Find shared and dataset-specific markers
def getFactorMarkers(liger_object, dataset1 = None, dataset2 = None, factor_share_thresh = 10,
                     dataset_specificity = None, log_fc_thresh = 1, umi_thresh = 30,
                     frac_thresh = 0, pval_thresh = 0.05, num_genes = 30, print_genes = False):
    pass

#######################################################################################
#### Conversion/Transformation

# Create a Seurat object containing the data from a liger object
# TO-DO names function
def ligerToSeurat(liger_object, nms = names(object@H), renormalize = True, use_liger_genes = True,
                  by_dataset = False):
    pass

# Create liger object from one or more Seurat objects
def seuratToLiger(liger_object, combined_seurat = False, names = "use-projects", meta_var = None,
                  assays_use = None, raw_assay = "RNA", remove_missing = True, renormalize = True,
                  use_seurat_genes = True, num_hvg_info = None, use_idents = True, use_tsne = True,
                  cca_to_H = False):
    pass

# Construct a liger object with a specified subset
def subsetLiger(liger_object, clusters_use = None, cells_use = None, remove_missing = True):
    pass

# Construct a liger object organized by another feature
def reorganizeLiger(liger_object, by_feature, keep_meta = True, new_label = "orig.dataset"):
    pass

# Convert older liger object into most current version (based on class definition)
def convertOldLiger(liger_object, override_raw = False):
    pass









































































