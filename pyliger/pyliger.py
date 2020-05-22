# -*- coding: utf-8 -*-
import re
import time
import warnings
import scipy.io
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
from scipy.sparse import csr_matrix, isspmatrix
from anndata import AnnData
import matplotlib.pyplot as plt



## The LIGER Class
class Liger(object):
    """ Main LIGER class
    
    The liger object is created from two or more single cell datasets. To construct a
    liger object, the user needs to provide at least two expression (or another
    single-cell modality) matrices. The class provides functions for data
    preprocessing, integrative analysis, and visualization.
    
    Args:
        
    Attributes:
        adata_list(list):
            List of AnnData objects, one per experiment/dataset (genes by cells)
            In each AnnData objects, main matrix stores raw data and two addtional
            layers store normalized data and scaled data with keys norm_data and
            scale_data respectively. 
            H(list): 
            Cell loading factors (one matrix per dataset, dimensions cells by k)
            W(): 
            Shared gene loading factors (k by genes)
            V(list): 
            Dataset-specific gene loading factors (one matrix per dataset, dimensions k by genes)
        cell_data(pd dataframe): 
            Dataframe of cell attributes across all datasets (nrows equal to total number
            cells across all datasets)
        var_genes(list): 
            Subset of informative genes shared across datasets to be used in matrix
            factorization
        H_norm(): 
            Normalized cell loading factors (cells across all datasets combined into single
            matrix)
        clusters(pd dataframe): 
            Joint cluster assignments for cells
        tsne_coords(): 
            Matrix of 2D coordinates obtained from running t-SNE on H_norm or H matrices
        alignment_clusters(): 
            Initial joint cluster assignments from shared factor alignment
        snf(list): 
            List of values associated with shared nearest factor matrix for use in clustering and
            alignment (out_summary contains edge weight information between cell combinations)
        agg_data(list): 
            Data aggregated within clusters
        parameters(list): 
            List of parameters used throughout analysis
        version(): 
            Version of package used to create object
    """
    
    __slots__ = ('adata_list', 'cell_data', 'var_genes', 'H_norm', 'clusters', 
                 'tsne_coords', 'alignment_clusters', 'agg_data', 'parameters', 
                 'snf', 'version')
    
    def __init__(self, adata_list):
        self.adata_list = adata_list
        
    def show(self):
        print("An object of class liger with {} datasets and {} total cells.".format(len(self.raw_data),len(self.cell_data)))


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
#def ligerToSeurat(liger_object, nms = names(object@H), renormalize = True, use_liger_genes = True,
#                  by_dataset = False):
#    pass

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



