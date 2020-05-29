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
    """Main LIGER class
    
    The liger object is created from two or more single cell datasets. To construct a
    liger object, the user needs to provide at least two expression (or another
    single-cell modality) matrices. The class provides functions for data
    preprocessing, integrative analysis, and visualization.
        
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



