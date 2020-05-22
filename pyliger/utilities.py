#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from anndata import AnnData

def MergeSparseDataAll(adata_list, library_names = None):
    """ Function to merage all sparse data into a single one
    
    Function takes in a list of DGEs, with gene rownames and cell colnames, 
    and merges them into a single DGE.
    Also adds library_names to cell_names if expected to be overlap (common with 10X barcodes)
    
    Args:
        adata_list(list): 
            List of AnnData objects which store expression matrices (gene by cell).
        library_names(list):  optional, (default None)
    
    Return:
        M():
            AnnData object store expression matrix across datasets
    """
    """
    col_offset = 0
    allGenes = 
    allCells = []
    
    for i in range(len(adata_list)):
        curr = datalist[i]
        curr_s = 
        
        # Now, alter the indexes so that the two 3-column matrices can be properly merged.
        # First, make the current and full column numbers non-overlapping.
        curr_s[] = 
        
        # Update full cell names
        if library_names is not None:
            cellnames = 
        else:
            cellnames = 
            
        allCells = 
        
        # Next, change the row (gene) indexes so that they index on the union of the gene sets,
        # so that proper merging can occur.
        idx = 
        newgenescurr = 
        curr_s
        
        # Now bind the altered 3-column matrices together, and convert into a single sparse matrix.
        if not full_mat:
            full_mat = curr_s
        else:
            full_mat = 
        
        col_offset = len(allCells)
        
    M = 
    """
    merged_adata = AnnData()
    
    for adata in adata_list:
        merged_adata = merged_adata.concatenate(adata.T, join='inner')
    
    return merged_adata.T

def refine_clusts_knn(H, clusts, k, eps=0.1):
    for i in range(len(Hs)):
        clusts_H = clusters
        H_knn = 
        new_clusts = 
        
    return clusts