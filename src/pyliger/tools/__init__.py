from ._metrics import calc_dataset_specificity
from ._umap import run_umap
from ._marker import get_factor_markers
from ._wilcoxon import run_wilcoxon, _wilcoxon
from ._quantile_norm import quantile_norm
from ._go import run_GO_analysis

#######################################################################################
#### Conversion/Transformation TODO

# Create a Seurat object containing the data from a liger object
# TO-DO names function
# def ligerToSeurat(liger_object, nms = names(object@H), renormalize = True, use_liger_genes = True,
#                  by_dataset = False):
#    pass

# Create liger object from one or more Seurat objects
def seuratToLiger(liger_object, combined_seurat=False, names="use-projects", meta_var=None,
                  assays_use=None, raw_assay="RNA", remove_missing=True, renormalize=True,
                  use_seurat_genes=True, num_hvg_info=None, use_idents=True, use_tsne=True,
                  cca_to_H=False):
    pass


# Construct a liger object with a specified subset
def subsetLiger(liger_object, clusters_use=None, cells_use=None, remove_missing=True):
    pass


# Construct a liger object organized by another feature
def reorganizeLiger(liger_object, by_feature, keep_meta=True, new_label="orig.dataset"):
    pass


# Convert older liger object into most current version (based on class definition)
def convertOldLiger(liger_object, override_raw=False):
    pass
