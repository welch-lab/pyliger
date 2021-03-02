from ._metrics import calc_dataset_specificity
from ._umap import run_umap
from ._marker import get_factor_markers
from ._wilcoxon import run_wilcoxon, _wilcoxon



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


"""
def _remove_missing_obs(liger_object,
                        slot_use='raw_data',
                        use_rows=True):

    num_samples = len(liger_object.adata_list)

    removed = str(np.where(slot_use in ['raw_data', 'norm_data'] and use_rows == True, 'cells', 'genes'))
    expressed = str(np.where(removed == 'cells', ' any genes', ''))

    for i in range(num_samples):
        data_type = liger_object.adata_list[i].uns['sample_name']
        if slot_use == 'raw_data':
            filter_data = liger_object.adata_list[i].X
        elif slot_use == 'scale_data':
            filter_data = liger_object.adata_list[i].layers['scale_data']

        if use_rows:
            missing = np.array(np.sum(filter_data, axis=1)).flatten() == 0
        else:
            missing = np.array(np.sum(filter_data, axis=0)).flatten() == 0
        if np.sum(missing) > 0:
            print('Removing {} {} not expressing{} in {}.'.format(np.sum(missing), removed, expressed, data_type))
            if use_rows:
                # show gene name when the total of missing is less than 25
                if np.sum(missing) < 25:
                    print(liger_object.adata_list[i].obs.index[missing])
                liger_object.adata_list[i] = liger_object.adata_list[i][~missing, :].copy()
            else:
                # show cell name when the total of missing is less than 25
                if np.sum(missing) < 25:
                    print(liger_object.adata_list[i].var.index[missing])
                liger_object.adata_list[i] = liger_object.adata_list[i][:, ~missing].copy()

    return liger_object
"""