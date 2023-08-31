import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import vstack

"""
The pyliger main class
"""


class Liger(object):
    """Main LIGER class

    The liger object is created from two or more single cell datasets. To construct a
    liger object, the user needs to provide at least two expression (or another
    single-cell modality) matrices. The class serves as a container for results generated
    from  data preprocessing, integrative analysis, and visualization.

    Attributes:
        adata_list(list):
            List of AnnData objects, one per experiment/dataset (cells by genes)
            In each AnnData objects, main matrix stores raw data and two addtional
            layers store normalized and scaled data with keys 'norm_data' and
            'scale_data' respectively.
            H(matrix):
            Cell loading factors (one matrix per dataset, dimensions cells by k)
            W(matrix):
            Shared gene loading factors (k by genes)
            V(matrix):
            Dataset-specific gene loading factors (one matrix per dataset, dimensions k by genes)
        cell_data(pd dataframe):
            Dataframe of cell attributes across all datasets (nrows equal to total number
            cells across all datasets)
        var_genes(list):
            Subset of informative genes shared across datasets to be used in matrix
            factorization
        H_norm(pd dataframe):
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

    __slots__ = (
        "adata_list",
        "cell_data",
        "var_genes",
        "tsne_coords",
        "alignment_clusters",
        "agg_data",
        "parameters",
        "snf",
        "version",
    )

    def __init__(self, adata_list=[]):
        self.adata_list = adata_list

    # @property
    # def adata_list(self):
    #    return self._adata_list
    # @adata_list.setter
    # def adata_list(self):

    @property
    def num_samples(self):
        return len(self.adata_list)

    @property
    def num_var_genes(self):
        return len(self.var_genes)

    @property
    def sample_names(self):
        return [adata.uns["sample_name"] for adata in self.adata_list]

    @property
    def H(self):
        """num_cells x k"""
        return [adata.obsm["H"] for adata in self.adata_list]

    @property
    def V(self):
        """num_genes x k"""
        # return [adata.varm['V'][adata.uns['var_gene_idx'], :] for adata in self.adata_list]
        return [adata.varm["V"] for adata in self.adata_list]

    @property
    def W(self):
        """W is the shared component. Only return one is enough."""
        return self.adata_list[0].varm["W"]

    def show(self):
        print(
            "An object of class liger with {} datasets and {} total cells.".format(
                self.num_samples
            ),
            len(self.cell_data),
        )

    def add_adata(self, new_arrive):
        if isinstance(new_arrive, list):
            for adata in new_arrive:
                if isinstance(new_arrive, AnnData):
                    self.pending_list.append(adata)
                else:
                    print("Invalid input. Input must be list of AnnData object")
        elif isinstance(new_arrive, AnnData):
            self.pending_list.append(new_arrive)
        else:
            print("Invalid input. Input must be AnnData object")

    def save_raw(self):
        for i in range(self.num_samples):
            self.adata_list[i].raw = self.adata_list[i]

    def find_dataset_idx(self, dataset_name):
        for i in range(self.num_samples):
            if dataset_name == self.adata_list[i].uns["sample_name"]:
                return i
            else:
                continue
        return "Dataset does not exist"

    def get_data(self, set_name, dataset_use="all", combine=False, use_var=False):
        """"""
        if dataset_use == "all":
            if set_name == "raw":
                data = [adata.X for adata in self.adata_list]
            else:
                data = [adata.layers[set_name] for adata in self.adata_list]

            if combine:
                data = vstack(data)

        else:
            if set_name == "raw":
                data = self.adata_list[dataset_use].X
            else:
                data = self.adata_list[dataset_use].layers[set_name]

        return data

    def get_obs(self, obs_name, return_values=False):
        obs_values = pd.concat([adata.obs[obs_name] for adata in self.adata_list])

        if return_values:
            return obs_values.values
        else:
            return obs_values

    def return_H(self, dataset_use="all"):
        H_list = []
        for adata in self.adata_list:
            if dataset_use == "all":
                H_list.append(adata.obsm["H"])
            elif dataset_use == adata.uns["sample_name"]:
                H_list.append(adata.obsm["H"])
            else:
                continue
        return H_list

    def return_raw(self, dataset_use="all"):
        for adata in self.adata_list:
            if dataset_use == "all":
                yield adata.X
            elif dataset_use == adata.uns["sample_name"]:
                yield adata.X

    def get_varm(self, var_name, dataset_use="all"):
        """

        Args:
            var_name:
            dataset_use:

        Returns:

        """
        if var_name == "W":
            adata = self.adata_list[0]
            var_values = adata.varm[var_name][adata.uns["var_gene_idx"], :]
        else:
            if dataset_use == "all":
                var_values = np.concatenate(
                    [
                        adata.varm[var_name][adata.uns["var_gene_idx"], :]
                        for adata in self.adata_list
                    ]
                )
            else:
                adata = self.adata_list[dataset_use]
                var_values = adata.varm[var_name][adata.uns["var_gene_idx"], :]

        return var_values

    def get_gene_values(
        self,
        gene,
        use_cols=False,
        methylation_indices=None,
        log2scale=False,
        scale_factor=10000,
    ):
        """"""

        if methylation_indices is None:
            methylation_indices = []

        gene_vals_total = []
        for idx, adata in enumerate(self.adata_list):
            gene_names = self.adata_list[idx].var.index
            if gene in gene_names:
                gene_vals = np.ravel(adata[:, gene].layers["norm_data"].toarray())
            else:
                gene_vals = np.zeros(adata.shape[0], dtype=int)

            if log2scale and idx not in methylation_indices:
                gene_vals = np.log2(gene_vals * scale_factor + 1)

            gene_vals_total.append(gene_vals)

        return np.concatenate(gene_vals_total)

    def save_obsm(self, obsm_value, obsm_name, dataset_use="all"):
        if dataset_use == "all":
            dataset_use = list(range(self.num_samples))

        idx = 0
        for i in dataset_use:
            self.adata_list[i].obsm[obsm_name] = obsm_value[
                idx : (idx + self.adata_list[i].shape[0])
            ]
            idx += self.adata_list[i].shape[0]

        return None

    def get_obsm(self, obsm_name, dataset_use="all"):
        """ """
        if dataset_use == "all":
            dataset_use = list(range(self.num_samples))

        obsm_values = []
        for i in dataset_use:
            obsm_values.append(self.adata_list[i].obsm[obsm_name])

        return obsm_values

    def save(self):
        """

        :return:
        """
        pass

    def load(self):
        pass
