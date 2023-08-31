import igraph as ig
import numpy as np
from annoy import AnnoyIndex
from numba import njit
from pynndescent import NNDescent
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors


def run_knn(H, k):
    """ """
    neigh = NearestNeighbors(n_neighbors=k, radius=0, algorithm="kd_tree")
    neigh.fit(H)
    H_knn = neigh.kneighbors(H, n_neighbors=k, return_distance=False)

    return H_knn


def run_ann(H, k, num_trees=None):
    # implementation using annoy library
    num_observations = H.shape[0]

    # decide number of trees
    if num_trees is None:
        if num_observations < 100000:
            num_trees = 10
        elif num_observations < 1000000:
            num_trees = 20
        elif num_observations < 5000000:
            num_trees = 50
        else:
            num_trees = 100

    # build knn graph
    t = AnnoyIndex(k, "angular")
    for i in range(num_observations):
        t.add_item(i, H[i])
    t.build(num_trees)

    # create knn indices matrices
    H_knn = np.vstack([t.get_nns_by_vector(H[i], k) for i in range(num_observations)])
    # for i in range(num_observations):
    #    if i == 0:
    #        H_knn = np.array(t.get_nns_by_vector(H[i], k))
    #    else:
    #        H_knn = np.vstack((H_knn, t.get_nns_by_vector(H[i], k)))

    return H_knn


def run_pynndescent(H, k, num_trees=None):
    num_observations = H.shape[0]

    index = NNDescent(H)
    H_knn = index.query(H, k=k)

    return H_knn


@njit
def cluster_vote(clusts, H_knn, k):
    """"""
    for i in range(H_knn.shape[0]):
        clust_counts = {}
        for j in range(k):
            if clusts[H_knn[i, j]] not in clust_counts:
                clust_counts[clusts[H_knn[i, j]]] = 1
            else:
                clust_counts[clusts[H_knn[i, j]]] += 1

        max_clust = -1
        max_count = 0
        for key, value in clust_counts.items():
            if value > max_count:
                max_clust = key
                max_count = value
            elif value == max_count:
                if key > max_clust:
                    max_clust = key
                    max_count = value
        clusts[i] = max_clust
    return clusts


def refine_clusts(H, clusts, k, use_ann, num_trees=None):
    """helper function for refining clusters related to function quantile_norm"""
    if use_ann:
        H_knn = run_ann(H, k, num_trees)
    else:
        H_knn = run_knn(H, k)

    clusts = cluster_vote(clusts, H_knn, k)

    return clusts


# @njit
def compute_snn(knn, prune):
    """helper function to compute the SNN graph"""
    # int for indexing
    knn = knn.astype(int)

    k = knn.shape[1]
    num_cells = knn.shape[0]

    rows = np.repeat(list(range(num_cells)), k)
    columns = knn.flatten()
    data = np.repeat(1, num_cells * k)
    snn = csr_matrix((data, (rows, columns)), shape=(num_cells, num_cells))

    snn = snn @ snn.transpose()

    rows, columns = snn.nonzero()
    data = snn.data / (k + (k - snn.data))
    data[data < prune] = 0

    return csr_matrix((data, (rows, columns)), shape=(num_cells, num_cells))


def build_igraph(snn):
    """"""
    sources, targets = snn.nonzero()
    weights = snn[sources, targets]

    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph()
    g.add_vertices(snn.shape[0])
    g.add_edges(list(zip(sources, targets)))
    g.es["weight"] = weights

    return g


def _assign_cluster(cluster, liger_object):
    """helper function to assign cluster to each dataset after community detection"""
    cluster = np.ravel(cluster)
    idx = 0
    for i, adata in enumerate(liger_object.adata_list):
        liger_object.adata_list[i].obs["cluster"] = cluster[
            idx : (idx + adata.shape[0])
        ]
        idx += adata.shape[0]
    return None
