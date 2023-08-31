import numpy as np
from plotnine import (
    aes,
    geom_point,
    geom_text,
    ggplot,
    ggtitle,
    guide_legend,
    guides,
    scale_color_hue,
    theme,
    theme_classic,
    xlab,
    ylab,
)


def plot_by_dataset_and_cluster(
    liger_object,
    clusters=None,
    title=None,
    pt_size=0.3,
    text_size=10,
    do_shuffle=True,
    rand_seed=1,
    axis_labels=None,
    do_legend=True,
    legend_size=7,
    return_plots=False,
    legend_text_size=12,
):
    """Plot t-SNE coordinates of cells across datasets

    Generates two plots of all cells across datasets, one colored by dataset and one colored by
    cluster. These are useful for visually examining the alignment and cluster distributions,
    respectively. If clusters have not been set yet (quantileAlignSNF not called), will plot by
    single color for second plot. It is also possible to pass in another clustering (as long as
    names match those of cells).

    Parameters
    ----------
    liger_object : TYPE
        DESCRIPTION.
    clusters : TYPE, optional
        DESCRIPTION. The default is None.
    title : TYPE, optional
        DESCRIPTION. The default is None.
    pt_size : TYPE, optional
        DESCRIPTION. The default is 0.3.
    text_size : TYPE, optional
        DESCRIPTION. The default is 3.
    do_shuffle : TYPE, optional
        DESCRIPTION. The default is True.
    rand_seed : TYPE, optional
        DESCRIPTION. The default is 1.
    axis_labels : TYPE, optional
        DESCRIPTION. The default is None.
    do_legend : TYPE, optional
        DESCRIPTION. The default is True.
    legend_size : TYPE, optional
        DESCRIPTION. The default is 5.
    return_plots : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    # tsne_coords = [adata.obs['tsne_coords'] for adata in liger_object.adata_list]
    tsne_df = liger_object.tsne_coords
    tsne_df["Cluster"] = np.asarray(
        np.concatenate(
            [adata.obs["cluster"].to_numpy() for adata in liger_object.adata_list]
        )
    )
    tsne_df["Cluster"] = tsne_df["Cluster"].astype("category")
    tsne_df["Dataset"] = np.concatenate(
        [
            np.repeat(adata.uns["sample_name"], adata.shape[0])
            for adata in liger_object.adata_list
        ]
    )

    if do_shuffle:
        tsne_df = tsne_df.sample(frac=1, random_state=rand_seed)

    p1 = (
        ggplot(data=tsne_df, mapping=aes(x="tsne1", y="tsne2", color="Dataset"))
        + geom_point(size=pt_size)
        + guides(color=guide_legend(override_aes={"size": legend_size}))
        + scale_color_hue(h=15 / 360.0, l=0.65, s=1.0, color_space="husl")
    )

    centers = (
        tsne_df.groupby("Cluster")
        .agg(tsne1=("tsne1", "median"), tsne2=("tsne2", "median"))
        .reset_index()
    )

    p2 = (
        ggplot(data=tsne_df, mapping=aes(x="tsne1", y="tsne2", color="Cluster"))
        + geom_point(size=pt_size)
        + geom_text(
            data=centers, mapping=aes(label="Cluster"), color="black", size=text_size
        )
        + guides(color=guide_legend(override_aes={"size": legend_size}))
        + scale_color_hue(h=15 / 360.0, l=0.65, s=1.0, color_space="husl")
    )

    if title:
        p1 = p1 + ggtitle(title[0])
        p2 = p2 + ggtitle(title[1])

    if axis_labels:
        p1 = p1 + xlab(axis_labels[0]) + ylab(axis_labels[1])
        p2 = p2 + xlab(axis_labels[0]) + ylab(axis_labels[1])

    p1 = p1 + theme_classic(legend_text_size)
    p2 = p2 + theme_classic(legend_text_size)

    if not do_legend:
        p1 = p1 + theme(legend_position="none")
        p2 = p2 + theme(legend_position="none")

    if return_plots:
        return [p1, p2]
    else:
        return None
