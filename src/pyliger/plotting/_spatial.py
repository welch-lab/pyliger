import matplotlib.pyplot as plt
import pandas as pd
from plotnine import (
    aes,
    element_blank,
    geom_point,
    ggplot,
    labs,
    scale_color_cmap,
    theme,
    theme_classic,
    xlim,
    ylim,
)

from pyliger.plotting._utilities import get_gene_values


def plot_spatial(
    liger_object,
    gene,
    log2scale=True,
    image_key="lowres",
    alpha=0.5,
    crop=True,
    zero_color="#F5F5F5",
    return_plots=False,
):
    """

    :param liger_object:
    :param image_res:
    :param crop:
    :param zero_color:
    :param return_plots:
    :return:
    """
    figs = []
    for adata in liger_object.adata_list:
        gene_vals = get_gene_values(liger_object, gene, log2scale=log2scale)
        figs = _plot_spatial_adata(adata, gene, image_key, gene_vals)

    if return_plots:
        return figs
    else:
        return None


def _plot_spatial_adata(adata, gene, image_key, gene_vals):
    """"""
    res_dict = {"lowres": "tissue_lowres_scalef", "hires": "tissue_hires_scalef"}

    ### 1. Extract image information
    pxl_in_fullres = adata.obsm["pxl_in_fullres"]
    scale = adata.uns["image"]["scale_factors"][res_dict[image_key]]
    height, width, _ = adata.uns["image"][image_key].shape

    pt_size = (
        adata.uns["image"]["scale_factors"]["spot_diameter_fullres"] * scale * 1 / 32
    )
    df = _coordinate_transfer(pxl_in_fullres, scale, height)
    df["gene"] = gene_vals

    option = "plasma"
    zero_color = "#F5F5F5"

    ### 2. Plotting
    ggp = (
        ggplot(data=df, mapping=aes(x="x", y="y", color="gene"))
        + geom_point(size=pt_size)
        + xlim(0, width)
        + ylim(0, height)
        + labs(colour=gene)
        + theme_classic(12)
        + theme(
            axis_line=element_blank(),
            axis_text_x=element_blank(),
            axis_text_y=element_blank(),
            axis_ticks=element_blank(),
            axis_title_x=element_blank(),
            axis_title_y=element_blank(),
        )
        + scale_color_cmap(cmap_name=option, na_value=zero_color)
    )

    # attach spatial image
    fig = ggp.draw(return_ggplot=False, show=False)  # get matplotlib Figure object
    plt.close()
    ax = fig.get_axes()[0]
    ax.imshow(
        adata.uns["image"][image_key], extent=[0, width, 0, height]
    )  # add spatial image layer

    return fig


def _coordinate_transfer(pxl_in_fullres, scale, height):
    """Convert jpg coordinate into xy for plotting"""
    coords = pxl_in_fullres * scale

    coords[:, 0] = height - coords[:, 0]
    coords[:, [1, 0]] = coords[:, [0, 1]]

    df = pd.DataFrame(coords)
    df.columns = ["x", "y"]

    return df
