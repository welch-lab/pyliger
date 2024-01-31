import warnings
from typing import Optional

import numpy as np
import pandas as pd
from plotnine import (
    aes,
    element_blank,
    geom_point,
    ggplot,
    ggtitle,
    labs,
    scale_color_cmap,
    scale_color_gradientn,
    theme,
    theme_classic,
    xlab,
    xlim,
    ylab,
    ylim,
)
from sklearn.preprocessing import scale

from pyliger.plotting._utilities import get_gene_values


def plot_gene(
    liger_object,
    gene,
    use_raw: bool = False,
    use_scaled: bool = False,
    scale_by="dataset",
    log2scale: Optional[bool] = None,
    methylation_indices=None,
    plot_by="dataset",
    set_dr_lims=False,
    pt_size=0.1,
    min_clip=None,
    max_clip=None,
    clip_absolute=False,
    points_only=False,
    option="plasma_r",
    cols_use=None,
    zero_color="#F5F5F5",
    axis_labels=None,
    do_legend=True,
    return_plots=False,
    keep_scale=False,
    legend_text_size=12,
):
    """Plot gene expression on dimensional reduction (t-SNE) coordinates

    Generates plot of dimensional reduction coordinates (default t-SNE) colored by expression of
    specified gene. Data can be scaled by dataset or selected feature column from cell_data (or across
    all cells). Data plots can be split by feature.

    Parameters
    ----------
    liger_object : TYPE
        Should call runTSNE before calling.
    gene
        Gene for which to plot expression.
    use_raw : TYPE, optional
        DESCRIPTION. The default is False.
    use_scaled : TYPE, optional
        DESCRIPTION. The default is False.
    scale_by : TYPE, optional
        DESCRIPTION. The default is 'dataset'.
    log2scale : TYPE, optional
        DESCRIPTION. The default is None.
    methylation_indices : TYPE, optional
        DESCRIPTION. The default is None.
    plot_by : TYPE, optional
        DESCRIPTION. The default is 'dataset'.
    set_dr_lims : TYPE, optional
        DESCRIPTION. The default is False.
    pt_size : TYPE, optional
        DESCRIPTION. The default is 0.1.
    min_clip : TYPE, optional
        DESCRIPTION. The default is None.
    max_clip : TYPE, optional
        DESCRIPTION. The default is None.
    clip_absolute : TYPE, optional
        DESCRIPTION. The default is False.
    points_only : TYPE, optional
        DESCRIPTION. The default is False.
    option : TYPE, optional
        DESCRIPTION. The default is 'plasma'.
    cols_use : TYPE, optional
        DESCRIPTION. The default is None.
    zero_color : TYPE, optional
        DESCRIPTION. The default is '#F5F5F5'.
    axis_labels : TYPE, optional
        DESCRIPTION. The default is None.
    do_legend : TYPE, optional
        DESCRIPTION. The default is True.
    return_plots : TYPE, optional
        DESCRIPTION. The default is False.
    keep_scale :
        Maintain min/max color scale across all plots when using plot.by (default FALSE)

    Returns
    -------
    None.
    """
    if plot_by != scale_by and use_scaled:
        warnings.warn(
            "Provided values for plot_by and scale_by do not match; results may not be very interpretable."
        )

    ### 1. Extract Gene Values
    if use_raw:
        if log2scale is None:
            log2scale = False

        # drop only outer level names
        gene_vals = get_gene_values(gene, data_use="raw", log2scale=log2scale)

    else:
        if log2scale is None:
            log2scale = True

        # rescale in case requested gene not highly variable
        if use_scaled:
            # check for feature
            if scale_by is not None and scale_by not in ["nUMI", "nGene", "dataset"]:
                raise ValueError(
                    "Please select existing feature in cell_data to scale_by, or add it before calling."
                )

            gene_vals = get_gene_values(liger_object, gene, log2scale=log2scale)

            # set up dataframe with groups
            gene_df = pd.DataFrame({"gene": gene_vals}, dtype=np.float64)

            if scale_by is None:
                gene_df["scale_by"] = np.repeat("none", gene_vals.shape[0])
            else:
                gene_df["scale_by"] = liger_object.get_obs(scale_by, return_values=True)

            # scale by selected feature
            gene_df1 = gene_df.groupby("scale_by")["gene"].transform(
                scale, with_mean=False
            )
            gene_vals = gene_df1["gene"]

            if log2scale:
                gene_vals = np.log2(10000 * gene_vals + 1)

        else:
            # using normalized data
            # indicate methylation indices here
            gene_vals = get_gene_values(
                liger_object,
                gene,
                methylation_indices=methylation_indices,
                log2scale=log2scale,
            )

    if np.sum(gene_vals) != 0:
        gene_vals[gene_vals == 0] = np.nan

    ### 2. Build dataframe for plotting
    # extract min and max expression values for plot scaling if keep_scale = T
    if keep_scale:
        max_exp_val = np.nanmax(gene_vals)
        min_exp_val = np.nanmin(gene_vals)

    umap_coords = np.concatenate(liger_object.get_obsm("umap_coords"))
    dr_df = pd.DataFrame(umap_coords, columns=["tsne1", "tsne2"])
    dr_df["gene"] = gene_vals

    # get dr limits for later
    lim1 = [dr_df["tsne1"].min(), dr_df["tsne1"].max()]
    lim2 = [dr_df["tsne2"].min(), dr_df["tsne2"].max()]

    if plot_by is not None:
        if plot_by not in ["nUMI", "nGene", "dataset"]:
            raise ValueError(
                "Please select existing feature in cell_data to plot_by, or add it before calling."
            )
        dr_df["plot_by"] = liger_object.get_obs(plot_by, return_values=True)
    else:
        dr_df["plot_by"] = np.repeat("none", gene_vals.shape[0])

    # expand clip values if only single provided
    num_levels = dr_df["plot_by"].nunique()
    if min_clip is None:
        min_clip = dict(zip(liger_object.sample_names, np.repeat(np.nan, num_levels)))

    if max_clip is None:
        max_clip = dict(zip(liger_object.sample_names, np.repeat(np.nan, num_levels)))

    ###!!!    #if min_clip is not None and
    ### 3. Create plot for each dataset
    p_list = {}
    for group_name, sub_df in dr_df.groupby("plot_by"):
        # maybe do quantile cutoff here
        if not clip_absolute:
            if not np.isnan(max_clip[group_name]):
                max_v = np.nanquantile(sub_df["gene"], q=max_clip[group_name])
            else:
                max_v = 1000
            if not np.isnan(min_clip[group_name]):
                min_v = np.nanquantile(sub_df["gene"], q=min_clip[group_name])
            else:
                min_v = 0
        else:
            max_v = max_clip[group_name]
            min_v = min_clip[group_name]

        sub_df.loc[(sub_df["gene"] > max_v) & (sub_df["gene"].notna()), "gene"] = max_v
        sub_df.loc[(sub_df["gene"] < min_v) & (sub_df["gene"].notna()), "gene"] = min_v

        ggp = (
            ggplot(data=sub_df, mapping=aes(x="tsne1", y="tsne2", color="gene"))
            + geom_point(size=pt_size)
            + labs(colour=gene)
        )

        if cols_use is not None:
            if keep_scale:
                ggp = ggp + scale_color_gradientn(
                    colors=cols_use,
                    na_value=zero_color,
                    limits=[min_exp_val, max_exp_val],
                )
            else:
                ggp = ggp + scale_color_gradientn(colors=cols_use, na_value=zero_color)
        else:
            if keep_scale:
                ggp = ggp + scale_color_cmap(
                    cmap_name=option,
                    na_value=zero_color,
                    limits=[min_exp_val, max_exp_val],
                )
            else:
                ggp = ggp + scale_color_cmap(cmap_name=option, na_value=zero_color)

        if set_dr_lims:
            ggp = ggp + xlim(lim1) + ylim(lim2)

        if plot_by is not None:
            base = sub_df["plot_by"].iloc[0]
        else:
            base = ""

        ggp = ggp + ggtitle(base)

        if axis_labels is not None:
            ggp = ggp + xlab(axis_labels[0]) + ylab(axis_labels[1])

        if not do_legend:
            ggp = ggp + theme(legend_position="none")

        if points_only:
            ggp = ggp + theme(
                axis_line=element_blank(),
                axis_text_x=element_blank(),
                axis_text_y=element_blank(),
                axis_ticks=element_blank(),
                axis_title_x=element_blank(),
                axis_title_y=element_blank(),
                legend_position="none",
                panel_background=element_blank(),
                panel_border=element_blank(),
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                plot_title=element_blank(),
            )
        p_list[sub_df["plot_by"].iloc[0]] = ggp + theme_classic(legend_text_size)

    # if plot_by == 'dataset':
    #    p_list = p_list[]

    if return_plots:
        return p_list
    else:
        for plot in p_list:
            plot.draw()
        return None


def plot_gene_dict(liger_object, gene_dict):
    return None


def plot_gene_spatial(
    liger_object,
    gene,
    use_raw: bool = False,
    use_scaled: bool = False,
    scale_by="dataset",
    log2scale: Optional[float] = None,
    methylation_indices=None,
    plot_by="dataset",
    set_dr_lims=False,
    pt_size=0.1,
    min_clip=None,
    max_clip=None,
    clip_absolute=False,
    points_only=False,
    option="plasma_r",
    cols_use=None,
    zero_color="#F5F5F5",
    axis_labels=None,
    do_legend=True,
    return_plots=False,
    keep_scale=False,
):
    if plot_by != scale_by and use_scaled:
        warnings.warn(
            "Provided values for plot_by and scale_by do not match; results may not be very interpretable."
        )

    ### 1. Extract Gene Values
    if use_raw:
        if log2scale is None:
            log2scale = False

        # drop only outer level names
        gene_vals = get_gene_values(gene, data_use="raw", log2scale=log2scale)

    else:
        if log2scale is None:
            log2scale = True

        # rescale in case requested gene not highly variable
        if use_scaled:
            # check for feature
            if scale_by is not None and scale_by not in ["nUMI", "nGene", "dataset"]:
                raise ValueError(
                    "Please select existing feature in cell_data to scale_by, or add it before calling."
                )

            gene_vals = get_gene_values(liger_object, gene, log2scale=log2scale)

            # set up dataframe with groups
            gene_df = pd.DataFrame({"gene": gene_vals}, dtype=np.float64)

            if scale_by is None:
                gene_df["scale_by"] = np.repeat("none", gene_vals.shape[0])
            else:
                gene_df["scale_by"] = liger_object.get_obs(scale_by, return_values=True)

            # scale by selected feature
            gene_df1 = gene_df.groupby("scale_by")["gene"].transform(
                scale, with_mean=False
            )
            gene_vals = gene_df1["gene"]

            if log2scale:
                gene_vals = np.log2(10000 * gene_vals + 1)

        else:
            # using normalized data
            # indicate methylation indices here
            gene_vals = get_gene_values(
                liger_object,
                gene,
                methylation_indices=methylation_indices,
                log2scale=log2scale,
            )

    if np.sum(gene_vals) != 0:
        gene_vals[gene_vals == 0] = np.nan

    figs = []
    for adata in liger_object.adata_list:
        figs = _plot_spatial_adata(adata, image_key)

    if return_plots:
        return figs
    else:
        return None
