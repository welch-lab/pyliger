import matplotlib as mpl
import numpy as np
import pandas as pd
from plotnine import (
    aes,
    annotate,
    coord_cartesian,
    element_blank,
    geom_point,
    ggplot,
    ggtitle,
    scale_color_cmap,
    theme,
    theme_bw,
    theme_classic,
    xlab,
    ylab,
)
from tqdm import tqdm

from pyliger.tools._marker import get_factor_markers
from pyliger.tools._metrics import calc_dataset_specificity


def plot_gene_loadings(
    liger_object,
    dataset1=None,
    dataset2=None,
    num_genes_show=12,
    num_genes=30,
    mark_top_genes=True,
    factor_share_thresh=10,
    log_fc_thresh=1,
    umi_thresh=30,
    frac_thresh=0,
    pval_thresh=0.05,
    do_spec_plot=True,
    max_val=0.1,
    pt_size=0.1,
    option="plasma_r",
    zero_color="#F5F5F5",
    return_plots=False,
    axis_labels=None,
    do_title=False,
):
    """Generate t-SNE plots and gene loading plots

    Plots t-SNE coordinates of all cells by their loadings on each factor. Underneath it displays the
    most highly loading shared and dataset-specific genes, along with the overall gene loadings
    for each dataset.

    It is recommended to call this function into a PDF due to the large number of
    plots produced.

    """
    ### Parameter setting
    if dataset1 is None or dataset2 is None:
        # dataset1 = liger_object.adata_list[0].uns['sample_name']
        # dataset2 = liger_object.adata_list[1].uns['sample_name']
        dataset1 = 0
        dataset2 = 1

    ### Extract Values
    H_aligned = liger_object.get_obsm("H_norm")
    W_orig = liger_object.W
    V1 = liger_object.V[dataset1]
    V2 = liger_object.V[dataset2]
    W = np.minimum(W_orig + V1, W_orig + V2)

    dataset_specificity = calc_dataset_specificity(
        liger_object, dataset1=dataset1, dataset2=dataset2, do_plot=do_spec_plot
    )

    factors_use = np.abs(
        dataset_specificity[2][dataset_specificity[2] <= factor_share_thresh]
    )

    markers = get_factor_markers(
        liger_object,
        dataset1=dataset1,
        dataset2=dataset2,
        factor_share_thresh=factor_share_thresh,
        num_genes=num_genes,
        log_fc_thresh=log_fc_thresh,
        umi_thresh=umi_thresh,
        frac_thresh=frac_thresh,
        pval_thresh=pval_thresh,
        dataset_specificity=dataset_specificity,
    )

    loadings_list = [V1, W, V2]

    names_list = [dataset1, "Shared", dataset2]
    tsne_coords = liger_object.tsne_coords
    return_plots = []

    ### Create plot
    for i in tqdm(factors_use):
        factorlab = "Factor " + str(i)
        tsne_df = tsne_coords.append({factorlab: H_aligned[:, i]})
        tsne_df[tsne_df[factorlab] == 0][factorlab] = np.nan
        factor_ds = factorlab + "Dataset Specificity: " + dataset_specificity[2][i]
        data_max = np.max(H_aligned[:, i])

        # plot t-SNE
        if max_val is not None:
            values = [0, max_val, 1]
        else:
            values = None

        def _rescale(x, _from):
            norm = mpl.colors.TwoSlopeNorm(
                vmin=_from[0],
                vcenter=_from[0] + max_val * (_from[1] - _from[0]),
                vmax=_from[1],
            )
            return norm(x)

        p1 = (
            ggplot(data=tsne_df, mapping=aes(x="tsne1", y="tsne2", color=factorlab))
            + geom_point(size=pt_size)
            + scale_color_cmap(cmap_name=option, na_value=zero_color, rescaler=_rescale)
            + theme_classic(12)
        )

        if axis_labels is not None:
            p1 = p1 + xlab(axis_labels[0]) + ylab(axis_labels[1])

        if do_title:
            p1 = p1 + ggtitle(factor_ds)

        # subset to specific factor and sort by p-value
        top_genes_V1 = markers[0][markers[0]["factor_num"] == i]
        top_genes_V1 = top_genes_V1.sort_values(by=["p_value"])["gene"]

        # don't sort for W
        top_genes_W = markers[1][markers[1]["factor_num"] == i]["gene"]
        top_genes_V2 = markers[2][markers[2]["factor_num"] == i]
        top_genes_V2 = top_genes_V2.sort_values(by=["p_value"])["gene"]

        top_genes_list = [top_genes_V1, top_genes_W, top_genes_V2]
        plot_list = []
        for idx, gene_list in enumerate(top_genes_list):
            # subset down to those which will be shown if sorting by p-val
            if len(gene_list) > num_genes_show:
                top_genes_list[idx] = gene_list[1:num_genes_show]

            top_genes = top_genes_list[idx]

            # make dataframe for cum gene loadings plot
            sorted = np.argsort(loadings_list[idx][:, i])

            # sort by loadings instead - still only showing num_genes_show
            # look through top num.genes in loadings
            top_loaded = sorted[len(sorted) : (len(sorted) - num_genes) : -1]
            top_genes = top_genes[top_loaded]
            if len(top_genes) == 0:
                top_genes = ["no genes"]

            gene_df = pd.DataFrame(
                {
                    "loadings": sorted,
                    "xpos": np.linspace(0, 1, num=len(sorted) + 1),
                    "top_k": 0,
                }
            )
            y_lim_text = gene_df["loadings"].max()

            # plot and annotate with top genes
            out_plot = (
                ggplot(data=gene_df, mapping=aes(x="xpos", y="loadings"))
                + geom_point(size=0.4)
                + theme_bw()
                + theme(
                    axis_ticks_direction_x=element_blank(),
                    axis_line_x=element_blank(),
                    axis_title=element_blank(),
                    axis_text_x=element_blank(),
                    panel_grid_major_x=element_blank(),
                    panel_grid_minor_x=element_blank(),
                )
                + ggtitle(names_list[idx])
                + annotate(
                    "text",
                    x=1.1,
                    y=np.linspace(y_lim_text, 0, num=num_genes_show),
                    label=top_genes,
                    hjust=0,
                    col="#8227A0",
                )
                + coord_cartesian(
                    xlim=[0, 1],  # this focuses the x-axis on the range of interest
                    clip="off",
                )
                + theme()
            )  # TODO
            if mark_top_genes:
                out_plot = out_plot + geom_point(
                    data=gene_df,
                    mapping=aes("xpos", "loadings"),
                    color="#8227A0",
                    size=0.5,
                )

            plot_list.append(out_plot)

        if not return_plots:
            return_plots[i].draw()

    if return_plots:
        return return_plots
