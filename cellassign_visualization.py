import re
import math
from pathlib import Path
import click
from math import ceil
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

from utils import parse_markers_type1,parse_markers_type2,parse_markers_by_title

# rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.lib.ggplot2 as ggplot2

def cellratio_ggplot2(celltype_df,x,y,fill,x_label,y_label,title,legend_name,ofile):

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_celltype_df = ro.conversion.py2rpy(celltype_df)
        gp = ggplot2.ggplot(r_celltype_df)
        pp = (gp +
            ggplot2.aes_string(x=x,y=y,fill=fill) +
            ggplot2.geom_bar(position="fill", stat="identity") +
            ggplot2.labs(x="",y=y_label,title=title) +
            ggplot2.scale_y_continuous(expand=ro.FloatVector([0.001,0.001]))+
            ggplot2.theme_bw() +
            ggplot2.theme(**{"axis.text":ggplot2.element_text(color="black"),
                            "axis.text.x":ggplot2.element_text(angle=45,hjust=1,vjust=1),
                                "plot.title" : ggplot2.element_text(hjust = 0.5)}) +
            ggplot2.guides(fill=ggplot2.guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = legend_name))
            )
        pp.save(ofile)

# show result of cell assignment
# at sample / celltype / cluster level

# umap sample/celltype/cluster/all genes for each celltype
# tsne(optional) sample/celltype/cluster/all genes for each celltype
# cellratio sample/celltype/cluster
# dotplot celltype
# heatmap celltype
# tracksplot celltype
# violin_plot all genes for each cluster/celltype
# matrixplot
# stacked-violin plot

def marker_gene_expr_plot(adata,celltype,gene_list,odir):
    parts_num = math.ceil(len(gene_list) / 4)
    for num in range(parts_num):
        start = num * 4
        end = min((num + 1) * 4, len(gene_list))
        sub_gene_list = gene_list[start:end]
        with rc_context({'figure.figsize':(5,5)}):
            axes = sc.pl.umap(adata,color=sub_gene_list,ncols=2,show=False)
            plt.suptitle(f"cell markers for {celltype}")
            celltype_str = celltype.replace(" ","-")
            plt.savefig(f"{odir}/umap_{celltype_str}_gene_part{num+1}.png",dpi=150,bbox_inches='tight')
            plt.close()
            if "X_tsne" in adata.obsm:
                axes = sc.pl.tsne(adata,color=sub_gene_list,ncols=2,show=False)
                plt.suptitle(f"cell markers for {celltype}")
                celltype_str = celltype.replace(" ","-")
                plt.savefig(f"{odir}/tsne_{celltype_str}_gene_part{num+1}.png",dpi=150,bbox_inches='tight')
                plt.close()
                
# def marker_gene_expr_plot(adata,celltype,gene_list,odir):
#     plt.close('all')
#     parts_num = math.ceil(len(gene_list) / 4)
#     for num in range(parts_num):
#         start = num * 4
#         end = min((num + 1) * 4, len(gene_list))
#         sub_gene_list = gene_list[start:end]
#         with rc_context({'figure.figsize':(15,15)}):
#             fig = plt.figure()
#             fig.suptitle(f"cell markers for {celltype}")
#             for idx,gene in enumerate(sub_gene_list):
#                 ax = plt.subplot(2,2,idx+1)
#                 sc.pl.umap(adata,color=gene,ax=ax,show=False)
#             plt.tight_layout()
#             celltype_str = celltype.replace(" ","-")
#             fig.savefig(f"{odir}/umap_{celltype_str}_gene_part{num+1}.png",dpi=150,bbox_inches='tight')
#             plt.close()
            
#             if "X_tsne" in adata.obsm:
#                 fig = plt.figure()
#                 fig.suptitle(f"cell markers for {celltype}")
#                 for idx,gene in enumerate(sub_gene_list):
#                     ax = plt.subplot(2,2,idx+1)
#                     sc.pl.tsne(adata,color=gene,ax=ax,show=False)
#                 plt.tight_layout()
#                 celltype_str = celltype.replace(" ","-")
#                 fig.savefig(f"{odir}/tsne_{celltype_str}_gene_part{num+1}.png",dpi=150,bbox_inches='tight')
#                 plt.close()
                     
def violin_plot(adata,celltype,gene_list,cluster_col,odir,ncol=1):
    plt.close('all')
    parts_num = math.ceil(len(gene_list) / 4)
    for num in range(parts_num):
        start = num * 4
        end = min((num + 1) * 4, len(gene_list))
        sub_gene_list = gene_list[start:end]
        # nrow = ceil(len(sub_gene_list)/float(ncol))
        nrow = 4
        fig = plt.figure(figsize=(15,15))
        fig.suptitle(f"cell markers for {celltype}")
        for idx,gene in enumerate(sub_gene_list):
            ax = plt.subplot(nrow,ncol,idx+1)
            sc.pl.violin(adata, gene, groupby=cluster_col,ax=ax, show=False)
        plt.tight_layout()
        celltype_str = celltype.replace(" ","-")
        fig.savefig(f"{odir}/violinplot_{celltype_str}_part{num+1}.png")


def cellassign_plot(
                    h5ad,
                    odir,
                    anno,
                    marker_file,
                    fast_celltype_mode:bool,
                    celltype_col:str = None,
                    cluster_col:str = None,
                    sample_col:str = None,
                    color_map:str = "viridis",
                    dpi:int = 150,
                    dendrogram:bool = False):
    """\b
    celltype assignment visualization
    Args:
    anno file format:
        csv file
        Title is required
        Barcode should be included
    """
    if not Path(odir).exists():
        Path(odir).mkdir(parents=True,exist_ok=True)

    #fig config
    sc.settings.figdir = str(Path(odir) / "celltype_plot")
    sc.set_figure_params(dpi_save=dpi,color_map=color_map)
    plt.rcParams["figure.autolayout"] = True

    adata = sc.read(h5ad)
    if anno:
        anno_df = pd.read_csv(anno)
        adata.obs = adata.obs.merge(anno_df,how="left",left_index=True,right_on="Barcode")

    ## get and filter markers
    filtered_markers_dict = parse_markers_by_title(adata,marker_file)

    # auto detect cluster_col and celltype_col
    if fast_celltype_mode:
        celltype_cols = [name for name in adata.obs.columns if re.search(r"celltype_level\d+",name)]
        if not celltype_cols:
            raise Exception("--fast_celltype_mode is indicated but no column named starts with celltype_level found!")
        celltype_cols.sort(key = lambda colname: int(colname.replace('celltype_level','')),reverse=True)
        celltype_col = celltype_cols[0]
        cluster_col = "cluster_detail"

    ## celltype umap to show

    if sample_col:
        with rc_context({'figure.figsize':(7,5)}):
            sc.pl.umap(adata,color=sample_col, 
                palette=sc.pl.palettes.godsnot_102,
                save="_sample_legend.png")

            if "X_tsne" in adata.obsm:
                with rc_context({'figure.figsize':(7,5)}):
                    sc.pl.tsne(adata,color=sample_col, 
                        palette=sc.pl.palettes.godsnot_102,
                        save="_sample_legend.png")
    
    plot_var_list = [(celltype_col,"celltype"),(cluster_col,"cluster")]
    for col,tag in plot_var_list:
        with rc_context({'figure.figsize':(5,5)}):
            sc.pl.umap(adata,color=col, 
                palette=sc.pl.palettes.godsnot_102,
                legend_fontsize=8, legend_fontoutline=2,
                legend_loc='on data', 
                save=f"_{tag}_label.png")
        with rc_context({'figure.figsize':(7,5)}):
            sc.pl.umap(adata,color=col, 
                palette=sc.pl.palettes.godsnot_102,
                save=f"_{tag}_legend.png")

            if "X_tsne" in adata.obsm:
                with rc_context({'figure.figsize':(5,5)}):
                    sc.pl.tsne(adata,color=col, 
                        palette=sc.pl.palettes.godsnot_102,
                        legend_fontsize=8, legend_fontoutline=2,
                        legend_loc='on data', 
                        save=f"_{tag}_label.png")
                with rc_context({'figure.figsize':(7,5)}):
                    sc.pl.tsne(adata,color=col, 
                        palette=sc.pl.palettes.godsnot_102,
                        save=f"_{tag}_legend.png")

    ## dotplot
    print(filtered_markers_dict)
    try:
        adata.obs[cluster_col] = adata.obs[cluster_col].astype("category")
    except:
        pass
    with rc_context({'figure.figsize':(40,40)}):
        if hasattr(adata.obsm,"X_pca") and dendrogram:
            sc.pl.dotplot(adata, 
                        filtered_markers_dict, 
                        groupby=cluster_col, 
                        dendrogram=True,
                        save="cellmarker_cluster.png")
            sc.pl.stacked_violin(adata, 
                                filtered_markers_dict,
                                groupby=cluster_col, 
                                dendrogram=True,
                                save="cellmarker_cluster_violin.png")
        else:
            sc.pl.dotplot(adata, 
                        filtered_markers_dict, 
                        groupby=cluster_col,
                        save="cellmarker_cluster.png")
            sc.pl.stacked_violin(adata, 
                                filtered_markers_dict,
                                groupby=cluster_col, 
                                dendrogram=False,
                                save="cellmarker_cluster_violin.png")

    ## violin plot        
    for celltype,gene_list in filtered_markers_dict.items():
            violin_odir = Path(odir) / "violin"
            violin_odir.mkdir(exist_ok=True)
            gene_exp_odir = Path(odir) / "genes_umap"
            gene_exp_odir.mkdir(exist_ok=True)
            violin_plot(adata,celltype,gene_list,cluster_col,violin_odir)
            marker_gene_expr_plot(adata,celltype,gene_list,gene_exp_odir)
            
    
    ## figure out width for heatmap and trackspot
    # adata.obs["cluster_col"].value_counts()
    cluster_per = adata.obs[cluster_col].value_counts(normalize=True).to_list()
    cluster_inch = int(11 / min(cluster_per) / dpi + 4)
    print(f"cluster_inch is {cluster_inch}")
    ## heatmap
    with rc_context({'figure.figsize':(15,cluster_inch),'ytick.labelsize':7}):
        sc.pl.heatmap(adata, 
            filtered_markers_dict, 
            groupby=cluster_col, 
            cmap='viridis', 
            dendrogram=dendrogram,
            show_gene_labels=True,
            figsize=(15,cluster_inch),
            save="_cellmarker.png")

    ##Trackspot
    with rc_context({'figure.figsize':(cluster_inch,15),'xtick.labelsize':7}):
        sc.pl.tracksplot(adata, 
            filtered_markers_dict, 
            groupby=cluster_col, 
            dendrogram=dendrogram,
            figsize=(cluster_inch,15),
            save="_cellmarker.png",
            title = f"Cell Assignment Tracksplot")

    if sample_col:
        ## cellratio
        celltype_df = adata.obs[[celltype_col,sample_col]]
        celltype_df["cells_ratio"] = 1
        cellratio_ggplot2(celltype_df,
            x=celltype_col,
            y="cells_ratio",
            x_label="Celltype",
            y_label="Cells Ratio",
            fill=sample_col,
            title="Sample Ratio for Each Celltype",
            legend_name="Sample",
            ofile=str(Path(odir) / "sample_ratio.png"))

        cellratio_ggplot2(celltype_df,
            x=sample_col,
            y="cells_ratio",
            x_label="Samples",
            y_label="Cells Ratio",
            fill=celltype_col,
            title="Celltype Ratio for Each Sample",
            legend_name="Celltype",
            ofile=str(Path(odir) / "celltype_ratio.png"))

if __name__ == "__main__":
    cellassign_plot()