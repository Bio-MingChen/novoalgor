from pathlib import Path
import click
from textwrap import dedent
from fast_cellscore import fast_celltype_v2
from cellassign_visualization import cellassign_plot
from specific_marker_global_filter import global_marker_gene_filter
from specific_marker_global_search_and_filter import marker_gene_search_and_filter
from tier_specific_marker_local_filter import local_marker_gene_filter
from get_top_n_marker import get_top_genes
from predicted_precision_plot import show_precision

def output_program_only_func(kwargs):

    with open(kwargs.get('output_program_only'),'w') as odata:
        scripts = dedent(f"""
        set -eo pipefail

        novoalgor fast-celltype \\
            {{h5ad}} \\
            {{marker_file}} \\
            --groupby seurat_clusters \\
            --odir {{odir}} \\
            --pcs {kwargs['pcs']} \\
            --n_top_gene {kwargs['n_top_gene']} \\
            --resolution {kwargs['resolution']} \\
            --subset_pcs {kwargs['subset_pcs']} \\
            --subset_n_top_gene {kwargs['subset_n_top_gene']} \\
            --subset_resolution {kwargs['subset_resolution']} \\
            --lowest_foldchange {kwargs['lowest_foldchange']} \\
            --max_marker_num {kwargs['max_marker_num']} \\
            --max_depth {kwargs['max_depth']} \\
            --max_fd {kwargs['max_fd']} \\
            --dropout {kwargs['dropout']} \\
            --undefined_rate {kwargs['undefined_rate']} \\
            --diff_method {kwargs['diff_method']} {'--correct_marker_num' if kwargs['correct_marker_num'] else ''}
        
        novoalgor cellvisual \\
            {{odir}}/predicted_root.h5ad \\
            {{odir}} \\
            --marker_file {{odir}}/valid_markers.xls \\
            --fast_celltype_mode

        """)
        odata.write(scripts)
        print(f"output program to {kwargs.get('output_program_only')}")

@click.group('novoalgor')
def cli():
    """\b
    Tools for cell marker search ,celltype assignment and visualization
    """

@cli.command('fast-celltype')
@click.argument("h5ad")
@click.argument("marker_file")
@click.option("--groupby", "-g", default="louvain", show_default=True, help="cluster group label")
@click.option("--run_basic_analysis", "-run", is_flag=True, help="run basic analysis")
@click.option("--use_raw", "-raw", is_flag=True, help="use adata.raw to compute gene expression")
@click.option("--pcs", "-p", default=50, type=int, show_default=True,help="number of PCs selected")
@click.option("--n_top_gene", "-n", default=2000, type=int, show_default=True, help="number of highly variable genes selected")
@click.option("--resolution", "-r", default=1, type=float, show_default=True,help="resolution used for clustering")
@click.option("--subset_pcs", "-sp", default=50, type=int, show_default=True,help="for re-cluster, number of PCs selected")
@click.option("--subset_n_top_gene", "-sn", default=2000, type=int, show_default=True,help="for re-cluster, number of highly variable genes selected")
@click.option("--subset_resolution", "-sr", default=1, type=float, show_default=True,help="for re-cluster, resolution used for clustering")
@click.option("--odir", '-o',default=".", help="output directory default is current directory")
@click.option("--lowest_foldchange","-l",default=0.1,show_default=True,help="lowerst foldchange for each gene")
@click.option("--max_marker_num","-mm",type=int,default=300,show_default=True,help="max number cell markers to use all for no limit")
@click.option("--max_depth","-md",type=int,default=-1,help="max depth to build tree")
@click.option("--max_fd","-fd",type=int,default=3,show_default=True,help="max log2foldchange,log2fd more than this will be set to this value")
@click.option("--dropout","-d",default=0.5,show_default=True,help="the gene's expression which the number of expressed cells below dropout rate will be set to 0")
@click.option("--undefined_rate","-u",type=float,default=0.7,show_default=True,help="cluster will be set to undefined if more than indicated percentage marker genes of celltype have no expression")
# @click.option("--fig_dir", "-f", help="output directory saving figures")
@click.option("--diff_method","-m",default="wilcoxon",show_default=True,help="different analysis [logreg,t-test,wilcoxon,t-test_overestim_var]")
@click.option("--cluster_method","-cl",default="louvain",type=click.Choice(["louvain","leiden"]),show_default=True,help='cluster methods to choose')
@click.option("--col_to_show_in_umap",'-c',help='one or more columns to show in final umap plot umap_root_predict_celltype.png, comma split')
@click.option("--correct_marker_num",'-cm',is_flag=True,help='correct different marker number effect by final score mutiplies length(longest markers)/length(current_markers)')
@click.option("--output_program_only","-program",help='output run program to a shell with all arguments')
@click.option("--dotsize","-ds",default=10,type=int,show_default=True,help='dot size of umap')
@click.option("--remove_batch_effect",is_flag=True,help='whether to remove batch effect')
@click.option("--batch_correct_method",default="harmony",type=click.Choice(["harmony","bbknn"]),show_default=True,help='batch effect method')
@click.option("--batch_col",default="orig.ident",show_default=True,help='column name of batch')
def fast_cellscore(**kwargs):
    """\b
    run fast-celltype algorithm to celltype assignment
    """
    if kwargs.get('output_program_only'):
        output_program_only_func(kwargs)
    else:
        del kwargs["output_program_only"]
        fast_celltype_v2(**kwargs)


@cli.command('cellvisual')
@click.argument("h5ad")
@click.argument('odir')
@click.option("--marker_file","-m",required=True,help="file contains markers")
@click.option('--fast_celltype_mode',"-f",is_flag=True,help="celltype_col and cluster_col will be auto detect if active")
@click.option('--celltype_col',help="designate which column to use as celltype in adata.obs")
@click.option('--cluster_col',help="designate which column to use as unsupervised cluster in adata.obs")
@click.option('--sample_col',help="designate which column to use as sample in adata.obs")
@click.option("--anno","-a",help="annotation need to add to adata.obs")
@click.option('--color_map',default="viridis",# 'plasma'
    help="matplotlib color map name,default is plasma,see more https://matplotlib.org/stable/tutorials/colors/colormaps.html")
@click.option("--dpi","-d",default=150,type=int,help="dpi for figure")
@click.option("--dendrogram","-e",is_flag=True,help="dendrogram or not for dotplot and heatmap")
def cellassign_visual(**kwargs):
    """\b
    visualization for celltype assignment 
    """
    cellassign_plot(**kwargs)

@cli.command('global-filter')
@click.argument("h5ad")
@click.argument("marker_file")
@click.option("--ref_celltype",'-r',help="name of reference celltype column",default='ref_celltype',show_default=True)
@click.option("--ofile","-o",help='output filename')
@click.option("--use_raw",'-u',is_flag=True,help='use adata.raw.X to compute DGE')
@click.option("--min_foldchange",'-fd',type=float,default=1,help='min foldchange to filter marker gene')
@click.option("--min_percent",'-p',default=0.5,type=float,help='min pts to filter marker gene')
@click.option('--max_num_marker','-m',default=-1,type=int,help='max number of markers to filter,default -1,all filtered markers will be outoputed')
def global_marker_filter(**kwargs):
    """\b
    filter specific markers from input marker file globally
    """
    global_marker_gene_filter(**kwargs)


@cli.command('global-search')
@click.argument("h5ad")
@click.option("--ref_celltype",'-r',help="name of reference celltype column",default='ref_celltype',show_default=True)
@click.option("--ofile","-o",help='output filename')
@click.option("--use_raw",'-u',is_flag=True,help='use adata.raw.X to compute DGE')
@click.option("--min_foldchange",'-fd',type=float,default=1,help='min foldchange to filter marker gene')
@click.option("--min_percent",'-p',default=0.5,type=float,help='min pts to filter marker gene')
@click.option("--top_n",'-t',default=100,type=int,help='top n DGE genes is selected to search specific marker')
@click.option('--max_num_marker','-m',default=-1,type=int,help='max number of markers to filter,default -1,all filtered markers will be outoputed')
def global_marker_search_and_filter(**kwargs):
    """\b
    search and filter specific markers globally
    """
    marker_gene_search_and_filter(**kwargs)

@cli.command('local-search')
@click.argument("h5ad")
@click.option("--ref_celltype",'-r',help="name of reference celltype column",default='ref_celltype',show_default=True)
@click.option("--tree_config",'-t',help="celltype relationship to build tree")
@click.option("--opre","-o",help='output filename')
@click.option("--use_raw",'-u',is_flag=True,help='use adata.raw.X to compute DGE')
@click.option("--min_foldchange",'-fd',type=float,default=1,help='min foldchange to filter marker gene')
@click.option("--min_percent",'-p',default=0.5,type=float,help='min pts to filter marker gene')
@click.option('--max_num_marker','-m',default=-1,type=int,help='max number of markers to filter,default -1,all filtered markers will be outoputed')
def local_marker_search_and_filter(**kwargs):
    """\b
    search and filter specific markers by input celltype tree locally
    """
    local_marker_gene_filter(**kwargs)

@cli.command('get-top-gene')
@click.argument("marker_file")
@click.option("--top_n",'-t',default=10,type=int,show_default=True,help="top n marker genes to retain")
@click.option("--ofile",'-o',default='filtered_marker_file.xls',show_default=True,help="output filename")
def get_top_n_genes(**kwargs):
    """\b
    get top n markers from marker file
    """
    get_top_genes(**kwargs)

@cli.command('show-precision')
@click.argument("tag")
@click.argument("h5ad")
@click.argument("ref_col")
@click.argument("predicted_col")
@click.argument("celltype_checkfile")
@click.option("--odir",'-o',default=Path('.').resolve(),show_default=True,help='output directory')
def show_precision_cmd(**kwargs):
    '''\b
    Compute precision for prediction results and output related figures
    Args:
    <tag>: label for your dateset like CBMC,Lung
    <h5ad>: h5ad contained reference and predicted results
    <ref_col>: column name of reference celltype
    <predicted_col>: predicted column name
    <celltype_checkfile>: format like below:
        Reference	Prediction
        Eryth	Erythrocytes
        CD14 Mono	CD14 Monocytes
        Mk	Megakaryocytes
        CD34	CD34+
        DC	cDC,pDC
        CD4 T	CD4 T-cells
        CD8 T	CD8 T-cells
        CD16 Mono	CD16 Monocytes
        B	B-cells
        NK	NK cells
        pDC	pDC
    You can indicate more than one predicted celltypes by comma split like DC showed above
    '''
    show_precision(**kwargs)

if __name__ == "__main__":
    cli()