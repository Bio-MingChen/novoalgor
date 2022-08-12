
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
from pysankey import sankey

mpl.use('Agg')

def parse_celltype_file(checkfile):
    ct_check_dict = {}
    with open(checkfile,'r',encoding='utf-8') as indata:
        title = indata.readline().rstrip('\n').split('\t')
        for line in indata:
            line_list = line.rstrip('\n').split('\t')
            ref_idx = title.index('Reference')
            predicted_idx = title.index('Prediction')
            ref_ct = line_list[ref_idx]
            predicted_ct = line_list[predicted_idx] 
            ct_check_dict[ref_ct] = [ct.strip() for ct in predicted_ct.split(',')]
    return ct_check_dict

def generate_umap_color(ct_check_dict,ct_set) -> dict:
    """
    Generate color dict by celltype relationship of reference and prediction
    Args:
        ct_check_dict: relationship of reference and prediction celltype
        ct_set: all celltypes for reference and prediction

    """
    color_dict = {}
    registered_set = set()
    umap_colors = sc.pl.palettes.godsnot_102[:len(ct_set)]
    for idx,(key,value_list) in enumerate(ct_check_dict.items()):
        color_dict[key] = umap_colors[idx]
        registered_set.add(key)
        for value in value_list:
            color_dict[value] = umap_colors[idx]
            registered_set.add(value)

    left_set = ct_set - registered_set

    if left_set:
        for idx,ct in enumerate(left_set):
            color_dict[ct] = umap_colors[len(ct_check_dict)+idx]

    return color_dict

def groupby_precision(df,predicted_col,check_dict):
    'df.name is the name of groupby'
    # print(df)
    boolean_series = df.apply(lambda row: row[predicted_col] in check_dict[df.name],axis=1)
    # boolean_series = df.apply(lambda row: print(row[predicted_col]),axis=1)
    total_len = len(boolean_series)
    right_len = len([1 for e in  boolean_series if e == True])
    return round(right_len/total_len,2)


def show_precision(tag,h5ad,ref_col,predicted_col,celltype_checkfile,odir):
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
    adata = sc.read(h5ad)
    sc.set_figure_params(figsize=(7,5))
    if isinstance(odir,str):
        odir = Path(odir)
    fig_dir = odir / "figures"
    sc.settings.figdir = str(fig_dir)

    ct_check_dict = parse_celltype_file(celltype_checkfile)
    ct_set = set()
    for ct in adata.obs[ref_col].unique():
        ct_set.add(ct)
    for ct in adata.obs[predicted_col].unique():
        ct_set.add(ct)

    color_dict = generate_umap_color(ct_check_dict,ct_set)
    print(color_dict)
    # umap plot
    sc.pl.umap(adata,color=[ref_col,predicted_col],legend_loc='on data',legend_fontsize='x-small',
        title=['Reference','fast-celltype'],palette=color_dict,
        save='_predicted.png')
    
    # only keep celltypes in celltype_checkfile to plot recall rate 
    adata = adata[adata.obs[ref_col].isin(list(ct_check_dict)),:]
    
    # precision computation and plot
    precision_df = adata.obs[[ref_col,predicted_col]].groupby(ref_col).apply(groupby_precision,predicted_col=predicted_col,check_dict=ct_check_dict).reset_index()
    precision_df.columns=[ref_col,'precision']
    total_boolean_series = adata.obs[[ref_col,predicted_col]].apply(lambda row:row[predicted_col] in ct_check_dict[row[ref_col]],axis=1)
    total_precision = len(total_boolean_series[total_boolean_series == True]) / len(total_boolean_series)

    # plot precision for each celltype and show the total precision
    fig, ax = plt.subplots()
    # mpl.style.use('ggplot')
    order = precision_df.sort_values(by='precision',ascending=False)[ref_col]

    sns.barplot(data=precision_df,y=ref_col,x='precision',order=order,orient='h',
        palette=color_dict,ax=ax)
    # ax.set_xticklabels(labels = order,rotation=90)
    ax.set(title=f'{tag} Total correctly assigned cells: {total_precision:.2%}',ylabel='',xlabel='Correctly assigned cells(%)')
    ax.set_xticks([0,0.25,0.5,0.75,1])
    ax.set_xticklabels([0,25,50,75,100])
    fig.savefig(str(fig_dir / "predicted_barplot.png"),transparent=False,dpi=150,bbox_inches='tight')
    print(f'Output figures to {fig_dir}')


    # Sankey plot
    ref_order = adata.obs[ref_col].value_counts().index.tolist()
    predicted_order = adata.obs[predicted_col].value_counts().index.tolist()
    adata.obs[ref_col] = pd.Categorical(adata.obs[ref_col],ref_order)
    adata.obs[predicted_col] = pd.Categorical(adata.obs[predicted_col],predicted_order)
    adata.obs = adata.obs.sort_values(by=[ref_col,predicted_col],ascending=True)

    ax = sankey(
    left=adata.obs[ref_col],
    right=adata.obs[predicted_col],
    colorDict=color_dict,
    fontsize=6,
    )
    ax.set(title=f'{tag} celltype prediction Sankey')
    plt.savefig(str(fig_dir /'predicted_sankey.png'),bbox_inches="tight", dpi=150)