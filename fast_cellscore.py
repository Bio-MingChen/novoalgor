import re
import json
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from textwrap import dedent
import click
import copy
import matplotlib
matplotlib.use('Agg')
#tree
from anytree import Node, RenderTree, LevelOrderGroupIter
from anytree.exporter import DotExporter
from anytree.search import find

import novotools.utils as tools
from utils import filter_genes_v2

def compute_dropout(col,dropout):
    cluster_dropout_rate = len(col[col == 0])/len(col)
    if cluster_dropout_rate >= dropout:
        return 0
    else:
        return 1

def exp_mean(df,marker_gene_adata,marker_gene_exp_df,groupby,dropout):
    """
    group marker_gene_exp_df by marker_gene_adata.obs[groupby]
    and compute dropout for each gene,if gene's dropout rate more than indicated, 
    its expression will be set to 0
    """
    group_marker_gene_exp_df = marker_gene_exp_df[marker_gene_adata.obs[groupby] == df[groupby][0]]
    dropout_col = group_marker_gene_exp_df.apply(compute_dropout,axis=0,dropout=dropout)
    return group_marker_gene_exp_df.mean(axis=0) * dropout_col

def get_fd(df,marker_genes):
    marker_df = df[["Gene","logfoldchanges"]][df["Gene"].isin(marker_genes)]
    marker_series = pd.Series(data=marker_df["logfoldchanges"].tolist(),index=marker_df["Gene"].tolist())
    return marker_series[marker_genes]

def compute_score(row,marker_dict,undefined_rate,correct_marker_num=False,longest_markers_num=None):
    celltype_score_dict = {}
    for celltype,markers in marker_dict.items():
        if markers:
            celltype_undefined_rate = len(row[markers][row[markers] == 0])/len(row[markers])
            if celltype_undefined_rate > undefined_rate:
                celltype_score_dict[celltype] = 0
            else:
                if correct_marker_num and longest_markers_num:
                    celltype_score_dict[celltype] = row[markers].sum() * (longest_markers_num/len(markers))
                else:
                    celltype_score_dict[celltype] = row[markers].sum()
        else:
            celltype_score_dict[celltype] = 0
            
    return pd.Series(data=celltype_score_dict,index=celltype_score_dict.keys())

def define_celltype(row,default_celltype):
    if row[row == max(row)][0] <= 0:
        return default_celltype
        # return "undefined"
    return row.index[row == max(row)][0]

def add_consistency_score(row,consistency_score_dict):
        key = str(row['celltype']) + "_" + str(row['groupby'])
        if consistency_score_dict.get(key):
            score = round(consistency_score_dict[key],2)
        else:
            score = 'NA'
        return score
    
def run_fast_celltype_algorithm(node,
                                adata,
                                groupby,
                                celltype_colname,
                                marker_dict,
                                negative_marker_dict,
                                valid_marker_dict,
                                data_dir,
                                odir,
                                lowest_foldchange,
                                use_raw=True,
                                run_basic_analysis=True,
                                n_comps=50,
                                n_top_genes=2000,
                                resolution=1,
                                diff_method="t-test",
                                max_fd=3,
                                dropout=0.7,
                                undefined_rate=0.8,
                                correct_marker_num=False,
                                cluster_method = 'louvain',
                                remove_batch_effect = False,
                                batch_col = "orig.ident",
                                batch_correct_method = "harmony",
                                ):
    """
    
    """ 
    if run_basic_analysis:
        if use_raw and adata.raw: # for children use_raw is required log-normalized data will be used to analyze
            adata = adata.raw.to_adata()
            print("use adata.raw to run...")
        else:
            print("use adata.X to run...")
            
        if hasattr(adata,"layers") and ('counts' in adata.layers):
            print('adata.layers["counts"] exists,so run flavor=seurat_v3')
            sc.pp.highly_variable_genes(adata,
                layer="counts",
                flavor="seurat_v3",
                n_top_genes=n_top_genes)
        else:
            print("run highly variable gene...")
            sc.pp.highly_variable_genes(adata,n_top_genes=n_top_genes)

        adata.raw = adata
        sc.pp.scale(adata)
        sc.tl.pca(adata, svd_solver='arpack',n_comps=n_comps)
        # block of integration
        if remove_batch_effect:
            if batch_correct_method == "harmony":
                sc.external.pp.harmony_integrate(adata,batch_col)
                sc.pp.neighbors(adata,use_rep="X_pca_harmony")
            elif batch_correct_method == "bbknn":
                sc.external.pp.bbknn(adata,batch_key=batch_col,n_pcs=n_comps)
        else:
            sc.pp.neighbors(adata)

        sc.tl.umap(adata)
        if cluster_method == "louvain":
            sc.tl.louvain(adata,resolution=resolution,key_added=groupby)
        elif cluster_method == "leiden":
            sc.tl.leiden(adata,resolution=resolution,key_added=groupby)

        sc.tl.rank_genes_groups(adata, groupby, method=diff_method)
        
    # In case of no DGE data exist
    if (not run_basic_analysis) and (not adata.uns.get("rank_genes_groups")):
        sc.tl.rank_genes_groups(adata, groupby, method=diff_method)
    longest_markers_num = max([len(i) for i in marker_dict.values()])
    marker_genes = list(set([e for i in marker_dict.values() for e in i]) | set([e for i in negative_marker_dict.values() for e in i]))
    ## maker sure use data matrix to run algorithm
    if use_raw or run_basic_analysis:
        marker_gene_adata = adata.raw.to_adata()[:,marker_genes]
    else:
        marker_gene_adata = adata[:,marker_genes]
        
    marker_gene_exp_df = marker_gene_adata.to_df()
    exp_mtx = marker_gene_adata.obs.groupby(groupby).apply(exp_mean,
                                                           marker_gene_adata=marker_gene_adata,
                                                           marker_gene_exp_df=marker_gene_exp_df,
                                                           groupby=groupby,
                                                           dropout=dropout)
    
    ## Different Gene Expreesion fold change
    DGE = adata.uns["rank_genes_groups"]
    cluster_category = marker_gene_adata.obs[groupby].cat.categories
    DGE_list = []
    for name_group,pvals_group,logfoldchanges_group in zip(DGE["names"],DGE["pvals_adj"],DGE["logfoldchanges"]):
        for i,n,p,l in zip(cluster_category,name_group,pvals_group,logfoldchanges_group):
            DGE_list.append([i,n,p,l])
    DGE_df = pd.DataFrame(DGE_list,columns=["Cluster","Gene","pvalue_adjust","logfoldchanges"])
    DGE_df["Cluster"] = DGE_df["Cluster"].astype("category")
    DGE_df["Cluster"] = DGE_df["Cluster"].cat.set_categories(cluster_category)
    # records with negative fold change and p value more than 0.05 have no effect to  result
    DGE_df.loc[DGE_df["pvalue_adjust"] > 0.05,"logfoldchanges"] = lowest_foldchange
    DGE_df.loc[DGE_df["logfoldchanges"] < 0,"logfoldchanges"] = lowest_foldchange
    DGE_df.loc[DGE_df["logfoldchanges"] > max_fd,"logfoldchanges"] = max_fd
    
    pd_mtx = DGE_df.groupby("Cluster").apply(get_fd,marker_genes=marker_genes)
    
    ## compute celltype score and assign celltype
    cellmarker_score = exp_mtx * pd_mtx
    celltype_score_mtx = cellmarker_score.apply(compute_score,
                                                axis=1,
                                                marker_dict=marker_dict,
                                                undefined_rate=undefined_rate,
                                                correct_marker_num=correct_marker_num,
                                                longest_markers_num=longest_markers_num)
    negative_celltype_score_mtx = cellmarker_score.apply(compute_score,
                                                        axis=1,
                                                        marker_dict=negative_marker_dict,
                                                        undefined_rate=1) # negative dropout should be 1
    default_celltype = "undefined" if node.name == "root" else node.name
    final_celltype_score_mtx = celltype_score_mtx - negative_celltype_score_mtx
    celltype_series = final_celltype_score_mtx.apply(define_celltype,axis=1,default_celltype=default_celltype)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(f"expression matrix is {exp_mtx}")
        print(f"foldchange matrix is {pd_mtx}")
        print(f"cellmarker_score is {cellmarker_score}")
        print(final_celltype_score_mtx)
        print(celltype_series)
    
    # compute consistency/conformity score
    cluster_to_celltype_dict = celltype_series.to_dict()
    consistency_score_dict = {}
    
    for cluster,celltype in cluster_to_celltype_dict.items():
        if marker_dict.get(celltype):
            markers = marker_dict[celltype]
            cluster_marker_gene_exp_df = marker_gene_exp_df[marker_gene_adata.obs[groupby] == cluster][markers]
            dropout_series = cluster_marker_gene_exp_df.apply(lambda col: len(col[col == 0])/len(col),axis=0)
            exp_series = cellmarker_score.loc[cluster,markers].map(lambda e:0 if e==0 else 1)
            consistency_score_dict[f"{celltype}_{cluster}"] = (dropout_series * exp_series).sum() / len(exp_series)
            cluster_markers = cellmarker_score.loc[cluster,markers].index
            exp_list = exp_series.tolist()
            valid_markers = [marker for idx,marker in enumerate(cluster_markers) if exp_list[idx]]
            if node.parent:
                child_cluster = f"{node.name}_{cluster}"
            else:
                child_cluster = cluster
            valid_marker_dict[(celltype,child_cluster)] = {
                "valid_markers":valid_markers,
                "valid_percentage":f'{len(valid_markers)}/{len(cluster_markers)}',
                "missed_markers": set(cluster_markers) - set(valid_markers),
                }
    print(json.dumps(consistency_score_dict,indent=2))
    celltype_series_index_df = celltype_series.reset_index()
    celltype_series_index_df.columns = ['groupby','celltype']
    celltype_series_index_df['consistency_score'] = celltype_series_index_df.apply(add_consistency_score,axis=1,consistency_score_dict=consistency_score_dict)
    # output data
    celltype_score_ofile = data_dir / f"{node.name.replace(' ','-')}_celltype_score.tsv"
    final_celltype_score_mtx.round(2).to_csv(str(celltype_score_ofile),sep="\t",index=True)
    celltype_series_ofile = data_dir/ f"{node.name.replace(' ','-')}_celltype.tsv"
    celltype_series_index_df.to_csv(str(celltype_series_ofile),sep="\t",index=False)
    
    # cast to adata.obs
    adata.obs[celltype_colname] = adata.obs[groupby].map(celltype_series.to_dict())

    # output child data and figures
    if node.depth != 0 :
        # adata.write(odir / f'{node.name.replace(" ","-")}_predicted.h5ad')
        sc.pl.umap(adata,
            color=[node.groupby,node.celltype],
            legend_fontsize='x-small', legend_fontoutline=1,legend_loc='on data', 
            ncols=2,save=f"_{node.name.replace(' ','-')}_predict_celltype.png")

    return adata

def iter_add_child(root,celltype_dict,celltype,parent,marker_list,negative_marker_list):
    if find(root,filter_= lambda node: node.name == celltype) is None:
        if (not parent) or parent == celltype:
            Node(celltype,parent=root,marker=marker_list,negative_marker=negative_marker_list)
        elif parent in celltype_dict:
            iter_add_child(root,celltype_dict,parent,celltype_dict[parent][0],celltype_dict[parent][1],celltype_dict[parent][2])
            parent_node = find(root,filter_= lambda node: node.name == parent)
            Node(celltype,parent=parent_node,marker=marker_list,negative_marker=negative_marker_list)
        else:
            raise Exception(f"{parent} is required for tree building!")


def parse_marker_file(marker_file,adata,max_marker_num):
    celltype_dict = {}
    with open(marker_file,'r') as indata:
        title = indata.readline()
        tp = tools.TitleParser(title)
        for line in indata:
            line_list = line.rstrip('\n').split("\t")
            celltype = tp.get_field(line_list,"Celltype",check=False)
            marker = tp.get_field(line_list,"Marker",check=False)
            if not marker:
                raise Exception(f"{celltype} has no marker found!")
            else:
                marker_list =marker.split(",")
                if max_marker_num != 'all':
                    marker_list = marker_list[:max_marker_num]
                    
            parent = tp.get_field(line_list,"Parent",check=False)
            negative_marker = tp.get_field(line_list,"Negative Marker",check=False)
            if negative_marker and negative_marker != "NA":
                negative_marker_list = negative_marker.split(",")
            else:
                negative_marker_list = []
                
            celltype_dict[celltype] = (parent,marker_list,negative_marker_list)
    # filter genes
    celltype_dict = filter_genes_v2(adata,celltype_dict,None)
    # build tree
    celltype_tree = Node("root")
    for celltype,(parent,marker_list,negative_marker_list) in  celltype_dict.items():
        iter_add_child(celltype_tree,celltype_dict,celltype,parent,marker_list,negative_marker_list)
        
    # print tree
    for pre, fill, node in RenderTree(celltype_tree):
        print("%s%s" % (pre, node.name))

    return celltype_tree

def get_marker_dict(children):
    marker_dict = {}
    negative_marker_dict = {}
    for child in children:
        marker_dict[child.name] = child.marker
        negative_marker_dict[child.name] = child.negative_marker
    
    return marker_dict,negative_marker_dict
        
def fast_celltype_v2(
    h5ad,
    marker_file,
    groupby,
    run_basic_analysis,
    use_raw,
    pcs,
    subset_pcs,
    n_top_gene,
    subset_n_top_gene,
    resolution,
    subset_resolution,
    odir,
    diff_method,
    cluster_method,
    dropout,
    lowest_foldchange,
    max_marker_num,
    max_depth,
    max_fd,
    undefined_rate,
    col_to_show_in_umap,
    correct_marker_num,
    dotsize,
    remove_batch_effect,
    batch_col,
    batch_correct_method,
    ):
    """\b
    Fast celltype algorithm version2
    """
    # create output directory
    odir = Path(odir).resolve() if odir else Path(".").resolve()
    odir.mkdir(parents=True,exist_ok=True)
    # fig_dir = Path(fig_dir) if fig_dir else odir / "figures"
    fig_dir = odir / "figures"
    sc.settings.figdir = str(fig_dir)
    sc.set_figure_params(figsize=(10,10),dpi_save=150)
    data_dir = odir / 'data'
    if not data_dir.exists():
        data_dir.mkdir(parents=True)
        
    # tree analysis
    adata = sc.read(h5ad)
    try:
        adata.uns['log1p']["base"] = None # https://bytemeta.vip/repo/scverse/scanpy/issues/2239
    except:
        pass
    adata.var_names_make_unique()
    if use_raw:
        adata = adata.raw.to_adata()

    celltype_tree = parse_marker_file(marker_file,adata,max_marker_num)
    
    # cut tree if necessay
    if max_depth and max_depth != -1:
        max_depth = 1 if max_depth < 1 else max_depth
        for pre,fill,node in RenderTree(celltype_tree):
            if node.depth == max_depth:
                node.children = []
    # valid markers
    valid_marker_dict = {}
    
    celltype_tree.adata = adata
    for level_group in LevelOrderGroupIter(celltype_tree):
        for node in level_group:
            if node.is_root:
                node.groupby = groupby
                node.celltype = f"celltype_level{node.depth}"
                marker_dict,negative_marker_dict = get_marker_dict(node.children)
                node.adata = run_fast_celltype_algorithm(
                                                        node,
                                                        node.adata,
                                                        groupby,
                                                        node.celltype,
                                                        marker_dict,
                                                        negative_marker_dict,
                                                        valid_marker_dict,
                                                        data_dir,
                                                        odir,
                                                        lowest_foldchange,
                                                        use_raw=False, # root use_raw has been handled at the top of function
                                                        run_basic_analysis=run_basic_analysis,
                                                        n_comps=pcs,
                                                        n_top_genes=n_top_gene,
                                                        resolution=resolution,
                                                        diff_method=diff_method,
                                                        max_fd=max_fd,
                                                        dropout=dropout,
                                                        undefined_rate=undefined_rate,
                                                        correct_marker_num=correct_marker_num,
                                                        cluster_method = cluster_method,
                                                        remove_batch_effect = remove_batch_effect,
                                                        batch_col = batch_col,
                                                        batch_correct_method = batch_correct_method,
                                                        )
                # create cluster
                hierarchy_cluster_dict = node.adata.obs[groupby].to_dict()
                #root celltype dict
                celltype_col_dict = {node.celltype: node.adata.obs[node.celltype].to_dict()}
                
            else:
                if not hasattr(node,"adata"):
                    parent_celltype = f"celltype_level{node.parent.depth}"
                    node.adata = node.parent.adata[node.parent.adata.obs[parent_celltype]==node.name,:]
                
                child_groupby = f"{groupby}_level{node.depth}"
                node.groupby = child_groupby
                node.celltype = f"celltype_level{node.depth}"
                
                if node.adata.shape[0] == 0:
                    node.have_data = False
                else:
                    node.have_data = True

                if node.have_data and node.children:
                    marker_dict,negative_marker_dict = get_marker_dict(node.children)
                    print(f"marker_dict is {marker_dict}")
                    print(f"negative_marker_dict is {negative_marker_dict}")
                    node.adata = run_fast_celltype_algorithm(
                                                            node,
                                                            node.adata,
                                                            child_groupby,
                                                            node.celltype,
                                                            marker_dict,
                                                            negative_marker_dict,
                                                            valid_marker_dict,
                                                            data_dir,
                                                            odir,
                                                            lowest_foldchange,
                                                            use_raw=True,
                                                            run_basic_analysis=True,
                                                            n_comps=subset_pcs,
                                                            n_top_genes=subset_n_top_gene,
                                                            resolution=subset_resolution,
                                                            diff_method=diff_method,
                                                            max_fd=max_fd,
                                                            dropout=dropout,
                                                            undefined_rate=undefined_rate,
                                                            correct_marker_num=correct_marker_num,
                                                            cluster_method = cluster_method,
                                                            )
                    # print(node.adata.obs.head())
                    # change cluster
                    child_cluster = node.adata.obs[child_groupby].to_dict()
                    for key,value in hierarchy_cluster_dict.items():
                        if key in child_cluster:
                            hierarchy_cluster_dict[key] = node.name + "-" + child_cluster[key]
                    # change celltype
                    if node.celltype not in celltype_col_dict:
                        celltype_col_dict[node.celltype] = copy.deepcopy(celltype_col_dict[parent_celltype])
                    
                    child_celltype = node.adata.obs[node.celltype].to_dict()
                    for key,value in child_celltype.items():
                        celltype_col_dict[node.celltype][key] = value
                else:
                    node.adata.obs[node.groupby] = node.adata.obs[node.parent.groupby]
                    node.adata.obs[node.celltype] = node.adata.obs[node.parent.celltype]
    
    # add hierarchy_cluster_dict to root
    hierarchy_cluster_col = pd.DataFrame.from_dict(hierarchy_cluster_dict,orient="index")
    hierarchy_cluster_col.columns = ["cluster_detail"]
    celltype_tree.adata.obs = celltype_tree.adata.obs.merge(hierarchy_cluster_col,left_index=True,right_index=True)
    # add celltype to root
    for celltype,value_dict in celltype_col_dict.items():
        if celltype != celltype_tree.celltype:
            value_col = pd.DataFrame.from_dict(value_dict,orient="index")
            value_col.columns = [celltype]
            celltype_tree.adata.obs = celltype_tree.adata.obs.merge(value_col,left_index=True,right_index=True)

    #output celltype result
    celltype_tree.adata.write(odir / "predicted_root.h5ad",compression='gzip')
    celltype_cols = [name for name in celltype_tree.adata.obs.columns if re.search(r"orig.ident|(celltype_level\d+)",name)]
    celltype_ofile = str(odir / "celltype.xls")
    out_meta = celltype_tree.adata.obs[celltype_cols + [celltype_tree.groupby,"cluster_detail"]]
    if "orig.ident" in out_meta.columns:
        out_meta.rename(columns={"orig.ident":"Sample"})
    out_meta.to_csv(celltype_ofile,index=True,index_label='Barcode')
    tree_max_depth = len([i for i in LevelOrderGroupIter(celltype_tree)])
    if max_depth == 1 or tree_max_depth == 2:
        plot_cols = [celltype_tree.groupby] + celltype_cols
    else:
        plot_cols = [celltype_tree.groupby,"cluster_detail"] + celltype_cols
    if col_to_show_in_umap:
        plot_cols += col_to_show_in_umap.split(",")
    # plot for scRNA seq pipeline
    for plot_variable in plot_cols:
        sc.pl.umap(celltype_tree.adata,
            color= plot_variable,
            palette=sc.pl.palettes.godsnot_102,
            legend_fontsize='x-small', legend_fontoutline=1,legend_loc='on data',size=dotsize,
            ncols=2,save=f"_root_{plot_variable}.png")

    # integrate all variables to one plot
    sc.pl.umap(celltype_tree.adata,
            color= plot_cols,
            palette=sc.pl.palettes.godsnot_102,
            legend_fontsize='x-small', legend_fontoutline=1,legend_loc='on data',size=dotsize,
            ncols=2,save="_root_predict_celltype.png")
    sc.pl.umap(celltype_tree.adata,
            color= plot_cols,
            palette=sc.pl.palettes.godsnot_102, size=dotsize,
            ncols=1,save="_root_predict_celltype_sidelegend.png")

    # output valid markers
    valid_marker_ofile = odir / "valid_markers.xls"
    with valid_marker_ofile.open("w") as odata:
        odata.write("Celltype\tCluster\tMarker\tScoredPercentage\tMissedMarkers\n")
        for (celltype,cluster),value_dict in valid_marker_dict.items():
            valid_markers = value_dict["valid_markers"]
            valid_percentage = value_dict["valid_percentage"]
            missed_markers = value_dict["missed_markers"]
            odata.write(f"{celltype}\t{cluster}\t{','.join(valid_markers)}\t{valid_percentage}\t{','.join(missed_markers)}\n")
    print(f"output valid markers to {valid_marker_ofile}")
    
if __name__ == "__main__":
    fast_celltype_v2()