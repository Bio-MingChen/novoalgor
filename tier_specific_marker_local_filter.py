import scanpy as sc
import pandas as pd
import numpy as np
import click

#tree
from anytree import Node, RenderTree, LevelOrderGroupIter
from anytree.exporter import DotExporter
from anytree.search import find

import novotools.utils as tools

def iter_add_child(root,celltype_dict,celltype,parent):
    if find(root,filter_= lambda node: node.name == celltype) is None:
        if (not parent) or parent == celltype:
            Node(celltype,parent=root)
        elif parent in celltype_dict:
            iter_add_child(root,celltype_dict,parent,celltype_dict[parent])
            parent_node = find(root,filter_= lambda node: node.name == parent)
            Node(celltype,parent=parent_node)
        else:
            raise Exception(f"{parent} is required for {celltype} tree building!")


def parse_marker_file(marker_file):
    celltype_dict = {}
    with open(marker_file,'r') as indata:
        title = indata.readline()
        tp = tools.TitleParser(title)
        for line in indata:
            line_list = line.rstrip('\n').split("\t")
            celltype = tp.get_field(line_list,"Celltype",check=False)              
            parent = tp.get_field(line_list,"Parent",check=False)
            celltype_dict[celltype] = parent
    # build tree
    celltype_tree = Node("root")
    for celltype,parent in celltype_dict.items():
        iter_add_child(celltype_tree,celltype_dict,celltype,parent)
        
    # print tree
    for pre, fill, node in RenderTree(celltype_tree):
        print("%s%s" % (pre, node.name))

    return celltype_tree

def compute_specific_markers(celltype,DGE_df,min_foldchange,min_percent,max_num_marker,celltype_len):
    """
    if gene expression in one group has min foldchange to all other groups
    and p value less than 0.05, mark it as specific marker else remove it
    """
    celltype_DGE_df = DGE_df.query(f'group == "{celltype}"')
    
    if celltype_DGE_df.shape[0] == 0:
        print(f'{celltype} has no DGE found!')
        return []

    candidate_markers = pd.unique(celltype_DGE_df.names)
    filtered_gene_list = []
    for gene in candidate_markers:
        marker_df = celltype_DGE_df.query(f'names == "{gene}"')
        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            # print(marker_df)
        if marker_df.shape[0] == celltype_len - 1:
            is_specific_marker = marker_df.apply(
                lambda row: (row['logfoldchanges']>=min_foldchange or pd.isna(row['logfoldchanges'])) and row['pct_nz_group'] >= min_percent,axis=1
                ).all()
            if is_specific_marker:
                with pd.option_context('display.max_rows', None, 'display.max_columns', None):
                    print(f"find specific marker:{gene}!")
                    # print(marker_df)
                fd_score = marker_df['logfoldchanges'].mean()
                print(f'foldchange score is {fd_score}')
                filtered_gene_list.append((gene,fd_score))

    if filtered_gene_list:
        filtered_gene_list.sort(key=lambda tup: tup[1],reverse=True)
        print(filtered_gene_list)
        # filtered_gene_list = [g for (g,s) in filtered_gene_list]
    else:
        print(f'{celltype} has no DGE gene match conditions!')

    if max_num_marker > 0 and len(filtered_gene_list) > max_num_marker:
        filtered_gene_list = filtered_gene_list[:max_num_marker]

    return filtered_gene_list
    

def get_one_to_one_DGE(adata,ref_celltype,use_raw,log2fc_min,groups):
    DGE_df_list = []
    for celltype in groups:
        print(f"run {celltype} DGE...")
        sc.tl.rank_genes_groups(adata,
                                groupby=ref_celltype,
                                use_raw=use_raw,
                                groups=groups,
                                reference=celltype,
                                pts=True)
        ref_df = sc.get.rank_genes_groups_df(adata,group=None,pval_cutoff=0.05,log2fc_min=log2fc_min)
        
        if not 'group' in ref_df.columns: # one to one result
            ref_df['group'] = groups[0]
            
        ref_df['reference'] = celltype
        DGE_df_list.append(ref_df)
    DGE_df = pd.concat(DGE_df_list,axis=0)
    return DGE_df

def create_new_level(row,level_group):
    for node in level_group:
        target_node = find(node,filter_= lambda n: n.name == row)
        if target_node:
            return node.name
    return row
    # raise Exception(f"No target node found for {row}")

def local_marker_gene_filter(h5ad,ref_celltype,tree_config,opre,use_raw,min_foldchange,min_percent,max_num_marker):
    adata = sc.read(h5ad)
    celltype_tree = parse_marker_file(tree_config)
    for depth,level_group in enumerate(LevelOrderGroupIter(celltype_tree)):
        if depth != 0:
            level_colname = f'ref_celltype_{depth}'
            print(level_group)
            adata.obs[level_colname] = adata.obs[ref_celltype].apply(create_new_level,level_group=level_group)
            for node in level_group:
                node.level_colname = level_colname
    print(adata.obs.head())
    adata.obs.to_csv(f"{opre}_meta.xls",sep=',',index=True)
    filtered_markers_dict = {}
    for depth,level_group in enumerate(LevelOrderGroupIter(celltype_tree)):
        for node in level_group:
            if not node.is_root:
                local_groups = [n.name for n in node.parent.children] # siblings do not contain itself
                if len(local_groups) > 1:
                    level_colname = node.level_colname
                    
                    DGE_df = get_one_to_one_DGE(adata,level_colname,use_raw,min_foldchange,local_groups)
                    if DGE_df.shape[0] > 0:
                        celltype_len = len(local_groups)
                        for celltype in local_groups:
                            specific_markers = compute_specific_markers(celltype,
                                                                        DGE_df,
                                                                        min_foldchange,
                                                                        min_percent,
                                                                        max_num_marker,
                                                                        celltype_len)
                            if specific_markers:
                                filtered_markers_dict[celltype] = (specific_markers,local_groups)
                                
    ofile = f"{opre}_celltype_prediction.xls"
    with open(ofile,'w') as odata:
        odata.write("Celltype\tParent\tLocalGroups\tMarker\tScore\n")
        for celltype,(markers,local_groups) in filtered_markers_dict.items():
            gene_markers = [g for g,s in markers]
            celltype_node = find(celltype_tree,filter_= lambda node: node.name == celltype)
            if celltype_node.parent.is_root:
                parent = celltype
            else:
                parent = celltype_node.parent.name
            odata.write(f"{celltype}\t{parent}\t{','.join(local_groups)}\t{','.join(gene_markers)}\t{markers}\n")


if __name__ == "__main__":
    marker_gene_filter_command()
