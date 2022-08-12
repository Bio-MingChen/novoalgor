import scanpy as sc
import pandas as pd
import numpy as np
import click


def compute_specific_markers(celltype,candidate_markers,DGE_df,min_foldchange,min_percent,max_num_marker,celltype_len):
    """
    if gene expression in one group has min foldchange to all other groups
    and p value less than 0.05, mark it as specific marker else remove it
    """
    filtered_gene_list = []
    for gene in candidate_markers:
        marker_df = DGE_df.query(f'group == "{celltype}" and names == "{gene}"')
        # marker_df = DGE_df[(DGE_df.group == celltype) & (DGE_df.names == gene)]
        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
            # print(marker_df)
        if marker_df.shape[0] == celltype_len - 1:
            is_specific_marker = marker_df.apply(
                lambda row: (row['logfoldchanges']>=min_foldchange or pd.isna(row['logfoldchanges'])) and row['pvals_adj']<0.05 and row['pct_nz_group'] >= min_percent,axis=1
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
    

def get_one_to_one_DGE(adata,ref_celltype,use_raw,log2fc_min):
    DGE_df_list = []
    celltypes = pd.unique(adata.obs[ref_celltype])
    for celltype in celltypes:
        print(f"run {celltype} DGE...")
        sc.tl.rank_genes_groups(adata,
                                groupby=ref_celltype,
                                use_raw=use_raw,
                                groups='all',
                                reference=celltype,
                                pts=True)
        ref_df = sc.get.rank_genes_groups_df(adata,group=None,pval_cutoff=0.05,log2fc_min=log2fc_min)
        ref_df['reference'] = celltype
        DGE_df_list.append(ref_df)
    DGE_df = pd.concat(DGE_df_list,axis=0)
    return DGE_df

def get_raw_DGE(adata,ref_celltype,use_raw,log2fc_min):
    sc.tl.rank_genes_groups(adata,
                            groupby=ref_celltype,
                            use_raw=use_raw,
                            groups='all',
                            reference='rest',
                            pts=True)
    raw_DGE_df = sc.get.rank_genes_groups_df(adata,group=None,pval_cutoff=0.05,log2fc_min=log2fc_min)
    return raw_DGE_df

def marker_gene_search_and_filter(h5ad,ofile,use_raw,ref_celltype,min_foldchange,min_percent,top_n,max_num_marker):
    adata = sc.read(h5ad)
    DGE_df = get_one_to_one_DGE(adata,ref_celltype,use_raw,min_foldchange)
    all_ref_celltypes = pd.unique(adata.obs[ref_celltype])
    celltype_len = len(all_ref_celltypes)
    filtered_markers_dict = {}
    print(all_ref_celltypes)
    raw_DGE_df = get_raw_DGE(adata,ref_celltype,use_raw,min_foldchange)
    for celltype in all_ref_celltypes:
        candidate_markers = pd.unique(raw_DGE_df.query(f"group == '{celltype}'").iloc[:top_n,:]['names'])
        print(f'candidate markers: {candidate_markers}')
        specific_markers = compute_specific_markers(celltype,
                                                    candidate_markers,
                                                    DGE_df,
                                                    min_foldchange,
                                                    min_percent,
                                                    max_num_marker,
                                                    celltype_len)
        filtered_markers_dict[celltype] = specific_markers
        
    with open(ofile,'w') as odata:
        odata.write("Celltype\tMarker\tScore\n")
        for celltype,markers in filtered_markers_dict.items():
            gene_markers = [g for g,s in markers]
            odata.write(f"{celltype}\t{','.join(gene_markers)}\t{markers}\n")

if __name__ == "__main__":
    marker_gene_filter()
