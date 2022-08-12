import scanpy as sc
import pandas as pd
import numpy as np
import click
import novotools.utils as tools

def filter_genes(adata,markers_dict):
    alarms = []
    filtered_markers_dict = {}
    for k,v in markers_dict.items():
        filtered_genes = []
        for g in v:
            if g not in adata.var_names:
                print(f"{g} is not in this dataset, remove it!")
            else:
                filtered_genes.append(g)
        filtered_markers_dict[k] = filtered_genes

    return filtered_markers_dict

def parse_markers_by_title(adata,marker_file):
    """
    parse markers file by title name Celltype and Marker
    """
    markers_dict = {}

    with open(marker_file,'r') as indata:
        title = indata.readline()
        tp = tools.TitleParser(title)
        for line in indata:
            line_list = line.rstrip('\n').split('\t')
            celltype = tp.get_field(line_list,"Celltype",check=False)
            markers = tp.get_field(line_list,"Marker",check=False)
            if markers:
                markers = [i.strip() for i in markers.split(",") if i.strip()]
            # if tp.have_title("Species"):
            #     species = tp.get_field(line_list,"Species",check=False)
            #     if species in ["Mouse","mouse"]:
            #         markers = [i.capitalize() for i in markers]

            markers_dict[celltype] = markers

    return filter_genes(adata,markers_dict)

    
def compute_specific_markers(celltype,markers,DGE_df,min_foldchange,min_percent,max_num_marker):
    """
    if gene expression in one group has min foldchange to all other groups
    and p value less than 0.05, mark it as specific marker else remove it
    """
    filtered_gene_list = []
    for gene in markers:
        marker_df = DGE_df.query(f'group == "{celltype}" and names == "{gene}"')
        is_specific_marker = marker_df.apply(
            lambda row: (row['logfoldchanges']>=min_foldchange or pd.isna(row['logfoldchanges'])) and row['pvals_adj']<0.05 and row['pct_nz_group'] >= min_percent,axis=1
            ).all()
        if is_specific_marker:
            with pd.option_context('display.max_rows', None, 'display.max_columns', None):
                print(f"find specific marker:{gene}!")
                print(marker_df)
            fd_score = marker_df['logfoldchanges'].mean()
            print(f'foldchange score is {fd_score}')
            if (gene,fd_score) not in filtered_gene_list:
                filtered_gene_list.append((gene,fd_score))

    if filtered_gene_list:
        print(filtered_gene_list)
        filtered_gene_list.sort(key=lambda tup: tup[1],reverse=True)
        # filtered_gene_list = [g for g,s in filtered_gene_list]
    else:
        print(f'{celltype} has no DGE gene match conditions!')
    if max_num_marker > 0 and len(filtered_gene_list) > max_num_marker:
        filtered_gene_list = filtered_gene_list[:max_num_marker]

    return filtered_gene_list
    

def get_one_to_one_DGE(adata,ref_celltype,use_raw):
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
        ref_df = sc.get.rank_genes_groups_df(adata,group=None)
        ref_df['reference'] = celltype
        DGE_df_list.append(ref_df)
    DGE_df = pd.concat(DGE_df_list,axis=0)
    return DGE_df

def global_marker_gene_filter(h5ad,marker_file,ofile,use_raw,ref_celltype,min_foldchange,min_percent,max_num_marker):
    adata = sc.read(h5ad)
    DGE_df = get_one_to_one_DGE(adata,ref_celltype,use_raw)
    all_ref_celltypes = pd.unique(adata.obs[ref_celltype])
    markers_dict = parse_markers_by_title(adata,marker_file)
    filtered_markers_dict = {}
    print(all_ref_celltypes)
    print(list(markers_dict))
    for celltype,markers in markers_dict.items():
        if celltype in all_ref_celltypes:
            specific_markers = compute_specific_markers(celltype,markers,DGE_df,min_foldchange,min_percent,max_num_marker)
            filtered_markers_dict[celltype] = specific_markers
        # else:
        #     filtered_markers_dict[celltype] = markers
        
    # with open(ofile,'w') as odata:
    #     odata.write("Celltype\tMarker\n")
    #     for celltype,markers in filtered_markers_dict.items():
    #         odata.write(f"{celltype}\t{','.join(markers)}\n")
    with open(ofile,'w') as odata:
        odata.write("Celltype\tMarker\tScore\n")
        for celltype,markers in filtered_markers_dict.items():
            gene_markers = [g for g,s in markers]
            odata.write(f"{celltype}\t{','.join(gene_markers)}\t{markers}\n")
            
if __name__ == "__main__":
    marker_gene_filter()
