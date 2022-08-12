import re
from collections import defaultdict
import novotools.utils as tools

def get_marker_meta(tag,line):
    tag_re = re.search(f'#{tag}:(.+)', line)
    if tag_re is None:
        raise Exception(f"{tag} is Required!")
    value = tag_re.group(1).strip().strip('"')
    return value

def filter_genes_v2(adata,markers_dict,alarm_file):
    alarms = []
    filtered_markers_dict = {}
    for k,(parent,p,n) in markers_dict.items():
        positive_filtered_genes = []
        negative_filtered_genes = []
        for gene_list,result_list in zip((p,n),(positive_filtered_genes,negative_filtered_genes)):
            for g in gene_list:
                if g not in adata.var_names:
                    print(f"{g} is not in this dataset, remove it!")
                    alarms.append(f"{g}不在数据集中，该marker基因在分析中将被移除！")
                else:
                    result_list.append(g)
        filtered_markers_dict[k] = (parent,positive_filtered_genes,negative_filtered_genes)
    if alarm_file:
        with open(alarm_file,'w') as odata:
            for alarm in alarms:
                odata.write(alarm)
    return filtered_markers_dict

def filter_genes(adata,markers_dict,alarm_file):
    alarms = []
    filtered_markers_dict = {}
    for k,v in markers_dict.items():
        filtered_genes = []
        for g in v:
            if g not in adata.var_names:
                print(f"{g} is not in this dataset, remove it!")
                alarms.append(f"{g}不在数据集中，该marker基因在分析中将被移除！")
            else:
                filtered_genes.append(g)
        filtered_markers_dict[k] = filtered_genes
    if alarm_file:
        with open(alarm_file,'w') as odata:
            for alarm in alarms:
                odata.write(alarm)
    return filtered_markers_dict

def parse_markers_by_title(adata,marker_file,alarm_file=None):
    """
    parse markers file by title name Celltype and Marker
    """
    markers_dict = defaultdict(list)

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
            #     species = tp.get_field(line_list,"Species")
            #     if species in ["Mouse","mouse"]:
            #         markers = [i.capitalize() for i in markers]

            markers_dict[celltype] = list(set(markers_dict[celltype] + markers))
    

    return filter_genes(adata,markers_dict,alarm_file)

def parse_markers_type2(adata,marker_file):
    """
    parse markers and remove genes which is not exist in adata.var_names
    Args:
    marker file format
    Species	Celltype	Marker	
    Mouse	HSC	COL14A1,DCN,Cygb,Lrat,Pdgfrb	
    Mouse	B cell	CD19,CD79A,FCMR,EBF1,MS4A1	
    Return:
    filtered markers dict
    """
    markers_dict = {}

    with open(marker_file,'r') as indata:
        title = indata.readline()
        for line in indata:
            line_list = line.strip().split('\t')
            markers = [i.strip() for i in line_list[2].split(",")]
            if line_list[0] in ["Mouse","mouse"]:
                markers = [i.capitalize() for i in markers]
            markers_dict[line_list[1]] = markers

    return filter_genes(adata,markers_dict)

def parse_markers_type1(adata,marker_file):
    markers_dict = {}
    meta_dict = {}

    with open(marker_file,'r') as indata:
        for line in indata:
            if line.startswith("#name"):
                meta_dict['name'] = get_marker_meta("name", line)
            elif line.startswith("#species"):
                meta_dict['species'] = get_marker_meta("species", line)
            elif line.startswith("#tissue"):
                meta_dict['tissue'] = get_marker_meta("tissue", line)
            else:
                celltype,markers = line.strip().split("\t")
                markers = [i.strip() for i in markers.split(",") if i]
                markers_dict[celltype] = markers
    
    return filter_genes(adata,markers_dict)