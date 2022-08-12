# This script is used to preprocess the all immune celltypes information get from celltypist
import click

import novotools.utils as tools

@click.command()
@click.argument("infile")
@click.argument("ofile")
@click.option("--only_curated","-only",is_flag=True,help="only output curated gene")
def preprocess_marker(infile,ofile,only_curated):
    """\b
    preprocess marker file
        combine the markers from article and celltypist model
        change title names to adapt to novoalgor format
    """
    with open(infile,'r') as indata,\
        open(ofile,'w') as odata:
        title = indata.readline()
        tp = tools.TitleParser(title)
        otitle = "Celltype\tParent\tMarker\tAnnotation\n"
        odata.write(otitle)
        for line in indata:
            line_list = line.strip().split("\t")
            parent = tp.get_field(line_list,"High-hierarchy cell types")
            parent = parent.replace("/","-")
            if parent == "B-cell lineage":
                parent = "B cells"
            if parent == "Cycling cells":
                continue
            if parent == "Erythroid":
                parent = "Erythrocytes"
            celltype = tp.get_field(line_list,"Low-hierarchy cell types")
            celltype = celltype.replace("/","-")
            cell_ontology = tp.get_field(line_list,"Cell Ontology ID")
            tissue = tp.get_field(line_list,"Tissues")
            dataset = tp.get_field(line_list,"Datasets")
            curated_markers = tp.get_field(line_list,"Curated markers")
            model_top10_markers = tp.get_field(line_list,"Top 10 important genes from the CellTypist model")
            if only_curated:
                markers = [g.strip() for g in curated_markers.split(',')]
            else:
                markers = [g.strip() for g in curated_markers.split(',')] + [g.strip() for g in model_top10_markers.split(',')]
                
            marker_str = ",".join(set(markers))
            annotation = f"{cell_ontology}|{tissue}|{dataset}"
            oline = f"{celltype}\t{parent}\t{marker_str}\t{annotation}\n"
            odata.write(oline)

        odata.write("T cells\tT cells\tTRBC2,CD3D,CD3E,CD3G,IL7R,GZMK,LEF1\tAdded by panglaodb\n")

if __name__ == "__main__":
    preprocess_marker()