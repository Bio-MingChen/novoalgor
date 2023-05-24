import click
import novotools.utils as tools

# @click.command()
# @click.argument("marker_file")
# @click.option("--top_n",'-t',default=10,type=int,show_default=True,help="top n marker genes to retain")
# @click.option("--ofile",'-o',default='filtered_marker_file.xls',show_default=True,help="output filename")


def get_top_genes(marker_file, top_n, ofile):
    """\b
    filter top n marker genes
    """
    with open(marker_file, 'r') as indata,\
            open(ofile, 'w') as odata:
        title = indata.readline()
        odata.write("Celltype\tParent\tMarker\tNegative Marker\n")
        tp = tools.TitleParser(title)
        for line in indata:
            line_list = line.strip().split("\t")
            celltype = tp.get_field(line_list, 'Celltype', check=False)
            marker = tp.get_field(line_list, 'Marker', check=False)
            parent = tp.get_field(line_list, 'Parent', check=False) or celltype
            negative_marker = tp.get_field(
                line_list, 'Negative Marker', check=False) or ""
            filtered_marker = marker.split(",")[:top_n]
            odata.write(
                f"{celltype}\t{parent}\t{','.join(filtered_marker)}\t{negative_marker}\n")


if __name__ == "__main__":
    get_top_genes()
