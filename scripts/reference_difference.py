import os
import sys
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(usage='reference_difference.py -r REFERENCE -g GFF -a ANTIGENS -o OUTPUT',
                                     description='''Searches for the reference genes that are absent in the found operons.''')

    parser.add_argument('-r', '--reference', default=None, nargs=1,
                        help='Reference genome annotation')
    parser.add_argument('-g', '--gff', default=None, nargs=1,
                        help='Bakta-annotated gff file')
    parser.add_argument('-a', '--antigens', default=None, nargs=1,
					    help='O-antiigens operons table')
    parser.add_argument('-o', '--output', default=None, nargs=1,
                        help='Output GFF filename')
    return parser.parse_args()

def read_gff(path_to_gff):
    """ function to read the .gff file """
    column_names=["chromosome","source","type","start","end","score","strand","phase","attributes"] # set the column names
    gff_data = pd.read_csv(path_to_gff, sep='\t', names=column_names) # read the file
    gff_data = gff_data[gff_data.chromosome.str.startswith("#") == False] # skip the commented lines
    gff_data = gff_data.astype({'start': 'int', 'end': 'int'}) # set the int type for start and end
    return gff_data.reset_index(drop=True) # reset index so that after commentaries skipping indexation will be ok
    
if __name__ == '__main__':
    reference_file = parse_args().reference[0]
    gff_file = parse_args().gff[0]
    antigens_file = parse_args().antigens[0]
    output_file = parse_args().output

    if not os.path.isfile(reference_file):
        print('Reference file not found')
        sys.exit(1)
    if not os.path.isfile(gff_file):
        print('GFF file not found')
        sys.exit(1)
    if not os.path.isfile(antigens_file):
        print('O-antigens file not found')
        sys.exit(1)

    if output_file is None:
        output_file = 'ref_tool_difference.tsv'
    else:
        output_file = output_file[0]
        if not output_file.endswith('.tsv'):
            output_file = output_file + '.tsv'
            
    if reference_file.endswith('.gz'):
    	print('Reference file should not be gzipped')
        sys.exit(1)
    
	# annotations and operon coordinates upload
	ref_annotation = read_gff(reference_file)
	tool_annotation = pd.read_csv(gff_file, sep="\t")
	operons = pd.read_csv(antigens_file, sep="\t")\
		.rename(columns={"chr": "chromosome"})

	# o-antigen operons genes from reference annotation
	ref_annot_genes = operons\
		.merge(ref_annotation, left_on='chromosome', right_on='chromosome')\
		.query('start_x <= start_y and end_x >= end_y')\
		.drop([
		"start_x", "end_x", "strand_x", "strand_y", "type_y", "color"
		], axis=1)

	# filtering by not matching coordinates
	coords = list(set(zip(ref_annot_genes.start_y, ref_annot_genes.end_y))\
			  .difference(set(zip(tool_annotation.start, tool_annotation.end))))
	filter_conditions = pd.Series([False] * len(ref_annot_genes), index=ref_annot_genes.index)

	for start, end in coords:
		filter_conditions |= ((ref_annot_genes.start_y == start) & (ref_annot_genes.end_y == end))

	ref_tool_difference = ref_annot_genes[filter_conditions][['chromosome', 'source', 'start_y', 'end_y', 
			                                              'score', 'phase', 'attributes']]\
		.rename(columns={"start_y": "start", "end_y": "end"})
		
	ref_tool_difference.to_csv(output_file, sep='\t', index=False, header=False)
