import os
import sys
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(usage='hypothetical_proteins_fix.py -g GFF -h HYPOTHETICAL -o OUTPUT',
                                     description='''Fuse GFF annotations from Operon Mapper and Bakta into one file.''')


    parser.add_argument('-g', '--gff', default=None, nargs=1,
                        help='Bakta-annotated gff file')
    parser.add_argument('-h', '--hypothetical', default=None, nargs=1,
                        help='Hypothetical proteins BLASP search results')
    parser.add_argument('-o', '--output', default=None, nargs=1,
                        help='Output GFF filename')
    return parser.parse_args()
    
    
if __name__ == '__main__':
    gff_file = parse_args().gff[0]
    hypothetical_file = parse_args().hypothetical[0]
    output_file = parse_args().output

    if not os.path.isfile(gff_file):
        print('GFF3 file not found')
        sys.exit(1)
    if not os.path.isfile(hypothetical_file):
        print('Hypothetical proteins file not found')
        sys.exit(1)

    if output_file is None:
        output_file = 'output_file_fixed_hypothetical.gff3'
    else:
        output_file = output_file[0]
        if not output_file.endsw# data upload
	o_ag_annot = pd.read_csv(gff_file, sep='\t') 
	hypothetical = pd.read_csv(hypothetical_file, sep='\t', header=None)\
		.rename(columns={0: "ID", 1: "Gene_name", 2: "Protein_name"}) 

	dict_hypothetical = hypothetical.to_dict(orient='records')

	for index, row in o_ag_annot.iterrows():
		for hypo in dict_hypothetical:
		    if hypo["ID"] in row['attributes']:
		        attrs = row['attributes'].split(";")
		        attrs[0] = f'Name={hypo["Gene_name"]}'
		        attrs[2] = f'product={hypo["Protein_name"]}'
		        o_ag_annot.at[index, 'attributes'] = ";".join(attrs)

	o_ag_annot.to_csv(output_file, sep='\t', index=False)ith('.gff3'):
            output_file = output_file + '.gff3'
            
	# data upload
	o_ag_annot = pd.read_csv(gff_file, sep='\t') 
	hypothetical = pd.read_csv(hypothetical_file, sep='\t', header=None)\
		.rename(columns={0: "ID", 1: "Gene_name", 2: "Protein_name"}) 

	dict_hypothetical = hypothetical.to_dict(orient='records')

	for index, row in o_ag_annot.iterrows():
		for hypo in dict_hypothetical:
		    if hypo["ID"] in row['attributes']:
		        attrs = row['attributes'].split(";")
		        attrs[0] = f'Name={hypo["Gene_name"]}'
		        attrs[2] = f'product={hypo["Protein_name"]}'
		        o_ag_annot.at[index, 'attributes'] = ";".join(attrs)

	o_ag_annot.to_csv(output_file, sep='\t', index=False)
