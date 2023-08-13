#!/usr/bin/env python3

import pandas as pd

# data upload
o_ag_annot = pd.read_csv("output_file.gff3", delimiter="\t") # FOR NIKITA: this input is the gff_fuse.py output
hypothetical = pd.read_csv("hypothetical_names.tsv", delimiter="\t", header=None)\
    .rename(columns={0: "ID", 1: "Gene_name", 2: "Protein_name"}) # FOR NIKITA: this input is Oxana's BLASTP using script output

dict_hypothetical = hypothetical.to_dict(orient='records')

for index, row in o_ag_annot.iterrows():
    for hypo in dict_hypothetical:
        if hypo["ID"] in row['attributes']:
            attrs = row['attributes'].split(";")
            attrs[0] = f'Name={hypo["Gene_name"]}'
            attrs[2] = f'product={hypo["Protein_name"]}'
            o_ag_annot.at[index, 'attributes'] = ";".join(attrs)

o_ag_annot.to_csv("output_file_fixed_hypothetical.gff3", sep='\t', index=False) # FOR NIKITA: this output is the final complete annotation