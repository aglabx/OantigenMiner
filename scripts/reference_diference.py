#!/usr/bin/env python3

import pandas as pd

def read_gff(path_to_gff):
    """ function to read the .gff file """
    column_names=["chromosome","source","type","start","end","score","strand","phase","attributes"] # set the column names
    gff_data = pd.read_csv(path_to_gff, sep='\t', names=column_names) # read the file
    gff_data = gff_data[gff_data.chromosome.str.startswith("#") == False] # skip the commented lines
    gff_data = gff_data.astype({'start': 'int', 'end': 'int'}) # set the int type for start and end
    return gff_data.reset_index(drop=True) # reset index so that after commentaries skipping indexation will be ok

# annotations and operon coordinates upload
ref_annotation = read_gff("GCF_000027025.1_ASM2702v1_genomic.gff") # FOR NIKITA that must be gunzipped reference genome annotation
tool_annotation = pd.read_csv("output_file.gff3", delimiter="\t") # FOR NIKITA that's gff_fuse.py combined annotation output
operons = pd.read_csv("o_antigen_operons.tsv", delimiter="\t")\
    .rename(columns={"chr": "chromosome"}) # FOR NIKITA that's your ✨fancy✨ extract_operons script output

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
    .rename(columns={"start_y": "start", "end_y": "end"}).to_csv("ref_tool_difference.tsv", sep="\t", index=False)
# FOR NIKITA: let ref_tool_difference.tsv be the specified by user filename
