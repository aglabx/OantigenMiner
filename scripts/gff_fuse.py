#!/usr/bin/env python3

import pandas as pd

def read_gff(path_to_gff):
    """ function to read the .gff file """
    column_names=["chromosome","source","type","start","end","score","strand","phase","attributes"] # set the column names
    gff_data = pd.read_csv(path_to_gff, sep='\t', names=column_names) # read the file
    gff_data = gff_data[gff_data.chromosome.str.startswith("#") == False] # skip the commented lines
    gff_data = gff_data.astype({'start': 'int', 'end': 'int'}) # set the int type for start and end
    return gff_data.reset_index(drop=True) # reset index so that after commentaries skipping indexation will be ok

# download id names from genomic data
with open('./gff_merge/contig_names.txt', 'r') as file:
    lines = [line.strip() for line in file]

# download gffs
bakta_inp = read_gff("./gff_merge/test_from_bakta.gff3")
operon_mapper_inp = read_gff("./gff_merge/o_antigen_operons.gff3")[["chromosome", "start", "end", "attributes"]]

# rename chromosome field in gff originating from bakta
contig_to_id = dict(zip(bakta_inp.chromosome.unique(), lines))
bakta_inp["chromosome"] = bakta_inp.chromosome.apply(lambda x: contig_to_id[x]) # rename chromosome field

# intersect gffs
merged = bakta_inp\
    .merge(operon_mapper_inp, left_on='chromosome', right_on='chromosome')\
    .query('start_x == start_y and end_x == end_y')\
    .drop(["start_y", "end_y"], axis=1)\
    .rename(columns={
        "start_x": "start",
        "end_x": "end"
    })

# make attributes
merged["attributes_raw"] = merged["attributes_x"].apply(lambda x: ";".join(x.split(";")[1:] + [x.split(";")[0]]))
merged["attributes"] = merged.attributes_raw + ";" + merged.attributes_y.apply(lambda x: x.split(";")[0])
merged = merged.drop(["attributes_x", "attributes_y", "attributes_raw"], axis=1)
merged.loc[merged.attributes.apply(lambda x: "transposase" in x.lower()), 'type'] = "Transposase"

# write gff file
merged.to_csv("output_file.gff3", sep='\t', index=False)
