import os
import sys
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(usage='gff_fuse.py -i INPUT -o OUTPUT',
                                     description='''Fuse GFF annotations from Operon Mapper and Bakta into one file.''')

    parser.add_argument('-c', '--contigs', default=None, nargs=1,
                        help='Contigs names file')
    parser.add_argument('-b', '--bakta', default=None, nargs=1,
                        help='Bakta GFF file')
    parser.add_argument('-a', '--antigens', default=None, nargs=1,
                        help='O-antigens operons GFF file')
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
    contigs_file = parse_args().contigs[0]
    bakta_file = parse_args().bakta[0]
    antigens_file = parse_args().antigens[0]
    output_file = parse_args().output

    if not os.path.isfile(contigs_file):
        print('Contigs file not found')
        sys.exit(1)
    if not os.path.isfile(bakta_file):
        print('Bakta file not found')
        sys.exit(1)
    if not os.path.isfile(antigens_file):
        print('O-antigens file not found')
        sys.exit(1)

    if output_file is None:
        output_file = os.path.splitext(input_file)[0] + '.gff3'
    else:
        output_file = output_file[0]
        if not output_file.endswith('.gff3'):
            output_file = output_file + '.gff3'

    with open(contigs_file, 'r') as file:
        lines = [line.strip() for line in file]

    bakta_inp = read_gff(bakta_file)
    operon_mapper_inp = read_gff(antigens_file)[["chromosome", "start", "end", "attributes"]]

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

    merged.to_csv(output_file, sep='\t', index=False)
