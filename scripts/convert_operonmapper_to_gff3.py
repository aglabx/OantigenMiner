import os
import sys
import argparse
import pandas as pd

pd.options.mode.chained_assignment = None # default='warn'

def parse_args():
    parser = argparse.ArgumentParser(usage='convert_operonmapper_to_gff3.py -i INPUT -o OUTPUT',
                                     description='''Convert Operon Mapper operons output list to gff3 format''')

    parser.add_argument('-i', '--input', default=None, nargs=1,
                        help='Input file (Operon Mapper result "list_of_operons" file)')
    parser.add_argument('-o', '--output', default=None, nargs=1,
                        help='GFF3 output filename')
    parser.add_argument('-s', '--seqidpath',  default='.', nargs='?',
                        help='seq_id')

    return parser.parse_args()


def convert_to_gff3(data: pd.DataFrame, seq_ids: dict) -> pd.DataFrame:
    data.rename(columns={'PosLeft': 'start',
                                 'postRight': 'end',
                                 'Strand': 'strand',
                                 'Type': 'type',
                                 }, inplace=True)

    data['Operon'] = data['Operon'].fillna(method='ffill')
    data[['Function', 'COGgene', 'strand', 'IdGene']] = data[['Function', 'COGgene', 'strand', 'IdGene']].fillna('.')
    data = data.groupby("Operon").apply(lambda group: group.iloc[1:]).reset_index(drop=True)

    data['attributes'] = ('operon=' + data['Operon'].astype(int).astype(str).str.strip() + ';' +
                             'coggene=' + data['COGgene'].str.strip(' ').str.strip() + ';' +
                             'function=' + data['Function'].str.strip(' ').str.replace(' ', '_').str.replace(';', '_').str.strip() + ';' +
                             'gene_name=' + data['IdGene'].str.strip()  + ';' +
                             'locus_tag=' + data['IdGene'].str.strip()
                            )
    data = data[~data['type'].isna()]
    data['score'] = '.'
    data['source'] = 'OperonMapper'
    data['seq_id'] = data.apply(lambda entry: seq_ids[entry['IdGene']], axis = 1)
    data['phase'] = 0
    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)
    data = data[["seq_id", "source","type","start","end","score","strand","phase","attributes"]]
    return data


def extract_protid2seqid(coords_file_path: str) -> dict:
    """
    reads information about seq_id from "OFRs_coordinates" file
    :param coords_file_path: (str) path to OperonMapper output file called "OFRs_coordinates"
    :return: (dict) of correspondence between internal OperonMapper
    output protein_id and sequence id: {prot_id1: seq_id1, ...}
    """
    prot_id2seq_id = {}
    with open(coords_file_path, 'rt') as gff_file:
        for line in gff_file:
            if not line.startswith('#'):
                entry = line.strip().split('\t')
                seqid = entry[0]
                protid = entry[-1].split(';')[0].split('=')[1]
                prot_id2seq_id[protid] = seqid
    return prot_id2seq_id


if __name__ == '__main__':
    input_file = parse_args().input[0]
    output_file = parse_args().output
    seqidpath = parse_args().seqidpath

    if not os.path.isfile(input_file):
        print('Input file not found')
        sys.exit(1)

    if output_file is None:
        output_file = os.path.splitext(input_file)[0] + '.gff3'
    else:
        output_file = output_file[0]
        if not output_file.endswith('.gff3'):
            output_file = output_file + '.gff3'

    orf_coords2seq_id = extract_protid2seqid(coords_file_path=seqidpath)

    operons_data = pd.read_csv(input_file, sep='\t')
    operons_data = convert_to_gff3(operons_data, seq_ids=orf_coords2seq_id)
    operons_data.to_csv(output_file, sep='\t', index=False, header=False)
    print('GFF3 written to', output_file, 'file')
