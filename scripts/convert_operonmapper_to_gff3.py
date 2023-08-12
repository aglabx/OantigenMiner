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
    parser.add_argument('-s', '--seqid',  default=None, nargs='?',
                        help='seq_id')

    return parser.parse_args()


def convert_to_gff3(data: pd.DataFrame, seq_id: str) -> pd.DataFrame:
    data['Operon'] = data['Operon'].fillna(method='ffill')
    data['Function'] = data['Function'].fillna('.')
    data = data.dropna(subset=data.columns.difference(['IdGene']))
    data.rename(columns={'PosLeft': 'start',
                                 'postRight': 'end',
                                 'Strand': 'strand',
                                 'Type': 'type',
                                 }, inplace=True)

    data['attributes'] = ('operon=' + data['Operon'].astype(int).astype(str).str.strip() + ';' +
                             'coggene=' + data['COGgene'].str.strip(' ').str.strip() + ';' +
                             'function=' + data['Function'].str.strip(' ').str.replace(' ', '_').str.replace(';', '_').str.strip() + ';' +
                             'gene_name=' + data['IdGene'].str.strip()  + ';' +
                             'locus_tag=' + data['IdGene'].str.strip()
                            )
    data['score'] = '.'
    data['source'] = 'OperonMapper'
    data['seq_id'] = seq_id
    data['phase'] = 0
    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)
    data = data[["seq_id", "source","type","start","end","score","strand","phase","attributes"]]
    return data


if __name__ == '__main__':
    input_file = parse_args().input[0]
    output_file = parse_args().output
    seq_id = parse_args().seqid

    if not os.path.isfile(input_file):
        print('Input file not found')
        sys.exit(1)

    if output_file is None:
        output_file = os.path.splitext(input_file)[0] + '.gff3'
    else:
        output_file = output_file[0]
        if not output_file.endswith('.gff3'):
            output_file = output_file + '.gff3'

    operons_data = pd.read_csv(input_file, sep='\t')
    operons_data = convert_to_gff3(operons_data, seq_id=seq_id)
    operons_data.to_csv(output_file, sep='\t', index=False, header=False)
    print('GFF3 written to', output_file, 'file')
