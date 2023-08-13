import os
import sys
import argparse
import pandas as pd
import genomenotebook as gn

pd.options.mode.chained_assignment = None  # default='warn'


def parse_args():
    parser = argparse.ArgumentParser(usage='extract_operons_from_operonmapper.py -i INPUT -o OUTPUT',
                                     description='''Extract operons data from Operon Mapper operons output list''')

    parser.add_argument('-i', '--input', default=None, nargs=1,
                        help='Input file (Operon Mapper GGF3-formated output)')
    parser.add_argument('-o', '--output', default=None, nargs=1,
                        help='tab-separated operons data filename')
    parser.add_argument('-c', '--color', default='lightblue', nargs='?',
                        help='color to plot operons')
    return parser.parse_args()


def get_operon_coords(one_operon: pd.DataFrame) -> pd.Series:
    if one_operon['strand'].iloc[0] == '+':
        start = one_operon['start'].min()
        end = one_operon['end'].max()
    else:
        start = one_operon['start'].max()
        end = one_operon['end'].min()

    return pd.Series({'start': start, 'end': end,
                      'strand': one_operon['strand'].iloc[0],
		      'chr': one_operon['seq_id'].iloc[0],
                      'N_genes': len(one_operon),
                      })


def extract_operons(operons_genes: pd.DataFrame, color: str) -> pd.DataFrame:
    operons = operons_genes.groupby('operon').apply(get_operon_coords).reset_index()
    operons['left'] = operons[['start', 'end']].min(axis=1)
    operons['right'] = operons[['start', 'end']].max(axis=1)
    operons['color'] = color
    operons['type'] = 'region'
    return operons

if __name__ == '__main__':
    input_file = parse_args().input[0]
    output_file = parse_args().output
    color = parse_args().color

    if not os.path.isfile(input_file):
        print('Input file not found')
        sys.exit(1)

    if output_file is None:
        output_file = os.path.splitext(input_file)[0] + '.tsv'
    else:
        output_file = output_file[0]
        if not output_file.endswith('.tsv'):
            output_file = output_file + '.tsv'

    operons_gff = gn.parse_gff(input_file)
    operons = extract_operons(operons_gff, color)
    operons.to_csv(output_file, sep='\t', index=False)
    print('Operons data written to', output_file, 'file')
