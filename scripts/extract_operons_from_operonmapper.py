import os
import sys
import argparse
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'


def parse_args():
    parser = argparse.ArgumentParser(usage='extract_operons_from_operonmapper.py -i INPUT -o OUTPUT',
                                     description='''Extract operons data from Operon Mapper operons output list''')

    parser.add_argument('-i', '--input', default=None, nargs=1,
                        help='Input file (Operon Mapper GGF3-formated output)')
    parser.add_argument('-o', '--output', default=None, nargs=1,
                        help='tab-separated operons data filename')
    return parser.parse_args()
