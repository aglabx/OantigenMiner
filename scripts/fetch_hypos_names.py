import argparse

from collections import defaultdict
from Bio import Entrez


def parse_args():
    parser = argparse.ArgumentParser(description=
                                     """Part annotation of hypothetical proteins in operons. 
                                     Fetch hypothetical proteins names for WP proteins
                                     """
                                     )

    parser.add_argument('in_tsv', default=None, nargs='?')

    parser.add_argument('email', default=None, nargs='?')
    parser.add_argument('out_file', default=None, nargs='?')
    return parser.parse_args()


def fetch_gene_product(protein_id: str, email: str) -> str:
    """
    fetches gene product from ncbi
    :param protein_id: (str) protein id
    :param email: (str) email is obligatory for Entrez work
    :return: (str) product title for protein id
    """
    Entrez.email = email
    with Entrez.esummary(db="protein", id=protein_id, rettype="docsum", retmode="text") as lhandle:
        record = Entrez.read(lhandle)
        if 'Title' in record[0].keys():
            seq_title = record[0]['Title']
        else:
            seq_title = None
    return seq_title


def read_tsv(tsv_path: str) -> dict:
    """
    reads blast table into dict: {query: [hit_1, hit_2, ...]}
    :param tsv_path: (str) path to input file
    :return: (dict) of blast outputs
    """
    blast_results = defaultdict(list)
    with open(tsv_path, 'rt') as in_file:
        for line in in_file:
            entry = line.strip().split('\t')
            query = entry[0]
            hit = entry[1]
            blast_results[query].append(hit)
    return blast_results

def fetch_gene_products(blast_results: dict, email: str) -> dict:
    """
    fetches blast hits for all blast results, based on known (WP_numbers)
    :param blast_results: (dict) with correspondence between query and hits
    :param email: (str) email is obligatory for Entrez work
    :return: (dict) with correspondence between query and blast results (gene title)
    """
    funcs = {}
    for query in blast_results.keys():
        for hit in blast_results[query]:
            if hit.startswith('WP_'):
                product = fetch_gene_product(protein_id=hit, email=email)
                if product is not None:
                    funcs[query] = product
                    break
                else:
                    funcs[query] = 'hypothetical protein'
    return funcs


def write_tsv(blast_prot_titles: dict, out_file: str) -> None:
    """
    writes files into tsv file
    :param blast_prot_titles: (dict) with correspondence between query and blast results (gene title)
    :param out_file: (str) path to output tab separated file
    :return: None
    """
    with open(out_file, 'wt') as out_file:
        for entry in blast_prot_titles.items():
            out_file.write('\t'.join(entry))
            out_file.write('\n')

if __name__ == '__main__':
    in_tsv = parse_args().in_tsv
    email = parse_args().email
    out_file = parse_args().out_file

    blast_results = read_tsv(tsv_path=in_tsv)
    fetched_blast_results = fetch_gene_products(blast_results=blast_results, email=email)
    write_tsv(blast_prot_titles=fetched_blast_results, out_file=out_file)
