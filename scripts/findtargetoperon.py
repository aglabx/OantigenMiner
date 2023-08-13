import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=
                                     """Part of finding o-antigen operons. 
                                    Extracts operons which have given genes of interest, found with HMM.
                                     """
                                     )
    parser.add_argument('in_gff', default=None, nargs='?')
    parser.add_argument('in_tsv', default=None, nargs='?')
    parser.add_argument('out_file', default=None, nargs='?')
    parser.add_argument('kegg_num_threshold', default=None, nargs='?')
    return parser.parse_args()


def read_tsv(tsv_path: str) -> set:
    """
    reads tsv file into set of locus_ids
    :param tsv_path: (str)
    :return: (str) of genes of interest
    """
    target_loci_ids = set()
    with open(tsv_path, 'rt') as tsv_file:
        for line in tsv_file:
            entry = line.strip().split('\t')
            target_locus_id = entry[0]
            target_loci_ids.add(target_locus_id)
    return target_loci_ids


def _split_attributes(attrs: str):
    pairs = attrs.split(';')
    kvs = (pair.split('=', maxsplit=1) for pair in pairs) 
    return {k:v for k, v in kvs}


def read_gff(gff_path: str):
    """
    reads gff file into pd dataframe
    :param gff_path: (str) path to gff file
    :return: pd.DataFrame from gff file
    """
    colnames = ['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    gff_data = pd.read_csv(gff_path, sep='\t', names=colnames, skiprows=1)
    gff_data['attribute_dict'] = gff_data['attributes'].apply(_split_attributes)
    
    norm_attribute = pd.json_normalize(gff_data.attribute_dict)
    gff_data = pd.concat([gff_data, norm_attribute], axis=1)
    # remove temp columns:
    gff_data = gff_data.drop(columns=['attribute_dict'])
    return gff_data


def extract_operons(gff_df, targets: set, thres: int):
    """
    extract operons
    :param gff_df: (pandas.DataFrame), containing information from gff file
    :param targets: (set) of loci_ids, genes, responsible for O-antigen biosynthesis
    :param thres: (int) minimal number of target genes in operon to consider operon to be the target one
    :return: (pandas.DataFrame) of operons with genes, responsible for O-antigen biosynthesis
    """
    gff_of_selected = gff_df[gff_df['locus_tag'].isin(targets)]

    grouped = gff_of_selected.groupby('operon').count()
    operons_interest = set(grouped.loc[grouped['locus_tag'] >= thres].index)
    gff_of_selected = gff_df[gff_df['operon'].isin(operons_interest)]

    return gff_of_selected


def write_gff(gff_df, out_path: str) -> None:
    """
    writes pandas.DataFrame into gff file
    :param gff_df: (pandas.DataFrame) needed to be written into gff file
    :param out_path: (str) path to output file
    :return: None
    """
    gff_df = gff_df[['seq_id', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]
    gff_df.to_csv(out_path, sep='\t', index=False, header=False)


if __name__ == '__main__':
    in_gff = parse_args().in_gff
    in_tsv = parse_args().in_tsv
    out_file = parse_args().out_file
    kegg_num_threshold = int(parse_args().kegg_num_threshold)

    targets = read_tsv(tsv_path=in_tsv)
    gff = read_gff(gff_path=in_gff)
    operons = extract_operons(gff_df = gff, targets = targets, thres=kegg_num_threshold)

    write_gff(gff_df=operons, out_path=out_file)
