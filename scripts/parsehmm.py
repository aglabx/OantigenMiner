import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=
                                     """Part of finding domains using hmms pipeline. 
                                     Converts \s separated file into tab separated file (tsv) and
                                     adds genes info.
                                     """
                                     )
    parser.add_argument('in_file', default=None, nargs='?')
    parser.add_argument('decoder', default=None, nargs='?')
    parser.add_argument('out_file', default=None, nargs='?')
    return parser.parse_args()


def parse_hmm(hmm_data: str) -> list:
    """
    reads hmmsearch results and convert it into tsv file
    :param hmm_data: (str) path to hmmsearch results
    :return: (list) of hmmsearch hits
    """
    lines = []
    with open(hmm_data, 'r') as table:
        for line in table:
            if not line.startswith('#'):
                lines.append(line.strip().split())
    entries = []
    for line in lines:
        entry = line[:18]
        description = ' '.join(line[18:])
        entry.append(description)
        entries.append(entry)
    return entries


def read_decoder(dec_path: str) -> dict:
    """
    extracts dict of correspondence between kegg id (ex. K00973) and gene name (ex. rfbA)
    :param dec_path: (str) path to tsv table with information about correspondence between
     kegg id (ex. K00973) and gene name (ex. rfbA)
    :return: (dict) of correspondence between kegg id (ex. K00973) and
     gene name (ex. rfbA) {kegg_id: gene_name}
    """
    kegg2gene = {}
    with open(dec_path, 'rt') as tsv_file:
        for line in tsv_file:
            if not line.startswith('kegg'):
                entry = line.strip().split('\t')
                kegg_id = entry[0]
                gene_name = entry[2]
                kegg2gene[kegg_id] = gene_name
    return kegg2gene


def write_tsv(entries: list, out_path: str) -> None:
    """
    writes entries into file in .tsv format (tab-separated table)
    :param entries: (list) of filtered hmmsearch hits (see prettify_tsv)
    :param out_path: (str) path to output file
    :return: None
    """
    with open(out_path, 'wt') as out_file:
        for entry in entries:
            out_file.write('\t'.join(entry))
            out_file.write('\n')


def prettify_tsv(entries: list, decoder: dict) -> list:
    """
    adds gene names to table and removes "useless columns"
    :param entries: (list) of hmmsearch hits
    :param decoder: (dict) of correspondence between kegg id (ex. K00973) and
     gene name (ex. rfbA) {kegg_id: gene_name}
    :return: (list) of filtered hmmsearch hits
    """
    filtered_entries = []
    for entry in entries:
        filtered_entry = [entry[0], entry[2], decoder[entry[2]], entry[4]]
        filtered_entries.append(filtered_entry)
    return filtered_entries


if __name__ == '__main__':
    in_file = parse_args().in_file
    decoder_tsv = parse_args().decoder
    out_file = parse_args().out_file

    kegg2gene = read_decoder(dec_path=decoder_tsv)

    hmms = parse_hmm(hmm_data=in_file)
    hmms_prettified = prettify_tsv(entries=hmms, decoder=kegg2gene)

    write_tsv(entries=hmms_prettified, out_path=out_file)