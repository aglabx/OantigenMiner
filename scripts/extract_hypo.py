import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser


def parse_args():
    parser = argparse.ArgumentParser(description=
                                     """Part annotation of hypothetical proteins in operons. 
                                     Extracts sequences for hypothetical proteins from Prodigal&Bacta pipe
                                     """
                                     )
    parser.add_argument('in_gff', default=None, nargs='?')
    parser.add_argument('in_faa', default=None, nargs='?')
    parser.add_argument('out_file', default=None, nargs='?')
    return parser.parse_args()

def find_hypo(gff_path: str) -> set:
    """
    finds all hypothetical proteins IDS, based on "hypothetical protein" keyword.
    :param gff_path: (str) path to bacta resulting gff file
    :return: (set) of protein IDS, corresponding to hypothetical proteins
    """
    prot_ids = set()
    with open(gff_path, 'rt') as gff_file:
        for line in gff_file:
            if 'hypothetical protein' in line:
                attributes = line.strip().split('\t')[-1]
                attributes_dict = {i.split('=')[0]:
                                    i.split('=')[1] for i in attributes.split(';')}

                prot_id = attributes_dict['locus_tag']
                prot_ids.add(prot_id)
    return prot_ids


def read_fasta(fasta_path: str) -> dict:
    """
    reads fasta file into dict
    :param fasta_path: (str) path to fasta file
    :return: (dict) of correspondence between protID and sequence
    """
    records = {}
    with open(fasta_path, 'rt') as fasta:
        for record in SimpleFastaParser(fasta):
            seq_id = record[0].split(' ')[0]
            sequence = record[1]
            records[seq_id] = sequence
    return records


def extract_hypo(prot_ids: set, sequences: dict) -> list:
    """
    Extracts sequences for hypothetical proteins from Prodigal&Bacta pipe
    :param prot_ids: (set) of protein IDS, corresponding to hypothetical proteins
    :param sequences: (dict) of correspondence between prot_id and sequence
    :return:
    """
    targets = []
    for prot_id in prot_ids:
        target_seq = prot_id, sequences[prot_id]
        targets.append(target_seq)
    return targets


def write_fasta(sequences_entries: list, out_fa_path: str) -> None:
    """
    writes entries into fasta file
    :param sequences_entries: (list) of tuples: (seqid, sequence) to write
    :param out_fa_path: (str) path to output file
    :return: None
    """
    with open(out_fa_path, 'wt') as out_fa:
        for line in sequences_entries:
            out_fa.write(f'>{line[0]}\n{line[1]}\n')



if __name__ == '__main__':
    # gff_path = 'salmonella_output_file.gff3'
    in_gff = parse_args().in_gff
    in_faa = parse_args().in_faa
    out_file = parse_args().out_file

    hypo_ids = find_hypo(gff_path=in_gff)
    fasta_all = read_fasta(fasta_path=in_faa)
    fasta_hypos = extract_hypo(prot_ids=hypo_ids, sequences=fasta_all)
    write_fasta(sequences_entries=fasta_hypos, out_fa_path=out_file)
