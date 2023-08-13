from findtargetoperon import read_gff, write_gff
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import math


def remove_transposaes(seq_record: SeqRecord, insertions: pd.DataFrame):
    seq = seq_record.seq

    pos = 0
    delta = 0

    deletions = []  # new index, removed string, delta_at_current_pos
    cleaned = []

    for _, insertion in insertions.iterrows():
        start = insertion["start"]
        end = insertion["end"]

        if start < pos:
            start = pos

        clean = seq[pos:start]

        cleaned.append(clean)
        transposone_seq = seq[start:end]
        assert len(transposone_seq) == end - start
        new_delta = delta + len(transposone_seq)

        deletions.append((
            start - delta, 
            transposone_seq, 
            delta,
            insertion['score'],
            insertion['strand']
        ))

        delta = new_delta
        pos = end

    cleaned.append(seq[pos : len(seq)])

    clean_seq = Seq("".join(str(s) for s in cleaned))

    sr = seq
    sc = clean_seq

    assert len(sr) - len(sc) - delta == 0

    return clean_seq, deletions


def _splice(args):
    secs = list(SeqIO.parse(args.fasta, "fasta"))
    annotation = read_gff(args.gff)
    annotation = annotation.sort_values(by="start")
    insertions = annotation["type"] == "insertion_sequence"
    annotation = annotation[insertions]

    annotation = annotation.reset_index(drop=True)

    sequences = []
    dels = []

    for seq_record in secs:
        this_seq = annotation["seq_id"] == seq_record.id
        cs, ds = remove_transposaes(seq_record, annotation[this_seq])

        s = seq_record[:1]
        s.seq = cs
        s.description = s.description + ", transposon free"

        sequences.append(s)
        dels.append(ds)

        # Just checking if we'll be able to rebuild sequence later
        reb = rebuild(cs, ds)
        assert len(seq_record.seq) == len(reb)
        assert seq_record.seq == reb

    SeqIO.write(sequences, args.output, format="fasta")
    print(f"{args.output}")
    for sequence, ds in zip(sequences, dels):
        file = f"{args.output}_{sequence.id}.rebuild.csv"
        with open(file, "w") as out:
            for idx, seq, delta, score, strand in ds:
                print(f"{idx},{seq},{delta},{score},{strand}", file=out)
        print(file)


def rebuild(seq: Seq, deletions):
    out = []
    pos = 0
    for insert_idx, insertion, *_ in deletions:
        prev = seq[pos:insert_idx]
        out.append(prev)
        out.append(insertion)
        pos = insert_idx

    out.append(seq[pos : len(seq)])

    return Seq("".join(str(s) for s in out))


def _remap_genes(frm, to, delta, annotation):
    
    start_filter = (annotation['start'] >= frm) & (annotation['start'] < to)
    annotation.loc[start_filter, 'nstart'] = annotation.loc[start_filter, 'start'] + delta
    
    end_filter = (annotation['end'] >= frm) & (annotation['end'] < to)
    annotation.loc[end_filter, 'nend'] = annotation.loc[end_filter, 'end'] + delta

    return annotation
    

def mktrns(id, start, end, score, strand, attribs):
    return {
        'seq_id': id,
        'source': 'transposon_cutter',
        'type': 'insertion_sequence',
        'nstart': start,
        'nend': end,
        'score': score,
        'strand': strand,
        'phase': 0,
        'attributes': attribs,
        'gene_name': f'TRANSP_{str(start).zfill(5)}'
    }
    

def _encode(d: dict):
    res = []
    for x, y in d.items():
        res.append(f'{x}={y}')
    return ';'.join(res)


def fill_operons(annotation: pd.DataFrame):
    operon = float('nan')

    indexes = []
    buffer = []
    
    for i, row in annotation.iterrows():
        op = float(annotation.at[i, 'operon'])
        
        if row['type'] == 'insertion_sequence':
            annotation.at[i, 'attributes'] = 'function=transposasa'
        
        if not math.isnan(op):
            if buffer:
                if op == operon:
                    indexes.extend(buffer)
                buffer = []

            operon = op
        else:
            if row['type'] == 'insertion_sequence':
                buffer.append((i, operon))
      
    for idx, operon in indexes:
        annotation.at[idx, 'attributes'] += ';' + _encode({'operon': int(operon)}) 
        
    return annotation         
        

def reindex_annotation(rec: SeqRecord, buildfile: pd.DataFrame, annotation: pd.DataFrame):
    seq = rec.seq
    out = []
    pos = 0
    inserts = []

    for _, (insert_idx, insertion, delta, score, strand) in buildfile.iterrows():
        prev = seq[pos:insert_idx]
        out.append(prev)
        out.append(insertion)
        
        _remap_genes(pos, insert_idx, delta, annotation)
        basis = insert_idx + delta

        inserts.append(
            mktrns(
                rec.id, 
                basis, 
                basis + len(insertion), 
                score,
                strand or '+',
                ""
            )
        )
        
        pos = insert_idx
        delta += len(insertion)

    _remap_genes(pos, len(seq), delta, annotation)
    out.append(seq[pos : len(seq)])
    
    for _, row_data in enumerate(inserts):
        annotation.loc[len(annotation)] = row_data

    annotation['start'] = annotation['nstart']
    annotation['end'] = annotation['nend']
    
    annotation = annotation.sort_values(by='start')
    annotation.reset_index()
    
    annotation = fill_operons(annotation)

    return Seq("".join(str(s) for s in out)), annotation


def _rebuild(args):
    buildfile = pd.read_csv(args.restore, 
                            names=['insert_index', 'insertion', 'delta', 'score', 'strand'])
    sequence = list(SeqIO.parse(args.fasta, 'fasta'))
    assert len(sequence) == 1
    sequence = sequence[0]
    
    operons = read_gff(args.gff)
    
    sq, reindexed = reindex_annotation(sequence, buildfile, operons)
    
    orig = list(SeqIO.parse(args.validation, 'fasta'))[0].seq
    assert sq == orig
    print('Rebuild successful!')
    
    write_gff(reindexed, args.output)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    subs = parser.add_subparsers(required=True)

    splice = subs.add_parser("cut")
    splice.add_argument("fasta", help="complete genome")
    splice.add_argument("gff", help="transposon annotation (with insertion_sequence)")
    splice.add_argument(
        "-o", "--output", help="out transposone-free .fna + .fna.restore.csv"
    )
    splice.set_defaults(func=_splice)

    rebuild = subs.add_parser("rebuild")
    rebuild.add_argument("fasta", help="transposone-free genome")
    rebuild.add_argument("restore", help="restore file (was generated by `splice`)")
    rebuild.add_argument("--gff", required=True, help="annotation to be changed")
    rebuild.add_argument("--validation", required=True, help="original fasta to validate")
    rebuild.add_argument("-o", "--output", help="new annotation and rebuilt genome")
    rebuild.set_defaults(func=_rebuild)

    args = parser.parse_args()

    return args.func(args)


if __name__ == "__main__":
    main()
