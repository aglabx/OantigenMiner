import argparse
from scripts.write_config import write_config


def setupdb(args):
    # TODO: move all db creations here (e.g., `bakta_loaddb.sh`, `Snakefile: rule operon_mapping`),
    # TODO: write subprocess to run command: `snakemake --use-conda --conda-create-envs-only -c 1`
    pass


def search(args):
    write_config()
    # TODO: write subprocess to run command: `snakemake -c {threads} --use-conda --report`
    pass


def clear_conda(args):
    # TODO: remove folder '.snakemake/conda'
    pass


def main():
    parser = argparse.ArgumentParser()
    subs = parser.add_subparsers(required=True)

    setupdb = subs.add_parser("setup", help='setup all databases, packages and '
                                            'conda environments')
    setupdb.add_argument(
        "-o", "--output", default='.', help="path to databases"
    )
    setupdb.set_defaults(func=setupdb)

    search = subs.add_parser("search", help='run WOOF to find O-antigen operons')
    search.add_argument("-f", "--fasta", required=True, help="input genome")
    search.add_argument("-t", "--threads", required=True, default=1, help="maximum number of threads")
    search.add_argument("-e", "--email", required=True, help="email, "
                                                             "required for call of online tools (e.g. OperonMapper)")
    search.add_argument("-o", "--output", required=True, help="path to output folder")
    search.add_argument("--dbs", required=True, help="path to databases (see `woof.py setup`)")
    search.add_argument("-g", "--gff", default=None, help="annotation of genome in gff format "
                                                          "(otherwise annotation with Bakta will be performed)")
    search.add_argument("--min-genes", default=3, help="minimal number of found O-antigen "
                                                       "formation related genes in operon to consider "
                                                       "the operon to be O-antigen operon, (default: 3)")
    search.add_argument("--max-hmm-eval", default=3e-20, help="maximal e-value of found hmm match "
                                                              "(default: 3e-20)")
    # TODO: add transcriptome, add reference gff3
    search.set_defaults(func=search)

    clearconda = subs.add_parser("clear-conda", help='remove all conda envs, created by WOOF')
    clearconda.set_defaults(func=clear_conda)

    args = parser.parse_args()

    return args.func(args)


if __name__ == "__main__":
    main()
