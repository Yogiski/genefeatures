import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from genefeatures import gtf_tools as gt
# Add the parent directory to the Python path


def construct_queries(mutations: list) -> list:
    queries = []
    for gene, trans, _ in mutations:
        if trans:
            queries.append(
                {"attributes": {"gene_id": gene, "transcript_id": trans}}
            )
        else:
            queries.append(
                {"attributes": {"gene_id": gene, "tag": "MANE_Select"}}
            )
    return queries


def read_mutations(mutations: str) -> list:
    with open(mutations, "r") as f:
        res = map(lambda x: tuple(x.rstrip("\n").split(",")), f)
        gene_trnscrpt_mut = list(res)
    return gene_trnscrpt_mut


def main(fasta: str, gtf: str, model: str, mutations: str) -> None:
    muts = read_mutations(mutations)
    queries = construct_queries(muts)
    print(f"Parsing GTF: {gtf}")
    gtf = gt.parse_gtf(gtf)
    print("Done!")
    print("Submitting Queries")
    for q in queries:
        gtf.query(q, return_records=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='make_genefeatures_fasta',
        description=(
            'generate transcript sequence fasta'
            'for a given cell model and it\'s mutations'
        )
    )
    parser.add_argument(
        'fasta',
        help="fasta file for genome sequence"
    )
    parser.add_argument(
        'gtf',
        help="gtf file for sequence element annotations"
    )
    parser.add_argument(
        'model',
        help="cell model id, used for naming conventions"
    )
    parser.add_argument(
        "mutations",
        help=(
            "3 column csv or txt file specifying genes, transcript, and "
            "mutations needed for genefeature dataset"
        )
    )
    args = parser.parse_args()
    main(args.fasta, args.gtf, args.model, args.mutations)
