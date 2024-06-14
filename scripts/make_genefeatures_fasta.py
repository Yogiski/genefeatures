import argparse
from genefeatures import gtf_tools as gt
# from genefeatures import fasta_tools as ft


def main(fasta: str, gtf: str, model: str, mutations: str) -> None:
    gtf = gt.parse_gtf(gtf)
    raise NotImplementedError


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
        "genes",
        help=(
            "3 column csv or txt file specifying genes, transcript, and "
            "mutations needed for genefeature dataset"
        )
    )
    args = parser.parse_args()

    main(*args)
