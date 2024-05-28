import re
from typing import Tuple


class SequenceVariationParser:

    def __init__(self):

        self.seq_var_patterns = {
            "genomic": re.compile(r"g\.(.+)"),
            "cDNA": re.compile(r"c\.(.+)"),
            "mitochondrial": re.compile(r"m\.(.+)"),
            "RNA": re.compile(r"r\.(.+)"),
            "protein": re.compile(r"p\.(.+)")
        }
        self.dna_patterns = {
            "subs": re.compile(r"(\d+)([ACGT])>([ACGT])"),
            "indel": re.compile(
                r"(\d+)_(\d+)delins([ACGT]+)|"
                r"(\d+)_(\d+)del([ACGT]+)ins([ACGT]+)"
            ),
            "point_del": re.compile(r"(\d+)del(?!ins)([ACGT]*)"),
            "range_del": re.compile(r"(\d+)_(\d+)del(?!ins)([ACGT]*)"),
            "ins": re.compile(r"(\d+)_(\d+)ins([ACGT]+)"),
            "dup": re.compile(r"(\d+)_(\d+)dup([ACGT]*)"),
            "inv": re.compile(r"(\d+)_(\d+)inv(\d*)")
        }

    def match_variation_pattern(
        self,
        variation: str
    ) -> Tuple[str, str]:

        for var_type, pattern in self.seq_var_patterns.items():
            match = pattern.match(variation)
            if match:
                return var_type, match.group(1)
        raise ValueError(f"Unrecognized variation format: {variation}")

    def match_dna_change_pattern(self, change: str) -> Tuple[str, tuple]:

        for key, pattern in self.dna_patterns.items():
            match = pattern.match(change)
            if match:
                return key, match.groups()
        raise ValueError(f"Unrecognized dna change format: {change}")
