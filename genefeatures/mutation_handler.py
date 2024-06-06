from .sequence_index import SequenceIndex
from Bio.Seq import Seq
from typing import Tuple


class MutationHandler:

    def __init__(
        self,
        seq_index: SequenceIndex
    ):
        self.seq_index = seq_index

    # mutation methods
    @staticmethod
    def mutate_sequence(
        sequence: Seq, start: int, end: int, ref: str, alt: str
    ) -> Seq:
        if sequence[start:end] != ref:
            raise ValueError(
                f"Reference base(s) {ref} do not match "
                f"at position {start+1}-{end+1} "
                f"found {sequence[start:end]}"
            )
        return sequence[:start] + alt + sequence[end:]

    def dna_snv(self, seq: Seq, groups: Tuple[str, str, str]) -> Seq:
        pos, ref, alt = groups
        pos = self.seq_index[pos]
        end = pos + 1
        return self.mutate_sequence(seq, pos, end, ref=ref, alt=alt)

    def dna_point_deletion(self, seq: Seq, groups: Tuple[str, str]) -> Seq:
        pos, ref = groups
        pos = self.seq_index[pos] - 1
        if ref == "":
            ref = seq[pos]
        end = pos + 1
        return self.mutate_sequence(seq, pos, end, ref, "")

    def dna_range_deletion(
        self, seq: Seq, groups: Tuple[str, str, str]
    ) -> Seq:

        start, end, ref = groups
        start, end = self.seq_index[start] - 1, self.seq_index[end]
        if ref == "":
            ref = seq[start:end]
        return self.mutate_sequence(seq, start, end, ref, "")

    def dna_insertion(self, seq: Seq, groups: Tuple[str, str, str]) -> Seq:
        start, end, alt = groups
        start, end = self.seq_index[start] + 1, self.seq_index[end] + 1
        end = start
        ref = seq[start:end]
        return self.mutate_sequence(seq, start, end, ref, alt)

    def dna_duplication(self, seq: Seq, groups: Tuple[str, str]) -> Seq:
        start, end = groups
        start, end = self.seq_index[start], self.seq_index[end] + 1
        ref = seq[start:end]
        alt = ref + ref
        return self.mutate_sequence(seq, start, end, ref, alt)

    def dna_inversion(self, seq: Seq, groups: Tuple[str, str, str]) -> tuple:
        start, end, length = groups
        start, end = self.seq_index[start], self.seq_index[end]+1
        ref = seq[start:end]
        alt = ref[::-1]
        return self.mutate_sequence(seq, start, end, ref, alt)

    def dna_indel(
        self, seq: Seq, groups: Tuple[str, str, str, str, str, str]
    ) -> Seq:

        if groups[3] is None:
            start, end, alt = groups[:3]
            start, end = self.seq_index[start], self.seq_index[end] + 1
            ref = seq[start:end]
        elif groups[0] is None:
            start, end, ref, alt = groups[3:]
            start, end = self.seq_index[start], self.seq_index[end] + 1
        else:
            raise ValueError(
                "indel groups formatted incorrectly. "
                "must be tuple with len of 6 "
                "either first 3 entries or last 4 entries must be None "
                f"got {groups}"
            )
        return self.mutate_sequence(seq, start, end, ref, alt)
