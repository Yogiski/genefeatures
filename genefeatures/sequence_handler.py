from Bio.Seq import Seq, reverse_complement


class SequenceHandler:

    def __init__(
        self,
        seq_index: dict = None,
        coding_index: dict = None,
        codon_index: dict = None,
        strand: str = None
    ):
        self.seq_index = seq_index
        self.coding_index = coding_index
        self.codon_index = codon_index
        self.strand = strand

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

    def get_mutated_sequences(
        self,
        coding: Seq,
        full: Seq,
        pos: int,
        end: int,
        ref: str = "",
        alt: str = ""
    ) -> tuple:

        mutated_code = self.mutate_sequence(coding, pos, end, ref, alt)

        if self.strand == "-":
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)
            full_pos = self.coding_index[end] + 1
            full_end = self.coding_index[pos] + 1
        else:
            full_pos = self.coding_index[pos]
            full_end = self.coding_index[end]

        mutated_full = self.mutate_sequence(
            full, full_pos, full_end, ref, alt
        )
        return mutated_code, mutated_full

    def dna_snv(self, coding: Seq, full: Seq, groups: tuple) -> tuple:
        pos, ref, alt = groups
        pos = int(pos) - 1
        end = pos + 1
        mutated = self.get_mutated_sequences(
            coding, full, pos, end, ref=ref, alt=alt
        )
        return mutated

    def dna_point_deletion(
        self, coding: Seq, full: Seq, groups: tuple
    ) -> tuple:
        pos, ref = groups
        pos = int(pos) - 1
        if ref == "":
            ref = coding[pos]
        end = pos + 1
        mutated = self.get_mutated_sequences(
            coding, full, pos, end, ref=ref
        )
        return mutated

    def dna_range_deletion(
        self, coding: Seq, full: Seq, groups: tuple
    ) -> tuple:
        start, end, ref = groups
        start, end = int(start) - 1, int(end)
        if ref == "":
            ref = coding[start:end]

        mutated = self.get_mutated_sequences(
            coding, full, start, end, ref=ref
        )
        return mutated

    def dna_insertion(self, coding: Seq, full: Seq, groups: tuple) -> tuple:
        start, end, alt = groups
        start, end = int(start), int(end)
        end = start
        ref = coding[start:end]
        mutated = self.get_mutated_sequences(
            coding, full, start, end, ref=ref, alt=alt
        )
        return mutated

    def dna_duplication(self, coding: Seq, full: Seq, groups: tuple) -> tuple:
        start, end = groups
        start, end = int(start) - 1, int(end)
        ref = coding[start:end]
        alt = ref + ref
        mutated = self.get_mutated_sequences(
            coding, full, start, end, ref=ref, alt=alt
        )
        return mutated

    def dna_inversion(self, coding: Seq, full: Seq, groups: tuple) -> tuple:
        start, end, length = groups
        start, end = int(start) - 1, int(end)
        ref = coding[start:end]
        alt = ref[::-1]
        mutated = self.get_mutated_sequences(
            coding, full, start, end, ref=ref, alt=alt
        )
        return mutated

    def dna_indel(self, coding: Seq, full: Seq, groups: tuple) -> tuple:

        if groups[3] is None:
            start, end, alt = groups[:3]
            start, end = int(start) - 1, int(end)
            ref = coding[start:end]
        elif groups[0] is None:
            start, end, ref, alt = groups[3:]
            start, end = int(start) - 1, int(end)
        else:
            raise ValueError(
                "indel groups formatted incorrectly. "
                "must be tuple with len of 6 "
                "either first 3 entries or last 4 entries must be None "
                f"got {groups}"
            )
        mutated = self.get_mutated_sequences(
            coding, full, start, end, ref=ref, alt=alt
        )
        return mutated
