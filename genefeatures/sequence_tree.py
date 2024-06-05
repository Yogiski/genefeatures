import os
from collections import OrderedDict
from Bio.Seq import Seq, reverse_complement
from intervaltree import IntervalTree, Interval
from typing import Type, TypeVar
from .variation_parser import SequenceVariationParser
from .fasta_tools import extract_sequence
from .gtf_tools import GtfGff

gtf = TypeVar("gtf", bound="GtfGff")


class SequenceTree:

    def __init__(
        self,
        seq_id=None,
        seqname=None,
        strand=None,
        interval: Interval | IntervalTree = None,
        fasta=None
    ):
        if isinstance(interval, IntervalTree):
            self.intervaltree = interval
        elif isinstance(interval, Interval) or interval is None:
            self.intervaltree = IntervalTree(interval)
        else:
            raise TypeError(
                "interval argument must be types:"
                "Interval, IntervalTree, or None;"
                f"got {interval} of type {type(interval)}"
            )
        # init public attrs
        self.seq_id = seq_id
        self.strand = strand
        self.seqname = seqname
        self.fasta = fasta
        self.mutations = []
        # init private attrs
        self._sequence = None
        self._seq_index = None
        self._coding_seq = None
        self._codon_index = None
        self._aa_seq = None
        self._svp = SequenceVariationParser()

    # init methods
    @classmethod
    def from_gtf_gff(cls, records: list[dict] | Type[gtf], seq_id: str = None):

        if isinstance(records, GtfGff):
            records = records.export_records()
        elif isinstance(records, list):
            pass
        else:
            raise TypeError(
                "records must be type list or GtfGff;"
                f"got type {type(records)}"
            )
        try:
            intervaltree = IntervalTree()
            seqname = records[0]["seqname"]
            strand = records[0]["strand"]

            for r in records:
                # make sure seqname is consistent
                if r["seqname"] != seqname:
                    raise ValueError(
                        "seqname inconsistent;"
                        f"found seqnames: {seqname} and {r["seqname"]}"
                    )
                else:
                    inter = Interval(r.pop("start"), r.pop("end"), data=r)
                    intervaltree.add(inter)
        except KeyError as e:
            raise KeyError(
                f"record incomplete, field {e} is missing"
            )

        seqtree = cls(
            seq_id=seq_id,
            interval=intervaltree,
            seqname=seqname,
            strand=strand
        )
        return seqtree

    def add_interval(self, interval: Interval, sequence: str = None):

        if not isinstance(interval, Interval):
            raise ValueError("Interval must be an instance of Interval")
        if not (isinstance(sequence, str) or sequence is None):
            raise ValueError("Sequence must be a string or None")

        self.intervaltree.add(interval)

    # saftey checking functions
    def _check_strand(self):
        # if none assume forward transcription
        if self.strand in (1, "1", "+", ".", None, 0):
            self.strand = "+"
        elif self.strand in (-1, "-1", "-"):
            self.strand = "-"
        else:
            raise ValueError(
                "Strand must be + or 1 for positive or - or -1 for negative"
                f"got {self.strand}"
            )

    def _check_seqnames(self):

        for i in self.intervaltree:
            if self.seqname is None:
                self.seqname = i.data["seqname"]
            if i.data["seqname"] != self.seqname:
                raise ValueError(
                    "seqnames inconsistent found seqnames"
                    f"{self.seqname} and {i.data["seqname"]}"
                )

    def _check_fasta(self, fasta=None):

        if fasta is None:
            fasta = self.fasta
        if fasta is None:
            raise ValueError(
                f"object {self} has not attribute fasta; "
                "fasta must be set to read sequence"
            )
        fa_exists = os.path.isfile(fasta)
        fai_exists = os.path.isfile(f"{fasta}.fai")
        if not (fa_exists or fai_exists):
            raise FileNotFoundError(
                "fasta file and fasta index file must be present to "
                f"read sequence\nno such file: {fasta}"
            )

    # full NT sequence methods
    def read_full_seq(self, fasta: str = None, inplace=False):
        self._check_seqnames()
        if fasta is None:
            fasta = self.fasta
        self._check_fasta(fasta)
        start = self.intervaltree.begin()
        end = self.intervaltree.end()
        seq = extract_sequence(
            fasta,
            self.seqname,
            start,
            end + 1
        )
        seq = Seq(seq)
        if inplace:
            self.set_full_seq(seq)
        else:
            return seq

    def set_full_seq(self, sequence=None):
        if sequence is None:
            sequence = self.read_full_seq()
        self._sequence = sequence

    def get_full_seq(self):
        if self._sequence is None:
            self.set_full_seq()
        return self._sequence

    # seq index methods, for index NT seq wrt genomic coords
    @staticmethod
    def _init_seq_index(start: int, end: int) -> dict:
        seq_index = dict(zip(range(start, end+1), range(0, (end + 1) - start)))
        return seq_index

    def set_seq_index(self, start: int = None, end: int = None):
        if start is None or end is None:
            start = self.intervaltree.begin()
            end = self.intervaltree.end() + 1
        self._seq_index = self._init_seq_index(start, end)

    def get_seq_index(self):
        if self._seq_index is None:
            self.set_seq_index()
        return self._seq_index

    @staticmethod
    def _get_sub_seq_idx(seq_idx, start, end, strand):
        if strand == "-":
            start += 1
            end += 2
        else:
            end += 1
        seq_start = seq_idx[start]
        seq_end = seq_idx[end]
        return seq_start, seq_end

    # coding sequence methods
    def _init_coding_seq(self):

        coding_seq = Seq("")
        seq_index = self.get_seq_index()
        full_seq = self.get_full_seq()
        coding_index = []
        for i in sorted(self.intervaltree):
            if i.data["feature"] == "CDS":
                seq_start, seq_end = self._get_sub_seq_idx(
                    seq_index,
                    i.begin,
                    i.end,
                    self.strand
                )
                coding_seq += full_seq[seq_start:seq_end]
                coding_index += list(range(seq_start, seq_end))

        if self.strand == "-":
            coding_seq = coding_seq.reverse_complement()
            coding_index = list(reversed(coding_index))
        # coding index translates coding position to full seq position
        self._coding_index = dict(
            zip(range(0, len(coding_index)), coding_index)
        )
        return coding_seq

    def set_coding_seq(self, coding_seq=None):
        self._check_strand()
        if coding_seq is None:
            coding_seq = self._init_coding_seq()
        self._coding_seq = coding_seq

    def get_coding_seq(self):
        if self._coding_seq is None:
            self.set_coding_seq()
        return self._coding_seq

    # amino acid sequence methods
    @staticmethod
    def translate_coding_seq(coding_seq):
        return coding_seq.translate()

    def set_aa_seq(self, aa_seq=None):
        if aa_seq is None:
            aa_seq = self.get_coding_seq().translate()
        self._aa_seq = aa_seq

    def get_aa_seq(self):
        if self._aa_seq is None:
            self.set_aa_seq()
        return self._aa_seq

    @staticmethod
    def _init_codon_index(aa_seq, nt_seq):
        codons = [nt_seq[start:start+3] for start in range(0, len(nt_seq), 3)]
        codon_index = OrderedDict(zip(aa_seq, codons))
        return codon_index

    def set_codon_index(self):
        nt_seq = self.get_coding_seq()
        aa_seq = self.translate_coding_seq(nt_seq)
        codon_index = self._init_codon_index(aa_seq, nt_seq)
        self._codon_index = codon_index

    def get_codon_index(self):
        if self._codon_index is None:
            self.set_codon_index()
        return self._codon_index

    def _dna_change(self, change):

        change_type, groups = self._svp.match_dna_change_pattern(change)
        if change_type == "subs":
            mutated_seq = self._dna_snv(groups)

        elif change_type == "point_del":
            mutated_seq = self._dna_point_deletion(groups)

        elif change_type == "range_del":
            mutated_seq = self._dna_range_deletion(groups)

        elif change_type == "ins":
            mutated_seq = self._dna_insertion(groups)

        elif change_type == "dup":
            mutated_seq = self._dna_duplication(groups)

        elif change_type == "inv":
            mutated_seq = self._dna_inversion(groups)

        elif change_type == "indel":
            mutated_seq = self._dna_indel(groups)

        self.set_coding_seq(mutated_seq[0])
        self.set_full_seq(mutated_seq[1])

    def _mutate_sequence(
        self,
        sequence: Seq,
        start: int,
        end: int,
        ref: str,
        alt: str
    ) -> Seq:
        if sequence[start:end] != ref:
            raise ValueError(
                f"Reference base(s) {ref} do not match "
                f"at position {start+1}-{end+1} "
                f"found {sequence[start:end]}"
            )
        return sequence[:start] + alt + sequence[end:]

    def _get_mutated_sequences(
        self,
        pos: int,
        end: int,
        ref: str = "",
        alt: str = ""
    ) -> tuple:

        code = self.get_coding_seq()
        mutated_code = self._mutate_sequence(code, pos, end, ref, alt)

        full = self.get_full_seq()
        if self.strand == "-":
            ref = reverse_complement(ref)
            alt = reverse_complement(alt)
            full_pos = self._coding_index[end] + 1
            full_end = self._coding_index[pos] + 1
        else:
            full_pos = self._coding_index[pos]
            full_end = self._coding_index[end]

        mutated_full = self._mutate_sequence(
            full, full_pos, full_end, ref, alt
        )
        return mutated_code, mutated_full

    def _dna_snv(self, groups: tuple) -> tuple:
        pos, ref, alt = groups
        pos = int(pos) - 1
        end = pos + 1
        return self._get_mutated_sequences(pos, end, ref=ref, alt=alt)

    def _dna_point_deletion(self, groups: tuple) -> tuple:
        pos, ref = groups
        pos = int(pos) - 1
        if ref == "":
            ref = self.get_coding_seq()[pos]
        end = pos + 1
        return self._get_mutated_sequences(pos, end, ref=ref)

    def _dna_range_deletion(self, groups: tuple) -> tuple:
        start, end, ref = groups
        start, end = int(start) - 1, int(end)
        if ref == "":
            ref = self.get_coding_seq()[start:end]
        return self._get_mutated_sequences(start, end, ref=ref)

    def _dna_insertion(self, groups: tuple) -> tuple:
        start, end, alt = groups
        start, end = int(start), int(end)
        end = start
        ref = self.get_coding_seq()[start:end]
        return self._get_mutated_sequences(start, end, ref=ref, alt=alt)

    def _dna_duplication(self, groups: tuple) -> tuple:
        start, end = groups
        start, end = int(start) - 1, int(end)
        ref = self.get_coding_seq()[start:end]
        alt = ref + ref
        return self._get_mutated_sequences(start, end, ref=ref, alt=alt)

    def _dna_inversion(self, groups: tuple) -> tuple:
        start, end, length = groups
        start, end = int(start) - 1, int(end)
        ref = self.get_coding_seq()[start:end]
        alt = ref[::-1]
        return self._get_mutated_sequences(start, end, ref=ref, alt=alt)

    def _dna_indel(self, groups: tuple) -> tuple:

        if groups[3] is None:
            start, end, alt = groups[:3]
            start, end = int(start) - 1, int(end)
            ref = self.get_coding_seq()[start:end]
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
        return self._get_mutated_sequences(start, end, ref=ref, alt=alt)
