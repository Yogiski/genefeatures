from Bio.Seq import Seq
from intervaltree import IntervalTree, Interval
from .fasta_tools import extract_sequence
from typing import Type, TypeVar
from .gtf_tools import GtfGff

gtf = TypeVar("gtf", bound="GtfGff")


class SequenceTree:

    def __init__(
        self,
        seq_id=None,
        seqname=None,
        strand=None,
        interval: Interval | IntervalTree = None,
        sequence=None
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

        self.sequence = sequence
        self.seq_id = seq_id
        self.strand = strand
        self.seqname = seqname

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

    def read_sequence(self, fasta: str):

        self._check_seqnames()
        start = self.intervaltree.begin()
        end = self.intervaltree.end()
        seq = extract_sequence(
            fasta,
            self.seqname,
            start,
            end
        )
        self.sequence = Seq(seq)
        seq_idx = self._make_seq_index(start, end)
        self._seq_idx = seq_idx

    @staticmethod
    def _make_seq_index(start: int, end: int) -> dict:
        seq_idx = dict(zip(range(start, end+1), range(0, (end + 1) - start)))
        return seq_idx

    def _check_seqnames(self):

        for i in self.intervaltree:
            if self.seqname is None:
                self.seqname = i.data["seqname"]
            if i.data["seqname"] != self.seqname:
                raise ValueError(
                    f"""
                    seqnames inconsistent
                    found seqnames {self.seqname} and {i.data["seqname"]}
                    """
                )

    def get_coding_sequence(self):

        if not hasattr(self, "coding_sequence"):
            coding_seq = Seq("")
            for i in sorted(self.intervaltree):
                if i.data["feature"] == "CDS":
                    coding_seq += self._get_sub_sequence(
                        self._seq_idx,
                        self.sequence,
                        i.begin,
                        i.end + 1
                    )
            self.coding_seq = coding_seq
        else:
            coding_seq = self.coding_seq
        return coding_seq

    def translate(self):
        coding_seq = self.get_coding_sequence()
        aa_seq = coding_seq.translate()
        self.aa_seq = aa_seq
        return aa_seq

    @staticmethod
    def _get_sub_sequence(seq_idx, sequence, start, end):
        seq_start = seq_idx[start]
        seq_end = seq_idx[end]
        return sequence[seq_start:seq_end]

    def mutate(self):
        pass
