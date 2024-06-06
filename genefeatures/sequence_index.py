from intervaltree import Interval, IntervalTree
from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass
class SequenceIndex:

    interval_tree: IntervalTree
    strand: str = field(init=False)
    genomic_range: str = field(init=False)
    seq_index: Dict[int, int] = field(init=False, default_factory=dict)
    code_index: Dict[str, int] = field(init=False, default_factory=dict)

    def __post_init__(self):
        self.strand = self.get_strand()
        self.genomic_range = self.get_genomic_range()
        self.genomic_idx = self._init_genomic_index()
        self.transcript_idx = self._init_transcript_index()

    def get_strand(self) -> str:
        """Get the strand from the interval tree"""
        if self.interval_tree:
            for interval in self.interval_tree:
                return interval.data["strand"]
        return ''

    def get_genomic_range(self) -> str:
        """Get the genomic range from the interval tree"""
        if self.interval_tree:
            start_pos = self.interval_tree.begin()
            end_pos = self.interval_tree.end()
            return f"{start_pos}-{end_pos}"
        return ''

    def _init_genomic_index(self) -> Dict[int, int]:
        """make dict that converts genomic range to python index"""
        gen_range = range(
            self.interval_tree.begin(),
            self.interval_tree.end() + 1
            )
        seq_range = range(len(gen_range))
        return dict(zip(gen_range, seq_range))

    def _sort_intervals_get_positions(
        self,
        filtered: IntervalTree,
        is_reverse: bool
    ) -> Tuple[List[Interval], int, int]:
        if is_reverse:
            sorted_intervals = sorted(
                filtered,
                key=lambda x: x.end,
                reverse=True
            )
            start_pos = self.interval_tree.end()
            end_pos = self.interval_tree.begin() - 1
        else:
            sorted_intervals = sorted(filtered)
            start_pos = self.interval_tree.begin()
            end_pos = self.interval_tree.end() + 1

        return sorted_intervals, start_pos, end_pos

    @staticmethod
    def _get_interval_indeces(
        interval: Interval,
        is_reverse: bool
    ) -> Tuple[int, int, int]:
        if is_reverse:
            return interval.end, interval.begin - 1, -1
        else:
            return interval.begin, interval.end + 1, 1

    @staticmethod
    def _extend_trans_indices(
        trans_indices: List[str],
        prefix: str,
        start: int,
        end: int,
        step: int
    ) -> None:
        trans_indices.extend([f"{prefix}{i}" for i in range(start, end, step)])

    def _init_transcript_index(self) -> Dict[str, int]:

        trans_indices = []
        cds_idx = 1
        is_reverse = self.strand == "-"
        filtered = filter(
            lambda x: x.data["feature"] == "CDS", self.interval_tree
        )
        sorted_intervals, start, end = self._sort_intervals_get_positions(
            filtered, is_reverse
        )

        for i, interval in enumerate(sorted_intervals):
            indices = self._get_interval_indeces(interval, is_reverse)
            # make pre start codon index
            if i == 0:
                last_indices = indices[0] - indices[2]
                self._extend_trans_indices(
                    trans_indices, "-", abs(start - indices[0]), 0, -1
                )
            # make intron index
            if indices[0] - indices[2] != last_indices:
                self._extend_trans_indices(
                    trans_indices,
                    f"{cds_idx}+",
                    1,
                    abs(indices[0] - last_indices) + 1,
                    1
                )
            # make CDS index
            seq_len = abs(indices[0] - indices[1])
            self._extend_trans_indices(
                trans_indices, "", cds_idx, cds_idx + seq_len, 1
            )
            cds_idx += seq_len
            last_indices = indices[1]
        # stop codon - post stop utr index
        self._extend_trans_indices(
            trans_indices, "*", 1, abs(end - indices[1]) + 1, 1
        )

        seq_idx = self.genomic_idx.values()
        transcript_idx = dict(zip(trans_indices, sorted(seq_idx)))
        return transcript_idx

    def __getitem__(self, key: str) -> int:
        return self.transcript_idx[key]
