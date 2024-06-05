from intervaltree import Interval, IntervalTree
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
# from functools import lru_cache


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
        start_pos = self.interval_tree.begin()
        end_pos = self.interval_tree.end() + 1
        gen_range = range(start_pos, end_pos)
        seq_range = range(end_pos - start_pos)
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
    def _handle_pre_start_utr(
        indices: Tuple[int, int, int],
        start_pos: int,
        trans_indices: List[str]
    ) -> None:
        seq_len = abs(start_pos - indices[0])
        trans_indices.extend(
            [f"-{i}" for i in range(seq_len, 0, -1)]
        )

    @staticmethod
    def _handle_intron(
        indices: Tuple[int, int, int],
        last_indices: int,
        cds_idx: int,
        trans_indices: List[str]
    ) -> None:
        seq_len = abs(indices[0] - last_indices)
        trans_indices.extend(
            [f"{cds_idx}+{i}" for i in range(1, seq_len + 1)]
        )

    @staticmethod
    def _handle_cds(
        indices: Tuple[int, int, int],
        cds_idx: int,
        trans_indices: List[str]
    ) -> int:
        seq_len = abs((indices[0] - indices[1]))
        trans_indices.extend(
            [f"{i}" for i in range(cds_idx, cds_idx + seq_len)]
        )
        return cds_idx + seq_len

    @staticmethod
    def _handle_post_stop_utr(
        indices: Tuple[int, int, int],
        end_pos: int,
        trans_indices: List[str]
    ) -> None:
        seq_len = abs(end_pos - indices[1])
        trans_indices.extend(
            [f"*{i}" for i in range(1, seq_len + 1)]
        )

    def _init_transcript_index(self) -> Dict[str, int]:
        trans_indices = []
        pre_start_codon = True
        cds_idx = 1
        is_reverse = self.strand == "-"
        filtered = filter(
            lambda x: x.data["feature"] == "CDS", self.interval_tree
        )
        sorted_intervals, start, end = self._sort_intervals_get_positions(
            filtered, is_reverse
        )
        for i in sorted_intervals:
            print(i.data["feature"], i.begin, i.end)

        for interval in sorted_intervals:
            indices = self._get_interval_indeces(interval, is_reverse)
            if pre_start_codon:
                last_indices = indices[0] - indices[2]
                pre_start_codon = False
                self._handle_pre_start_utr(indices, start, trans_indices)
            if indices[0] - indices[2] != last_indices:
                self._handle_intron(
                    indices, last_indices, cds_idx, trans_indices
                )
            cds_idx = self._handle_cds(indices, cds_idx, trans_indices)
            last_indices = indices[1]

        self._handle_post_stop_utr(indices, end, trans_indices)
        seq_idx = self.genomic_idx.values()
        transcript_idx = dict(zip(trans_indices, seq_idx))
        return transcript_idx
