import unittest
import pickle
from Bio.Seq import reverse_complement
from genefeatures import gtf_tools as gt
from genefeatures.fasta_tools import extract_sequence
from genefeatures.sequence_index import SequenceIndex


class TestSequenceIndex(unittest.TestCase):

    def setUp(self):

        # forward transcript testing
        gtf = gt.parse_gtf("tests/data/test_hs_grch38.gtf")
        records = gtf.query(
            {"attributes": {"transcript_id": "ENST00000511072"}},
            return_records=True
        )
        self.itree_for = gt.records_to_interval_tree(records)
        start = self.itree_for.begin()
        end = self.itree_for.end()
        seq = extract_sequence(
            "tests/data/trunc_hs.grch38.dna.chr1.fa",
            1,
            start,
            end
        )
        self.seq_for = seq
        self.si_for = SequenceIndex(self.itree_for)

        # reverse transcript testing
        with open("tests/data/kras_gtfgff.pkl", "rb") as f:
            kras_gtf = pickle.load(f)
        kras = kras_gtf.query(
            {"attributes": {"transcript_id": "ENST00000311936"}},
            return_records=True
        )
        self.itree_rev = gt.records_to_interval_tree(kras)
        start = self.itree_rev.begin() - 25195145
        end = self.itree_rev.end() - 25195145
        seq = extract_sequence(
            "tests/data/kras_plus_10kb.fa",
            12,
            start,
            end
        )
        self.seq_rev = reverse_complement(seq)
        self.si_rev = SequenceIndex(self.itree_rev)

# Test genomic indices
    def test_genomic_index_forward(self):
        start = self.itree_for.begin()
        end = self.itree_for.end()
        seq_len = end - start
        self.assertEqual(self.si_for.genomic_idx[start], 0)
        self.assertEqual(self.si_for.genomic_idx[end], seq_len)

    def test_genomic_index_reverse(self):
        start = self.itree_rev.begin()
        end = self.itree_rev.end()
        seq_len = end - start
        self.assertEqual(self.si_rev.genomic_idx[start], 0)
        self.assertEqual(self.si_rev.genomic_idx[end], seq_len)

# Test forward transcript index

    def test_code_index_value_type(self):
        self.assertIsInstance(
            self.si_for.transcript_idx["1"],
            int
        )
        self.assertIsInstance(
            self.si_rev.transcript_idx["1"],
            int
        )

    def test_index_length_forward_matches(self):
        self.assertEqual(
            len(self.si_for.transcript_idx),
            len(self.si_for.genomic_idx)
        )

    def test_code_index_forward_start_codon(self):
        start_idx = slice(
            self.si_for.transcript_idx["1"],
            self.si_for.transcript_idx["4"]
        )
        self.assertEqual(
                self.seq_for[start_idx],
                "ATG"
        )

    def test_code_index_forward_utr_start_location(self):
        self.assertEqual(
            self.si_for.transcript_idx["-1"],
            self.si_for.transcript_idx["1"] - 1
        )

    def test_code_index_forward_stop_codon(self):
        stop_idx = slice(
            self.si_for.transcript_idx["*1"],
            self.si_for.transcript_idx["*4"]
        )
        self.assertEqual(
                self.seq_for[stop_idx],
                "TGA"
        )

    def test_code_index_forward_seq_index_first_element(self):
        self.assertEqual(
            self.si_for.transcript_idx["-92"],
            0
        )

    def test_code_index_forward_seq_index_last_element(self):
        self.assertEqual(
            self.si_for.transcript_idx["*656"],
            len(self.seq_for) - 1
        )

# Test reverse transcript index

    def test_index_length_reverse_matches(self):
        self.assertEqual(
            len(self.si_rev.transcript_idx),
            len(self.si_rev.genomic_idx)
        )

    def test_code_index_reverse_start_codon(self):
        start_idx = slice(
            self.si_rev.transcript_idx["1"],
            self.si_rev.transcript_idx["4"]
        )
        self.assertEqual(
                self.seq_rev[start_idx],
                "ATG"
        )

    def test_code_index_reverse_stop_codon(self):
        stop_idx = slice(
            self.si_rev.transcript_idx["*1"],
            self.si_rev.transcript_idx["*4"]
        )
        self.assertEqual(
                self.seq_rev[stop_idx],
                "TAA"
        )

    def test_code_index_reverse_utr_start_location(self):
        self.assertEqual(
            self.si_rev.transcript_idx["-1"] + 1,
            self.si_rev.transcript_idx["1"]
        )

    def test_code_index_reverse_seq_index_first_element(self):
        self.assertEqual(
            self.si_rev.transcript_idx["-5545"],
            0
        )

    def test_code_index_reverse_seq_index_last_element(self):
        self.assertEqual(
            self.si_rev.transcript_idx["*4552"],
            len(self.seq_rev) - 1
        )

    def test_log_change(self):
        self.si_for.log_change("insertion", 10, 15, 5)
        self.si_for.log_change("deletion", 20, 25, -5)
        expected_log = [("insertion", 10, 15, 5), ("deletion", 20, 25, -5)]
        self.assertEqual(self.si_for.change_log, expected_log)

    def test_update_index_insertion(self):
        original_idx = self.si_for.transcript_idx.copy()
        self.si_for.update_index(10, 15, 5)
        for key, value in original_idx.items():
            if value >= 15:
                self.assertEqual(self.si_for.transcript_idx[key], value + 5)
            else:
                self.assertEqual(self.si_for.transcript_idx[key], value)

    def test_update_index_deletion(self):
        original_idx = self.si_for.transcript_idx.copy()
        self.si_for.update_index(10, 15, -5)
        for key, value in original_idx.items():
            if value >= 15:
                self.assertEqual(self.si_for.transcript_idx[key], value - 5)
            else:
                self.assertEqual(self.si_for.transcript_idx[key], value)

    def test_update_index_no_change(self):
        original_idx = self.si_for.transcript_idx.copy()
        self.si_for.update_index(10, 15, 0)
        self.assertEqual(self.si_for.transcript_idx, original_idx)
