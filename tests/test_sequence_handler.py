import unittest
import pickle
from Bio.Seq import Seq, reverse_complement
from genefeatures import fasta_tools as ft
from genefeatures import gtf_tools as gt
from genefeatures.mutation_handler import MutationHandler
from genefeatures.sequence_index import SequenceIndex


class TestSequenceHandler(unittest.TestCase):

    def setUp(self):

        # setup reverse strand sequence tree
        with open("tests/data/kras_gtfgff.pkl", "rb") as f:
            kras_gtf = pickle.load(f)

        kras = kras_gtf.query(
            {"attributes": {"transcript_id": "ENST00000311936"}},
            return_records=True
        )
        si = SequenceIndex(gt.records_to_interval_tree(kras))
        self.si = si
        start = self.si.interval_tree.begin()
        end = self.si.interval_tree.end()

        seq_start = start - 25195145
        seq_end = end - 25195145
        seq = ft.extract_sequence(
            "tests/data/kras_plus_10kb.fa",
            12,
            seq_start,
            seq_end
        )
        self.seq = reverse_complement(seq)
        # init object to test
        self.sh = MutationHandler(
            si
        )

    def test_mutate_sequence_sub(self):
        seq = Seq("ACGT")
        m = self.sh.mutate_sequence(seq, 1, 2, "C", "X")
        self.assertEqual(m, "AXGT")
        m = self.sh.mutate_sequence(seq, 1, 3, "CG", "XX")
        self.assertEqual(m, "AXXT")
        self.assertIsInstance(m, Seq)

    def test_mutate_sequence_ins(self):
        seq = Seq("ACGT")
        m = self.sh.mutate_sequence(seq, 1, 1, "", "X")
        self.assertEqual(m, "AXCGT")
        m = self.sh.mutate_sequence(seq, 1, 1, "", "XX")
        self.assertEqual(m, "AXXCGT")

    def test_mutate_sequence_del(self):
        seq = Seq("ACGT")
        m = self.sh.mutate_sequence(seq, 1, 3, "CG", "")
        self.assertEqual(m, "AT")

    def test_mutate_sequence_indel(self):
        seq = Seq("ACGT")
        m = self.sh.mutate_sequence(seq, 1, 3, "CG", "XXX")
        self.assertEqual(m, "AXXXT")

    def test_dna_change_snv(self):
        mt_seq = self.sh.dna_snv(
            self.seq, ("34", "G", "T")
        )
        idx = self.si["34"]
        self.assertEqual(mt_seq[idx:idx+3], "TGT")
        with self.assertRaises(ValueError):
            self.sh.dna_snv(self.seq, ("33", "G", "T"))

    def test_dna_point_deletion(self):
        mt_seq = self.sh.dna_point_deletion(
            self.seq, ("34", "")
        )
        idx = self.si["34"]
        self.assertEqual(mt_seq[idx:idx+3], "GTG")
        with self.assertRaises(ValueError):
            self.sh.dna_point_deletion(self.seq, ("34", "C"))

    def test_dna_range_deletion(self):
        mt_seq = self.sh.dna_range_deletion(self.seq, ("34", "36", ""))
        idx = self.si["34"]
        self.assertEqual(mt_seq[idx:idx+3], "GGC")
        with self.assertRaises(ValueError):
            self.sh.dna_range_deletion(
                self.seq, ("34", "36", "AAA")
            )

    def test_dna_insertion(self):
        mt_seq = self.sh.dna_insertion(
            self.seq, ("33", "34", "AAA")
        )
        idx = self.si["34"]
        self.assertEqual(mt_seq[idx:idx+3], "AAA")
        self.assertEqual(mt_seq[idx+3:idx+6], "GGT")

    def test_dna_duplication(self):
        mt_seq = self.sh.dna_duplication(self.seq, ("34", "36"))
        idx = self.si["34"]
        self.assertEqual(mt_seq[idx:idx+3], "GGT")
        self.assertEqual(mt_seq[idx+3:idx+6], "GGT")

    def test_dna_inv(self):
        mt_seq = self.sh.dna_inversion(self.seq, ("34", "36", "3"))
        idx = self.si["34"]
        self.assertEqual(mt_seq[idx:idx+3], "TGG")
        self.assertEqual(mt_seq[idx+3:idx+6], "GGC")

    def test_dna_indel(self):
        mt_seq = self.sh.dna_indel(
            self.seq, ("34", "36", "ACGT", None, None, None, None)
        )
        idx = self.si["34"]
        self.assertEqual(mt_seq[idx-3:idx], "GCT")
        self.assertEqual(mt_seq[idx:idx+4], "ACGT")
        self.assertEqual(mt_seq[idx+4:idx+7], "GGC")

        mt_seq = self.sh.dna_indel(
            self.seq, (None, None, None, "34", "36", "GGT", "ACGT")
        )
        self.assertEqual(mt_seq[idx-3:idx], "GCT")
        self.assertEqual(mt_seq[idx:idx+4], "ACGT")
        self.assertEqual(mt_seq[idx+4:idx+7], "GGC")
