import unittest
import pickle
from Bio.Seq import Seq
from genefeatures import fasta_tools as ft
from genefeatures.sequence_handler import SequenceHandler
from genefeatures.sequence_tree import SequenceTree


class TestSequenceHandler(unittest.TestCase):

    def setUp(self):

        self.fasta = "tests/data/trunc_hs.grch38.dna.chr1.fa"
        # setup reverse strand sequence tree
        with open("tests/data/kras_gtfgff.pkl", "rb") as f:
            kras_gtf = pickle.load(f)

        kras = kras_gtf.query(
            {"attributes": {"transcript_id": "ENST00000311936"}}
        )
        rev = SequenceTree.from_gtf_gff(kras)
        start = rev.intervaltree.begin()
        end = rev.intervaltree.end()

        seq_start = start - 25195146
        seq_end = end - 25195146
        seq = ft.extract_sequence(
            "tests/data/kras_plus_10kb.fa",
            12,
            seq_start,
            seq_end
        )
        rev.set_full_seq(seq)
        # save coding and full sequences for mutation tests
        self.full = rev.get_full_seq()
        self.code = rev.get_coding_seq()
        rev.set_codon_index()
        rev.set_aa_seq()
        # init object to test
        self.sh = SequenceHandler(
            rev._seq_index, rev._coding_index, rev._codon_index, rev.strand
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
        mt_coding, mt_full = self.sh.dna_snv(
            self.code, self.full, ("34", "G", "T")
        )
        self.assertEqual(mt_coding[33:36], "TGT")
        self.assertEqual(mt_full[40104:40107], "ACA")
        with self.assertRaises(ValueError):
            self.sh.dna_snv(self.code, self.full, ("33", "G", "T"))

    def test_dna_point_deletion(self):
        mt_coding, mt_full = self.sh.dna_point_deletion(
            self.code, self.full, ("34", "")
        )
        self.assertEqual(mt_coding[33:36], "GTG")
        self.assertEqual(mt_full[40103:40106], "CAC")
        with self.assertRaises(ValueError):
            self.sh.dna_point_deletion(self.code, self.full, ("34", "C"))

    def test_dna_range_deletion(self):
        mt_coding, mt_full = self.sh.dna_range_deletion(
            self.code, self.full, ("34", "36", "")
        )
        self.assertEqual(mt_coding[33:36], "GGC")
        self.assertEqual(mt_full[40101:40104], "GCC")
        with self.assertRaises(ValueError):
            self.sh.dna_range_deletion(
                self.code, self.full, ("34", "36", "AAA")
            )

    def test_dna_insertion(self):
        mt_coding, mt_full = self.sh.dna_insertion(
            self.code, self.full, ("33", "34", "AAA")
        )
        self.assertEqual(mt_coding[33:36], "AAA")
        self.assertEqual(mt_coding[36:39], "GGT")
        self.assertEqual(mt_full[40107:40110], "TTT")
        self.assertEqual(mt_full[40104:40107], "ACC")

    def test_dna_duplication(self):
        mt_coding, mt_full = self.sh.dna_duplication(
            self.code, self.full, ("34", "36")
        )
        self.assertEqual(mt_coding[33:36], "GGT")
        self.assertEqual(mt_coding[36:39], "GGT")
        self.assertEqual(mt_full[40104:40107], "ACC")
        self.assertEqual(mt_full[40107:40110], "ACC")

    def test_dna_inv(self):
        mt_coding, mt_full = self.sh.dna_inversion(
            self.code, self.full, ("34", "36", "3")
        )
        self.assertEqual(mt_coding[33:36], "TGG")
        self.assertEqual(mt_coding[36:39], "GGC")
        self.assertEqual(mt_full[40104:40107], "CCA")
        self.assertEqual(mt_full[40101:40104], "GCC")

    def test_dna_indel(self):
        mt_coding, mt_full = self.sh.dna_indel(
            self.code,
            self.full,
            ("34", "36", "ACGT", None, None, None, None)
        )
        self.assertEqual(mt_coding[30:33], "GCT")
        self.assertEqual(mt_coding[33:37], "ACGT")
        self.assertEqual(mt_coding[37:40], "GGC")
        self.assertEqual(mt_full[40104:40108], "ACGT")
        self.assertEqual(mt_full[40101:40104], "GCC")
        mt_coding, mt_full = self.sh.dna_indel(
            self.code,
            self.full,
            (None, None, None, "34", "36", "GGT", "ACGT")
        )
        self.assertEqual(mt_coding[30:33], "GCT")
        self.assertEqual(mt_coding[33:37], "ACGT")
        self.assertEqual(mt_coding[37:40], "GGC")
        self.assertEqual(mt_full[40104:40108], "ACGT")
        self.assertEqual(mt_full[40101:40104], "GCC")
