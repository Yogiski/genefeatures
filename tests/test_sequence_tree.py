import unittest
from Bio.Seq import Seq
import pickle
from genefeatures.sequence_tree import SequenceTree
from genefeatures import gtf_tools as gt
from genefeatures import fasta_tools as ft


class TestSequenceTree(unittest.TestCase):

    def setUp(self):
        gtf = gt.parse_gtf("tests/data/test_hs_grch38.gtf")
        self.gtf = gtf.query(
            {
                "attributes": {"transcript_id": "ENST00000511072"}
            }
        )
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
        rev.set_seq_index(start=start, end=end)
        self.rev = rev

    def test_from_gtf(self):
        recs = self.gtf.export_records()
        st = SequenceTree.from_gtf_gff(recs)
        self.assertIsInstance(st, SequenceTree)
        st = SequenceTree.from_gtf_gff(self.gtf)
        self.assertIsInstance(st, SequenceTree)

    def test_from_gtf_missing_start(self):
        with self.assertRaises(KeyError):
            recs = self.gtf.export_records()
            recs[0].pop("start")
            SequenceTree.from_gtf_gff(recs)

    def test_from_gtf_missing_seqname(self):
        with self.assertRaises(KeyError):
            recs = self.gtf.export_records()
            recs[0].pop("seqname")
            SequenceTree.from_gtf_gff(recs)

    def test_from_gtf_missing(self):
        with self.assertRaises(KeyError):
            recs = self.gtf.export_records()
            recs[0].pop("strand")
            SequenceTree.from_gtf_gff(recs)

    def test_from_gtf_multiple_seqnames(self):
        with self.assertRaises(ValueError):
            recs = self.gtf.export_records()
            recs[0]["seqname"] = "6"
            SequenceTree.from_gtf_gff(recs)

    def test_init_seq_index(self):
        st = SequenceTree()
        seq_idx = st._init_seq_index(500, 510)
        self.assertEqual(list(range(0, 11)), list(seq_idx.values()))

    def test_get_sequence(self):
        st = SequenceTree.from_gtf_gff(self.gtf)
        seq = st.read_full_seq(self.fasta)
        self.assertIsInstance(seq, Seq)
        self.assertIn("A" or "C" or "G" or "T", seq)

    def test_get_coding_seq(self):
        st = SequenceTree.from_gtf_gff(self.gtf)
        st.fasta = self.fasta
        coding = st.get_coding_seq()
        self.assertTrue(coding.startswith("ATG"))
        self.assertEqual(len(coding) % 3, 0)

    def test_forward_translate(self):
        st = SequenceTree.from_gtf_gff(self.gtf)
        st.fasta = self.fasta
        aa = st.get_aa_seq()
        self.assertTrue(aa.startswith("M"))
        self.assertNotIn("*", aa)

    def test_reverse_translate(self):
        coding = self.rev.get_coding_seq()
        self.assertTrue(coding.startswith("ATG"))
        self.assertEqual(len(coding) % 3, 0)
        aa_seq = self.rev.get_aa_seq()
        self.assertTrue(aa_seq.startswith("M"))

    def test_codon_index(self):
        coding = self.rev.get_coding_seq()
        aa_seq = self.rev.get_aa_seq()
        codon_index = self.rev._init_codon_index(aa_seq, coding)
        self.assertEqual(codon_index[aa_seq[0]], "ATG")
        nt_codons = [str(i) for i in codon_index.values()]
        self.assertEqual(len("".join(nt_codons)) % 3, 0)

    def test_match_dna_change_pattern(self):
        # substitution
        change_type, groups = self.rev._match_dna_change_pattern("34G>T")
        self.assertEqual(change_type, "subs")
        # single NT deletion
        change_type, groups = self.rev._match_dna_change_pattern("34del")
        self.assertEqual(change_type, "point_del")
        self.assertEqual(groups[1], "")
        change_type, groups = self.rev._match_dna_change_pattern("34delG")
        self.assertEqual(change_type, "point_del")
        self.assertEqual(groups[1], "G")
        # multiple NT deletion
        change_type, groups = self.rev._match_dna_change_pattern("34_36del")
        self.assertEqual(change_type, "range_del")
        self.assertEqual(groups[2], "")
        change_type, groups = self.rev._match_dna_change_pattern("34_36delGGT")
        self.assertEqual(change_type, "range_del")
        self.assertEqual(groups[2], "GGT")
        # insertion
        change_type, groups = self.rev._match_dna_change_pattern("33_34insAAA")
        self.assertEqual(change_type, "ins")
        # duplication
        change_type, groups = self.rev._match_dna_change_pattern("34_36dupGGT")
        self.assertEqual(change_type, "dup")
        change_type, groups = self.rev._match_dna_change_pattern("34_36dup")
        self.assertEqual(change_type, "dup")
        change_type, groups = self.rev._match_dna_change_pattern("34_36inv")
        self.assertEqual(change_type, "inv")
        self.assertEqual(len(groups), 3)
        change_type, groups = self.rev._match_dna_change_pattern("34_36inv3")
        self.assertEqual(change_type, "inv")
        change_type, groups = self.rev._match_dna_change_pattern(
            "34_36delinsTAA"
        )
        self.assertEqual(change_type, "indel")
        self.assertIsNone(groups[3])
        self.assertIsNone(groups[4])
        self.assertIsNone(groups[5])
        self.assertIsNone(groups[6])
        change_type, groups = self.rev._match_dna_change_pattern(
            "34_36delGGTinsTAA"
        )
        self.assertEqual(change_type, "indel")
        self.assertIsNone(groups[0])
        self.assertIsNone(groups[1])
        self.assertIsNone(groups[2])

    def test_dna_change_snv(self):
        self.rev.get_coding_seq()
        mt_coding, mt_full = self.rev._dna_snv(("34", "G", "T"))
        self.assertEqual(mt_coding[33:36], "TGT")
        self.assertEqual(mt_full[40104:40107], "ACA")
        with self.assertRaises(ValueError):
            self.rev._dna_snv(("33", "G", "T"))

    def test_dna_point_deletion(self):
        self.rev.get_coding_seq()
        mt_coding, mt_full = self.rev._dna_point_deletion(("34", ""))
        self.assertEqual(mt_coding[33:36], "GTG")
        self.assertEqual(mt_full[40103:40106], "CAC")
        with self.assertRaises(ValueError):
            self.rev._dna_point_deletion(("34", "C"))

    def test_dna_range_deletion(self):
        self.rev.get_coding_seq()
        mt_coding, mt_full = self.rev._dna_range_deletion(("34", "36", ""))
        self.assertEqual(mt_coding[33:36], "GGC")
        self.assertEqual(mt_full[40101:40104], "GCC")
        with self.assertRaises(ValueError):
            self.rev._dna_range_deletion(("34", "36", "AAA"))

    def test_dna_insertion(self):
        self.rev.get_coding_seq()
        mt_coding, mt_full = self.rev._dna_insertion(("33", "34", "AAA"))
        self.assertEqual(mt_coding[33:36], "AAA")
        self.assertEqual(mt_coding[36:39], "GGT")
        self.assertEqual(mt_full[40107:40110], "TTT")
        self.assertEqual(mt_full[40104:40107], "ACC")

    def test_dna_duplication(self):
        self.rev.get_coding_seq()
        mt_coding, mt_full = self.rev._dna_duplication(("34", "36"))
        self.assertEqual(mt_coding[33:36], "GGT")
        self.assertEqual(mt_coding[36:39], "GGT")
        self.assertEqual(mt_full[40104:40107], "ACC")
        self.assertEqual(mt_full[40107:40110], "ACC")

    def test_dna_inv(self):
        self.rev.get_coding_seq()
        mt_coding, mt_full = self.rev._dna_inversion(("34", "36", "3"))
        self.assertEqual(mt_coding[33:36], "TGG")
        self.assertEqual(mt_coding[36:39], "GGC")
        self.assertEqual(mt_full[40104:40107], "CCA")
        self.assertEqual(mt_full[40101:40104], "GCC")

    def test_dna_indel(self):
        self.rev.get_coding_seq()
        mt_coding, mt_full = self.rev._dna_indel(
            ("34", "36", "ACGT", None, None, None, None)
        )
        self.assertEqual(mt_coding[30:33], "GCT")
        self.assertEqual(mt_coding[33:37], "ACGT")
        self.assertEqual(mt_coding[37:40], "GGC")
        self.assertEqual(mt_full[40104:40108], "ACGT")
        self.assertEqual(mt_full[40101:40104], "GCC")
        mt_coding, mt_full = self.rev._dna_indel(
            (None, None, None, "34", "36", "GGT", "ACGT")
        )
        self.assertEqual(mt_coding[30:33], "GCT")
        self.assertEqual(mt_coding[33:37], "ACGT")
        self.assertEqual(mt_coding[37:40], "GGC")
        self.assertEqual(mt_full[40104:40108], "ACGT")
        self.assertEqual(mt_full[40101:40104], "GCC")
