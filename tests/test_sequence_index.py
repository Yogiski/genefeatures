import unittest
import pickle
from Bio.Seq import reverse_complement
from genefeatures import gtf_tools as gt
from genefeatures.fasta_tools import extract_sequence
from genefeatures.sequence_tree import SequenceTree
from genefeatures.sequence_index import SequenceIndex


class TestSequenceIndex(unittest.TestCase):

    def setUp(self):

        gtf = gt.parse_gtf("tests/data/test_hs_grch38.gtf")
        records = gtf.query(
            {"attributes": {"transcript_id": "ENST00000511072"}}
        )
        self.itree_for = SequenceTree.from_gtf_gff(records).intervaltree
        start = self.itree_for.begin()
        end = self.itree_for.end()
        seq = extract_sequence(
            "tests/data/trunc_hs.grch38.dna.chr1.fa",
            1,
            start,
            end
        )
        print(seq[92:95])
        print(seq[-656:-653])
        self.seq_for = seq
        self.si_for = SequenceIndex(self.itree_for)

        with open("tests/data/kras_gtfgff.pkl", "rb") as f:
            kras_gtf = pickle.load(f)
        kras = kras_gtf.query(
            {"attributes": {"transcript_id": "ENST00000311936"}}
        )
        self.itree_rev = SequenceTree.from_gtf_gff(kras).intervaltree
        start = self.itree_rev.begin() - 25195145
        end = self.itree_rev.end() - 25195145
        seq = extract_sequence(
            "tests/data/kras_plus_10kb.fa",
            12,
            start,
            end
        )
        self.seq_rev = reverse_complement(seq)
        print(self.seq_rev[5545:5548])
        print(self.seq_rev[41132:41135])

        self.si_rev = SequenceIndex(self.itree_rev)

    def test_genomic_index_forward(self):
        start = self.itree_for.begin()
        end = self.itree_for.end()
        seq_len = end - start
        self.assertEqual(self.si_for.genomic_idx[start], 0)
        self.assertEqual(self.si_for.genomic_idx[end], seq_len)

    def test_code_index_forward(self):
        self.assertIsInstance(
            self.si_for.transcript_idx["1"],
            int
        )
        self.assertEqual(
            self.si_for.transcript_idx["-1"],
            self.si_for.transcript_idx["1"] - 1
        )
        self.assertEqual(
            len(self.si_for.transcript_idx),
            len(self.si_for.genomic_idx)
        )
        self.assertEqual(
            self.si_for.transcript_idx["-92"],
            0
        )
        self.assertEqual(
            self.si_for.transcript_idx["1"],
            92
        )
        self.assertEqual(
            self.si_for.transcript_idx["*656"],
            len(self.seq_for) - 1
        )
        start_idx = slice(
            self.si_for.transcript_idx["1"],
            self.si_for.transcript_idx["4"]
        )
        self.assertEqual(
                self.seq_for[start_idx],
                "ATG"
        )
        stop_idx = slice(
            self.si_for.transcript_idx["*1"],
            self.si_for.transcript_idx["*4"]
        )
        print(stop_idx)
        self.assertEqual(
                self.seq_for[stop_idx],
                "TGA"
        )

    def test_genomic_index_reverse(self):
        start = self.itree_rev.begin()
        end = self.itree_rev.end()
        seq_len = end - start
        self.assertEqual(self.si_rev.genomic_idx[start], 0)
        self.assertEqual(self.si_rev.genomic_idx[end], seq_len)

    def test_code_index_reverse(self):
        self.assertIsInstance(
            self.si_rev.transcript_idx["1"],
            int
        )
        self.assertEqual(
            self.si_rev.transcript_idx["-1"] + 1,
            self.si_rev.transcript_idx["1"]
        )
        self.assertEqual(
            len(self.si_rev.transcript_idx),
            len(self.si_rev.genomic_idx)
        )
        self.assertEqual(
            self.si_rev.transcript_idx["-5545"],
            0
        )
        self.assertEqual(
            self.si_rev.transcript_idx["1"],
            5545
        )
        self.assertEqual(
            self.si_rev.transcript_idx["*4552"],
            len(self.seq_rev) - 1
        )
        self.assertEqual(
            self.seq_rev[self.si_rev.transcript_idx["*4549"]:],
            "TGAC"
        )
        self.assertEqual(
            self.seq_rev[
                self.si_rev.transcript_idx["-5545"]:
                self.si_rev.transcript_idx["-5541"]
            ],
            "CTAG"
        )
        start_idx = slice(
            self.si_rev.transcript_idx["1"],
            self.si_rev.transcript_idx["4"]
        )
        self.assertEqual(
                self.seq_rev[start_idx],
                "ATG"
        )
        stop_idx = slice(
            self.si_rev.transcript_idx["*1"],
            self.si_rev.transcript_idx["*4"]
        )
        self.assertEqual(
                self.seq_rev[stop_idx],
                "TAA"
        )
