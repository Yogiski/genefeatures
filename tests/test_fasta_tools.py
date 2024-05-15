import unittest
from genefeatures import fasta_tools as ft

class TestFastaTools(unittest.TestCase):

    def setUp(self):
        self.fasta = "tests/data/trunc_hs.grch38.dna.chr1.fa"

    def test_fast_extract(self):
        seqname = 1
        start = 3069260
        stop = 3069262
        start_codon = ft.extract_sequence(self.fasta, seqname, start, stop)
        self.assertEqual(start_codon, "ATG")

        tstart = 3069211
        tstop = 3434342
        transcript = ft.extract_sequence(self.fasta, seqname, tstart, tstop)
        nucleotides = set(transcript)
        self.assertFalse("U" in nucleotides)
        self.assertIn("G", nucleotides)
        self.assertIn("C", nucleotides)
        self.assertIn("A", nucleotides)
        self.assertIn("T", nucleotides)

if __name__ == '__main__':
    unittest.main()