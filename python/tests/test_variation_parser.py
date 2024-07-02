import unittest
from genefeatures.variation_parser import SequenceVariationParser


class TestSequenceVariationParser(unittest.TestCase):

    def setUp(self):
        self.svp = SequenceVariationParser()

    def test_match_variation_geno(self):
        var_type, groups = self.svp.match_variation_pattern("g.76A>T")
        self.assertEqual(var_type, "genomic")
        self.assertEqual(groups, "76A>T")

    def test_match_variation_cdna(self):
        var_type, groups = self.svp.match_variation_pattern("c.76A>T")
        self.assertEqual(var_type, "cDNA")
        self.assertEqual(groups, "76A>T")

    def test_match_variation_mito(self):
        var_type, groups = self.svp.match_variation_pattern("m.76A>T")
        self.assertEqual(var_type, "mitochondrial")
        self.assertEqual(groups, "76A>T")

    def test_match_variation_rna(self):
        var_type, groups = self.svp.match_variation_pattern("r.76a>u")
        self.assertEqual(var_type, "RNA")
        self.assertEqual(groups, "76a>u")

    def test_match_variation_protein(self):
        var_type, group = self.svp.match_variation_pattern("p.G12C")
        self.assertEqual(var_type, "protein")
        self.assertEqual(group, "G12C")

    def test_match_variation_pattern_raises(self):
        with self.assertRaises(ValueError):
            self.svp.match_variation_pattern("K76A")
        with self.assertRaises(ValueError):
            self.svp.match_variation_pattern("KRAS.K76A")
        with self.assertRaises(ValueError):
            self.svp.match_variation_pattern("KRAS.p.K76A")

    def test_match_dna_change_pattern(self):
        # substitution
        change_type, groups = self.svp.match_dna_change_pattern("34G>T")
        self.assertEqual(change_type, "subs")
        # single NT deletion
        change_type, groups = self.svp.match_dna_change_pattern("34del")
        self.assertEqual(change_type, "point_del")
        self.assertEqual(groups[1], "")
        change_type, groups = self.svp.match_dna_change_pattern("34delG")
        self.assertEqual(change_type, "point_del")
        self.assertEqual(groups[1], "G")
        # multiple NT deletion
        change_type, groups = self.svp.match_dna_change_pattern("34_36del")
        self.assertEqual(change_type, "range_del")
        self.assertEqual(groups[2], "")
        change_type, groups = self.svp.match_dna_change_pattern("34_36delGGT")
        self.assertEqual(change_type, "range_del")
        self.assertEqual(groups[2], "GGT")
        # insertion
        change_type, groups = self.svp.match_dna_change_pattern("33_34insAAA")
        self.assertEqual(change_type, "ins")
        # duplication
        change_type, groups = self.svp.match_dna_change_pattern("34_36dupGGT")
        self.assertEqual(change_type, "dup")
        change_type, groups = self.svp.match_dna_change_pattern("34_36dup")
        self.assertEqual(change_type, "dup")
        change_type, groups = self.svp.match_dna_change_pattern("34_36inv")
        self.assertEqual(change_type, "inv")
        self.assertEqual(len(groups), 3)
        change_type, groups = self.svp.match_dna_change_pattern("34_36inv3")
        self.assertEqual(change_type, "inv")
        change_type, groups = self.svp.match_dna_change_pattern(
            "34_36delinsTAA"
        )
        self.assertEqual(change_type, "indel")
        self.assertIsNone(groups[3])
        self.assertIsNone(groups[4])
        self.assertIsNone(groups[5])
        self.assertIsNone(groups[6])
        change_type, groups = self.svp.match_dna_change_pattern(
            "34_36delGGTinsTAA"
        )
        self.assertEqual(change_type, "indel")
        self.assertIsNone(groups[0])
        self.assertIsNone(groups[1])
        self.assertIsNone(groups[2])
