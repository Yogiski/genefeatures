import unittest
from intervaltree import Interval, IntervalTree
from genefeatures import gtf_tools as gt
from genefeatures import gene_feature as gfeat

class TestGeneFeature(unittest.TestCase):

    def setUp(self):
        gtf = gt.parse_gtf("tests/data/test_hs_grch38.gtf")
        self.gtf = gtf.query({"attributes": {"gene_name": "PRDM16"}})
        self.gf = gfeat.GeneFeature(self.gtf)
    
    def test_add_record(self):
        self.assertIsInstance(self.gf, gfeat.GeneFeature)
        empty = gfeat.GeneFeature()
        recs = self.gtf.export_records()

        gf = empty.add_records(recs)
        gf = gfeat.GeneFeature(recs[0])
    
    def test_get_interval_attr(self):
        inter = self.gf.locations.all_intervals.pop()
        feat = self.gf._get_interval_attr(inter, "feature")
        self.assertIsInstance(feat, str)
        self.assertIn(feat, ["exon", "cds", "start_codon", "stop_codon", "gene", "transcript"])
    
    def test_partition_transcripts(self):
        gf = self.gf
        gf.partition_transcripts()
        ids = gf.transcript_ids
        self.assertIsInstance(gf.transcripts[ids[0]], IntervalTree)
        self.assertNotEqual(gf.transcripts[ids[0]], gf.transcripts[ids[1]])
        self.assertTrue(False)

if __name__ == '__main__':
    unittest.main()