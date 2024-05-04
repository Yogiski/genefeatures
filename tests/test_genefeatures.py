import unittest
from genefeatures import gtf_tools as gf

class TestGtfGff(unittest.TestCase):


    def setUp(self):

        self.empty = gf.GtfGff()
        self.gtf = gf.parse_gtf("tests/data/test_hs_grch38.gtf")


    def test_init(self):

        self.assertEqual(len(self.empty), 0)
        self.assertTrue(isinstance(self.empty, gf.GtfGff))
    

    def test_parse(self):

        self.assertEqual(type(self.gtf), gf.GtfGff)
        self.assertGreater(len(self.gtf), 0)
        full = gf.parse_gtf("tests/data/test_hs_grch38.gtf", gtf = self.empty)
        self.assertGreater(len(full), 0)


    def test_getitem(self):

        self.assertTrue(isinstance(self.gtf[0], dict))
        self.assertTrue(isinstance(self.gtf[[1, 3, 4]], list))
        self.assertTrue(isinstance(self.gtf[0:10], list))
        idx = "exon"
        with self.assertRaises(TypeError) as context:
            self.gtf[idx]
        print(context.exception)
        self.assertEqual(str(context.exception), f"Expected types int, slice, or list; got '{idx}' of type: {type(idx)}")


    def test_get_records_by_feature(self):

        records = self.gtf.get_records_by_feature("CDS")
        self.assertTrue(isinstance(records, list))
        feature = [d["feature"] for d in records]
        self.assertEqual(len(set(feature)), 1)


    def test_get_records_by_seqname(self):

        fail_case = self.gtf.get_records_by_seqname("does_not_exist")
        self.assertEqual(len(fail_case), 0)
        records = self.gtf.get_records_by_seqname("1")
        self.assertGreater(len(records), 0)
        self.assertTrue(isinstance(records, list))
        seqname = [d["seqname"] for d in records]
        self.assertEqual(len(set(seqname)), 1)
        int_seqname = self.gtf.get_records_by_seqname(1)
        self.assertGreater(len(int_seqname), 0)


    def test_get_records_by_attribute(self):

        records = self.gtf.get_records_by_attribute("gene_name", "PRDM16")
        self.assertGreater(len(records), 0)
        feature = [d["feature"] for d in records]
        self.assertGreater(len(set(feature)), 0)
    

if __name__ == '__main__':
    unittest.main()
