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
        self.assertEqual(
            str(context.exception),
            f"Expected types int, slice, or list; got '{idx}' of type: {type(idx)}"
        )
    def test_get_records(self):

        gtf = self.gtf

        single = gtf._get_records(gtf._record_hashes[0])
        self.assertTrue(isinstance(single, list))
        self.assertEqual(len(single), 1)

        multi = gtf._get_records(gtf._record_hashes[0:10])
        self.assertTrue(isinstance(multi, list))
        self.assertEqual(len(multi), 10)
    
    def test_lookup_hash(self):

        gtf = self.gtf
        # feature look up single
        single = gtf._lookup_hash(gtf.feature_index, "CDS")
        self.assertTrue(isinstance(single, list))
        self.assertGreater(len(single), 1)

        # feature look up multiple
        multiple = gtf._lookup_hash(gtf.feature_index, ["CDS", "exon"])
        self.assertTrue(isinstance(multiple, list))
        self.assertGreater(len(multiple), 2)

    def test_get_records_by_feature(self):

        # get one feature type
        records = self.gtf.get_records_by_feature("CDS")
        feature = [d["feature"] for d in records]
        self.assertEqual(len(set(feature)), 1)

        # get mulitple feature types
        records = self.gtf.get_records_by_feature(["CDS", "exon"])
        feature = [d["feature"] for d in records]
        self.assertEqual(len(set(feature)), 2)

    def test_get_records_by_seqname(self):

        records = self.gtf.get_records_by_seqname("1")
        self.assertGreater(len(records), 0)

        self.assertTrue(isinstance(records, list))
        seqname = [d["seqname"] for d in records]
        self.assertGreater(len(seqname), 0)

    def test_get_records_by_attribute(self):

        records = self.gtf.get_records_by_attribute({"exon_number": "1"})
        self.assertGreater(len(records), 0)
        feature = [d["feature"] for d in records]
        self.assertGreater(len(set(feature)), 1)
        records = self.gtf.get_records_by_attribute({"exon_number": ["1", "2"]})
        self.assertGreater(len(records), 0)

if __name__ == '__main__':
    unittest.main()
