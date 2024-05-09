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
        self.assertEqual(len(set(seqname)), 1)

        # take int as arg
        int_seqname = self.gtf.get_records_by_seqname(1)
        self.assertGreater(len(int_seqname), 0)

    def test_get_records_by_attribute(self):

        records = self.gtf.get_records_by_attribute({"gene_name": "PRDM16"})
        self.assertGreater(len(records), 0)
        feature = [d["feature"] for d in records]
        self.assertGreater(len(set(feature)), 0)
    
        records = self.gtf.get_records_by_attribute({"transcript_id": "ENST00000511072"})
        ids = [r["attributes"]["transcript_id"] for r in records]
        self.assertEqual(len(set(ids)), 1)

    
    def test_gtf_gff_from_records(self):
        new_gtf = self.gtf.gtf_gff_from_records(self.gtf[0])
        self.assertEqual(len(new_gtf), 1)
        self.assertTrue(isinstance(new_gtf, gf.GtfGff))
        new_gtf = self.gtf.gtf_gff_from_records(self.gtf[0:10])
        self.assertEqual(len(new_gtf), 10)
    
    def test_process_query_single_operator(self):

        gtf = self.gtf
        condition = {"AND":
            {"feature": "start_codon", "attributes": {"transcript_id": "ENST00000511072"}}
        }
        hashes = gtf._process_query(condition)
        records = gtf._get_records(hashes)
        feat = [r["feature"] for r in records]
        self.assertEqual(len(set(feat)), 1)
        attrs = [r["attributes"]["transcript_id"] for r in records]
        self.assertEqual(len(set(attrs)), 1)

        condition = {"AND":
            {"feature": "start_codon", "attributes": {"transcript_id": ["ENST00000511072", "ENST00000378391"]}}
        }
        hashes = gtf._process_query(condition)
        records = gtf._get_records(hashes)
        feat = [r["feature"] for r in records]
        self.assertEqual(len(set(feat)), 1)
        attrs = [r["attributes"]["transcript_id"] for r in records]
        self.assertEqual(len(set(attrs)), 2)
    
    def test_process_query_or_logic(self):

        gtf = self.gtf
        condition = {
            "OR":[
                {"AND": {"feature": "start_codon", "attributes": {"transcript_id": "ENST00000511072"}}},
                {"AND": {"feature": "start_codon", "attributes": {"transcript_id": "ENST00000378391"}}}
            ]
        }
        hashes = gtf._process_query(condition)
        records = gtf._get_records(hashes)
        feat = [r["feature"] for r in records]
        self.assertEqual(len(set(feat)), 1)
        attrs = [r["attributes"]["transcript_id"] for r in records]
        self.assertEqual(len(set(attrs)), 2)

        condition = {
            "AND": {
                "feature": "start_codon",
                "OR": [
                    {"attributes": {"transcript_id": "ENST00000511072"}},
                    {"attributes": {"transcript_id": "ENST00000378391"}}
                ]
            }
        }
        hashes = gtf._process_query(condition)
        records = gtf._get_records(hashes)
        feat = [r["feature"] for r in records]
        attrs = [r["attributes"]["transcript_id"] for r in records]
        self.assertEqual(len(set(feat)), 1)
        self.assertEqual(len(set(attrs)), 2)

    def test_process_query_not_logic(self):

        gtf = self.gtf
        condition = {
            "feature": "start_codon",
            "NOT": {"attributes": {"transcript_id": "ENST00000511072"}}
        }
        hashes = gtf._process_query(condition)
        records = gtf._get_records(hashes)
        feat = [r["feature"] for r in records]
        self.assertNotIn("ENST00000511072", feat)
        attrs = [r["attributes"]["transcript_id"] for r in records]
        self.assertEqual(len(set(feat)), 1)
        self.assertEqual(len(set(attrs)), 2)


if __name__ == '__main__':
    unittest.main()
