import unittest
from mutation2seq import genefeatures as gf

class TestGtfGff(unittest.TestCase):

    def test_init(self):
        gtf = gf.GtfGff()
        self.assertEqual(len(gtf), 0)

if __name__ == '__main__':
    unittest.main()
