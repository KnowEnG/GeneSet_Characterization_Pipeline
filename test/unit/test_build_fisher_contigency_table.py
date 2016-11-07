import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import numpy as np

class TestBuild_fisher_contingency_table(TestCase):
    def setUp(self):
        self.overlap_count = 1
        self.user_count = 2
        self.gene_count = 3
        self.count = 4

    def tearDown(self):
        del self.overlap_count
        del self.user_count
        del self.gene_count
        del self.count

    def test_build_fisher_contigency_table(self):
        ret = tl.build_fisher_contingency_table(self.overlap_count, self.user_count,
                                               self.gene_count, self.count)
        self.assertEqual(np.array_equal(ret, np.matrix([[1, 1], [2, 0]])), True)

if __name__ == '__main__':
    unittest.main()
