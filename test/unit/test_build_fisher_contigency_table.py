import unittest
from unittest import TestCase
import knpackage.toolbox as kn
import numpy as np

class TestBuild_fisher_contigency_table(TestCase):
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
        ret = kn.build_fisher_contigency_table(self.overlap_count, self.user_count,
                                               self.gene_count, self.count)
        self.assertEqual(ret.all(), np.matrix([[1, 1], [2, 0]]).all())

if __name__ == '__main__':
    unittest.main()