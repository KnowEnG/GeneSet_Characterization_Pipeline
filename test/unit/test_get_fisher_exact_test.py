import os
import unittest
from unittest import TestCase
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import geneset_characterization_toolbox as tl

class TestGet_fisher_exact_test(TestCase):
    def setUp(self):
        self.prop_gene_network_sparse = csr_matrix([[1, 1], [0, 1], [0, 1]])
        self.sparse_dict = {0: 'p1', 1: 'p2'}
        self.spreadsheet_df = pd.DataFrame({'GS1': [1, 1, 0]},
                                           index=['G1', 'G2', 'G3'])

    def tearDown(self):
        del self.prop_gene_network_sparse
        del self.sparse_dict
        del self.spreadsheet_df


    def test_get_fisher_exact_test(self):
        ret = tl.get_fisher_exact_test(self.prop_gene_network_sparse,
                                           self.sparse_dict, self.spreadsheet_df)
        for i in range(2):
            ret[i][2] = format(ret[i][2], '.4f')
        data = [['GS1', 'p1', '0.1761', 3, 2, 1, 1], ['GS1', 'p2', '-0.0000', 3, 2, 3, 2]]
        self.assertEqual(ret, data)

if __name__ == '__main__':
    unittest.main()