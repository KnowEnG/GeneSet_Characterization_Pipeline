import os
import unittest
from unittest import TestCase
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import geneset_characterization_toolbox as tl

class TestPerform_fisher_exact_test(TestCase):
    def setUp(self):
        self.prop_gene_network_sparse = csr_matrix([[1, 1], [0, 1], [0, 1]])
        self.sparse_dict = {0: 'p1', 1: 'p2'}
        self.spreadsheet_df = pd.DataFrame({'GS1': [1, 1, 0]},
                                           index=['G1', 'G2', 'G3'])
        self.results_dir = "unit_test_run_dir/results"

    def tearDown(self):
        del self.prop_gene_network_sparse
        del self.sparse_dict
        del self.spreadsheet_df
        del self.results_dir

    def test_perform_fisher_exact_test(self):
        ret = tl.perform_fisher_exact_test(self.prop_gene_network_sparse,
                                           self.sparse_dict, self.spreadsheet_df,
                                           self.results_dir)
        ret.columns = np.arange(ret.shape[1])
        ret.index = np.arange(ret.shape[0])
        ret[6] = ret[6].round(4)
        data = [['GS1', 'p1', 3, 2, 1, 1, 0.6667], ['GS1', 'p2', 3, 2, 3, 2, 1.0000]]
        res = pd.DataFrame(data)
        comp = ret.equals(pd.DataFrame(res))
        self.assertEqual(True, comp)

if __name__ == '__main__':
    unittest.main()