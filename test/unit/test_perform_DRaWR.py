import os
import unittest
from unittest import TestCase
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import geneset_characterization_toolbox as tl
import numpy as np


class TestPerform_DRaWR(TestCase):
    def setUp(self):
        self.network_sparse = csr_matrix(
            ([0.14999999999999999, 0.16250000000000001, 0.14999999999999999, 0.1875,
              0.16250000000000001, 0.1875, 0.083333333333333301, 0.083333333333333301,
              0.083333333333333301, 0.083333333333333301, 0.083333333333333301,
              0.083333333333333301, 0.083333333333333301, 0.083333333333333301,
              0.083333333333333301, 0.083333333333333301, 0.083333333333333301,
              0.083333333333333301],
             ([1, 2, 0, 2, 0, 1, 6, 5, 5, 6, 5, 6, 1, 2, 3, 0, 2, 4],
              [0, 0, 1, 1, 2, 2, 0, 1, 2, 2, 3, 4, 5, 5, 5, 6, 6, 6])), shape=(7, 7))
        self.spreadsheet_df = pd.DataFrame([[1, 1], [0, 1], [0, 1], [0, 1],
                                            [0, 1], [0, 0], [0, 0]], columns=['GS1', 'base'])

        self.spreadsheet_df.index = ['G1', 'G2', 'G3', 'G4', 'G5', 'P6', 'P7']
        self.len_gene_names = 5
        self.run_parameters = {'rwr_max_iterations': 500,
                          'rwr_convergence_tolerence': 0.0001, 'rwr_restart_probability': 0.5,
                          'results_directory': "unit_test_run_dir/results"}
    def tearDown(self):
        del self.network_sparse
        del self.spreadsheet_df
        del self.len_gene_names
        del self.run_parameters

    def test_perform_DRaWR(self):
        ret = tl.perform_DRaWR(self.network_sparse, self.spreadsheet_df,
                               self.len_gene_names, self.run_parameters)
        ret.index = np.arange(ret.shape[0])
        res = pd.DataFrame({'GS1': ['P7', 'P6']})
        comp = ret.equals(res)
        self.assertEqual(True, comp)


if __name__ == '__main__':
    unittest.main()