import os
import unittest
from unittest import TestCase
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import geneset_characterization_toolbox as tl
import numpy as np

class TestPerform_net_path(TestCase):
    def setUp(self):
        self.network_sparse = np.array([[0,0.5,1,0,1],[0.5,0,1,1,0],[1,1,0,1,1],[0,1,1,0,0],[1,0,1,0,0]])
        self.spreadsheet_df = pd.DataFrame([[0,1],[1,0],[1,1]], columns=['user1', 'user2'], index=['G1','G2','G3'])
        self.unique_gene_names = ['G1','G2','G3']
        self.pg_network_n1_names = ['P6', 'P7']
        self.run_parameters = {'rwr_max_iterations': 500,
                               'rwr_convergence_tolerence': 0.000001, 'rwr_restart_probability': 0.5, 'k_space': 2,
                               'results_directory': "unit_test_run_dir/results"}

    def tearDown(self):
        del self.network_sparse
        del self.spreadsheet_df
        del self.unique_gene_names
        del self.pg_network_n1_names
        del self.run_parameters

    def test_perform_net_path(self):
        ret = tl.get_net_path_results(self.spreadsheet_df, self.network_sparse, self.unique_gene_names,
                               self.pg_network_n1_names, self.run_parameters)
        ret.index = np.arange(ret.shape[0])
        res = pd.DataFrame({'user1': ['P6', 'P7'], 'user2': ['P7', 'P6']})
        comp = ret.equals(res)
        self.assertEqual(True, comp)


if __name__ == '__main__':
    unittest.main()
