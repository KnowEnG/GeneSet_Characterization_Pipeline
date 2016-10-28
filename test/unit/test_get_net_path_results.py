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
        self.smooth_rwr_matrix = np.array([[3,2],[4,1],[5,3],[7,7],[8,2]])
        self.gene_length = 3
        self.run_parameters = {'k_space': 2}

    def tearDown(self):
        del self.smooth_rwr_matrix
        del self.gene_length
        del self.run_parameters


    def test_perform_net_path(self):
        ret = tl.get_net_path_results(self.gene_length, self.smooth_rwr_matrix, self.run_parameters)
        res = np.array([[0.9316, 0.7976],[0.5239, 1.0000], [0.8925, 0.8517]])
        self.assertEqual(np.array_equal(res, np.round(ret, 4)), True)

if __name__ == '__main__':
    unittest.main()

