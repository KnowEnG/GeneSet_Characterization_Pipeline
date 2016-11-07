import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import numpy as np

class TestPerform_k_SVD(TestCase):
    def setUp(self):
        self.smooth_spreadsheet_matrix = np.array([[1, 5, 7], [2, 2, 4], [3, 4, 1]])
        self.k = 2

    def tearDown(self):
        del self.smooth_spreadsheet_matrix
        del self.k

    def test_perform_k_SVD(self):
        ret_U, ret_S = tl.calculate_k_SVD(self.smooth_spreadsheet_matrix, self.k)
        U = np.array([[-0.80896977, 0.39172363], [-0.44903089, 0.06944901], [-0.37939315, -0.91745814]])
        S = np.array([[3.24785452, 0], [0, 1.85452241]])
        self.assertEqual(np.array_equal(U, np.round(ret_U, 8)), True)
        self.assertEqual(np.array_equal(S, np.round(ret_S, 8)), True)

if __name__ == '__main__':
    unittest.main()
