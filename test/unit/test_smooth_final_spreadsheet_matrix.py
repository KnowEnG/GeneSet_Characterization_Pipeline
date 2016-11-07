import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import pandas as pd
import numpy as np

class TestSmooth_final_spreadsheet_matrix(TestCase):
    def setUp(self):
        self.final_rwr_matrix = np.array([[2, 0], [4, 5]])
        self.ret = np.array([[1.60943791, 0], [2.19722458, 2.39789527]])


    def tearDown(self):
        del self.final_rwr_matrix
        del self.ret

    def test_smooth_final_spreadsheet_matrix(self):
        result = tl.smooth_final_spreadsheet_matrix(self.final_rwr_matrix, 2)
        self.assertEqual(np.array_equal(self.ret, np.round(result, 8)), True)


if __name__ == '__main__':
    unittest.main()
