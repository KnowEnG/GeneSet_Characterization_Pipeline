import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import numpy as np
import pandas as pd


class TestRank_property(TestCase):
    def setUp(self):
        self.spreadsheet_df = pd.DataFrame([[1, 1], [1, 0], [0, 1]], index=['g1', 'g2', 'g3'], columns=['gs1', 'gs2'])
        self.cosine_matrix_df = pd.DataFrame([[-0.6, 0.3], [0.7, 0], [0.8, -0.2]], index=['g1', 'g2', 'g3'], columns=['p1', 'p2'])
        self.result = pd.DataFrame([['p2', 'p1'], ['p1', 'p2']], index=[0, 1], columns=['gs1', 'gs2'])
    def tearDown(self):
        del self.spreadsheet_df
        del self.cosine_matrix_df

    def test_rank_property(self):
        ret = tl.rank_property(self.spreadsheet_df, self.cosine_matrix_df)
        self.assertEqual(ret.equals(self.result), True)

if __name__ == '__main__':
    unittest.main()

