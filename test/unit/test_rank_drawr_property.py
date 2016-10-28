import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import numpy as np
import pandas as pd


class TestRank_property(TestCase):
    def setUp(self):
        self.final_spreadsheet_df = pd.DataFrame([[1, 1, 0.5], [1, 0, 0.5], [0, 1, 0.5]],
         index=['g1', 'p1', 'p2'], columns=['gs1', 'gs2', 'base'])
        self.pg_network_n1_names = ['p1', 'p2']
    def tearDown(self):
        del self.final_spreadsheet_df
        del self.pg_network_n1_names

    def test_rank_property(self):
        ret = tl.rank_drawr_property(self.final_spreadsheet_df, self.pg_network_n1_names)
        res = pd.DataFrame([[0.5, -0.5], [-0.5, 0.5]], index=['p1', 'p2'], columns=['gs1', 'gs2'])
if __name__ == '__main__':
    unittest.main()

