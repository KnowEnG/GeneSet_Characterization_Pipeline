import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import numpy as np
import pandas as pd


class TestRank_property(TestCase):
    def setUp(self):
        self.final_spreadsheet_df = pd.DataFrame([[1, 1.2, 0.5], [1, 0, 0.6], [0, 1, 0.5]],
         index=['G1', 'p1', 'p2'], columns=['gs1', 'gs2', 'base'])
        self.start_index = 1
        self.end_index = 3
        self.run_parameters = {"gene_names_map": "../../data/spreadsheets/TEST_spreadsheet_MAP.tsv"}
    def tearDown(self):
        del self.final_spreadsheet_df
        del self.start_index
        del self.end_index

    def test_rank_property(self):
        ret1 = tl.construct_drawr_result_df(self.final_spreadsheet_df, self.start_index, self.end_index, False, self.run_parameters)
        res1 = pd.DataFrame([['gs2', 'p2',0.5,1,0.5], ['gs1','p1',0.4,1,0.6],['gs1','p2',-0.5,0,0.5],['gs2','p1',-0.6,0,0.6]])
        res1.columns = ['user_gene_set', 'property_gene_set', 'difference_score', 'query_score', 'baseline_score']
        self.assertEqual(True, np.array_equal(ret1.values, res1.values))


        ret2 = tl.construct_drawr_result_df(self.final_spreadsheet_df, 0, 1, True, self.run_parameters)
        res2 = pd.DataFrame([['gs2', 'orig_G1',0.7,1.2,0.5], ['gs1','orig_G1',0.5,1,0.5]])
        res2.columns = ['user_gene_set', 'property_gene_set', 'difference_score', 'query_score', 'baseline_score']
        self.assertEqual(True, np.array_equal(ret2.values, res2.values))



if __name__ == '__main__':
    unittest.main()

