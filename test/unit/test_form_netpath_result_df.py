import os
import unittest
from unittest import TestCase
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import geneset_characterization_toolbox as tl
import numpy as np

class TestForm_netpath_result_df(TestCase):
    def setUp(self):
        self.cosine_matrix_df = pd.DataFrame([[1,2,3],[3,4,5],[7,8,9]], index=['g1','g2','g3'],columns=['p1','p2','p3'])
        self.spreadsheet_df = pd.DataFrame([[1,0],[0,1],[0,1]], index=['g1','g2','g3'],columns=['s1','s2'])

    def tearDown(self):
        del self.cosine_matrix_df
        del self.spreadsheet_df

    def test_perform_net_path(self):
        ret = tl.form_netpath_result_df(self.spreadsheet_df, self.cosine_matrix_df)
        ret.index = np.arange(ret.shape[0])
        res = pd.DataFrame([['s2','p3',14],['s2','p2',12],['s2','p1',10],['s1','p3',3],['s1','p2',2],['s1','p1',1]])
        res.columns = ['user_gene_set', 'property_gene_set', 'cosine_sum']
        print(res)
        print(ret)
        print(res==ret)
        print(res.columns==ret.columns)
        print(res.index==ret.index)
        comp = ret.equals(res)
        self.assertEqual(True, comp)


if __name__ == '__main__':
    unittest.main()