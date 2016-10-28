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
        self.cosine_matrix_df = pd.DataFrame([[1,2,3],[3,4,5],[7,8,9]], index=['g1','g2','g3'],columns=['p1','p2','p3'])
        self.spreadsheet_df = pd.DataFrame([[1,0],[0,1],[0,1]], index=['g1','g2','g3'],columns=['s1','s2'])

    def tearDown(self):
        del self.cosine_matrix_df
        del self.spreadsheet_df

    def test_perform_net_path(self):
        ret = tl.rank_netpath_property(self.spreadsheet_df, self.cosine_matrix_df)
        ret.index = np.arange(ret.shape[0])
        res = pd.DataFrame([['p3','p3'],['p2','p2'],['p1','p1']],columns=['s1','s2'])
        comp = ret.equals(res)
        self.assertEqual(True, comp)


if __name__ == '__main__':
    unittest.main()