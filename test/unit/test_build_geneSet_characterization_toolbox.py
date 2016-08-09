from unittest import TestCase
import unittest
import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
from geneset_characterization_toolbox import build_fisher_contigency_table, perform_fisher_exact_test

class testBuild_geneSet_characterization_toolbox(TestCase):
    def test_build_fisher_contingency_table(self):
        ret = build_fisher_contigency_table(1,2,3,4)
        self.assertEqual(ret.all(), np.matrix([[1, 1], [2, 0]]).all())

    def test_perform_fisher_exact_test(self):
        param1 = csr_matrix([[1,1], [0,1], [0,1]])
        param2 = {0: 'p1', 1: 'p2'}
        param3 = pd.DataFrame({'GS1':[1,1,0]}, index=['G1', 'G2', 'G3'])
        param4 = '/'
        ret = perform_fisher_exact_test(param1, param2, param3, param4)
        ret.columns = np.arange(7)
        data = [['GS1', 'p1', 3, 2, 1, 1, 0.667], ['GS1', 'p2', 3, 2, 3, 2, 1.0]]
        comp = ret.equals(pd.DataFrame(data))
        self.assertEqual(True, comp)

if __name__ == '__main__':
    unittest.main()