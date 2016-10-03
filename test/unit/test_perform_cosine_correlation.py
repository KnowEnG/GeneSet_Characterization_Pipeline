import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import numpy as np
import pandas as pd

class TestPerform_cosine_correlation(TestCase):
    def setUp(self):
        self.g_newspace_matrix = np.array([[6, 3], [-12, 5]])
        self.p_newspace_matrix = np.array([[-3, -2]])
        self.gene_names = ['g1', 'g2']
        self.property_name = ['p1']

    def tearDown(self):
        del self.g_newspace_matrix
        del self.p_newspace_matrix
        del self.gene_names
        del self.property_name

    def test_perform_cosine_correlation(self):
        cosine_matrix_df = tl.perform_cosine_correlation(
            self.g_newspace_matrix, self.p_newspace_matrix, self.gene_names, self.property_name)
        cosine_matrix_df = cosine_matrix_df.round(6)
        ret = pd.DataFrame([-0.992278, 0.554700], index=['g1', 'g2'], columns=['p1'])
        self.assertEqual(ret.equals(cosine_matrix_df), True)


if __name__ == '__main__':
    unittest.main()

