import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import numpy as np


class TestProject_matrix_to_new_space_and_split(TestCase):
    def setUp(self):
        self.U = np.array([[2, 3], [-4, 5], [-1, -2]])
        self.S_full_squared_matrix = np.array([[3, 0], [0, 1]])
        self.unique_gene_length = 2

    def tearDown(self):
        del self.unique_gene_length
        del self.U
        del self.S_full_squared_matrix

    def test_project_matrix_to_new_space_and_split(self):
        ret_g, ret_p = tl.project_matrix_to_new_space_and_split(self.U, self.S_full_squared_matrix, self.unique_gene_length)

        self.assertEqual(ret_g.all(), np.array([[6, 3], [-12, 5]]).all())
        self.assertEqual(ret_p.all(), np.array([[-3, -2]]).all())


if __name__ == '__main__':
    unittest.main()
