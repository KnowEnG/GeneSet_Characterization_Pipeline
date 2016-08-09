from unittest import TestCase
import unittest
import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
from sklearn.preprocessing import normalize
from geneset_characterization_toolbox import build_fisher_contigency_table
from geneset_characterization_toolbox import perform_fisher_exact_test
from geneset_characterization_toolbox import perform_DRaWR
from geneset_characterization_toolbox import append_column_to_spreadsheet
from geneset_characterization_toolbox import smooth_matrix_with_rwr

class testBuild_geneSet_characterization_toolbox(TestCase):
    def test_build_fisher_contingency_table(self):
        ret = build_fisher_contigency_table(1,2,3,4)
        self.assertEqual(ret.all(),np.matrix([[1,1],[2,0]]).all())

    def test_perform_fisher_exact_test(self):
        param1 = csr_matrix([[1,1],[0,1],[0,1]])
        param2 = {0: 'p1',1: 'p2'}
        param3 = pd.DataFrame({'GS1':[1,1,0]},index=['G1','G2','G3'])
        param4 = '/'
        ret = perform_fisher_exact_test(param1,param2,param3,param4)
        ret.columns = np.arange(7)
        data = [['GS1','p1',3,2,1,1,0.667],['GS1','p2',3,2,3,2,1.0]]
        comp = ret.equals(pd.DataFrame(data))
        self.assertEqual(True,comp)

    def test_append_column_to_spreadsheet(self):
        spreadsheet_df = pd.DataFrame([1,0,0,0,0,0,0],columns=['GS1'])
        len_gene = 5
        ret = append_column_to_spreadsheet(spreadsheet_df,len_gene)
        comp = ret.equals(pd.DataFrame({'GS1': [1,0,0,0,0,0,0],'base': [1,1,1,1,1,0,0]}))
        self.assertEqual(True,comp)

    def test_smooth_matrix_with_rwr(self):
        restart = np.transpose(np.array([[1,0,0,0,0,0,0],[1,1,1,1,1,0,0]]))
        network_sparse = csr_matrix(
                        ([0.14999999999999999,0.16250000000000001,0.14999999999999999,0.1875,
                          0.16250000000000001,0.1875,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301],
                        ([1,2,0,2,0,1,6,5,5,6,5,6,1,2,3,0,2,4],
                         [0,0,1,1,2,2,0,1,2,2,3,4,5,5,5,6,6,6])),shape=(7,7))
        run_parameters = {'number_of_iteriations_in_rwr': 500,'it_max': 10000,
                          'restart_tolerance': 0.0001,'restart_probability': 0.5}
        ret = smooth_matrix_with_rwr(normalize(restart,norm='l1',axis=0),
              normalize(network_sparse,norm='l1',axis=0),run_parameters)[0]
        ret = np.around(ret, decimals=3)
        res = np.transpose(np.array([[0.5646,0.142,0.1658,0.005,0.0132,0.0299,0.0794],
              [0.1824,0.1883,0.2106,0.1156,0.1157,0.0934,0.094]]))
        comp = (ret==np.round(res, decimals=3)).all()
        self.assertEqual(True,comp)

    def test_perform_DRaWR(self):
        network_sparse = csr_matrix(
                        ([0.14999999999999999,0.16250000000000001,0.14999999999999999,0.1875,
                          0.16250000000000001,0.1875,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301,0.083333333333333301,0.083333333333333301,
                          0.083333333333333301],
                        ([1,2,0,2,0,1,6,5,5,6,5,6,1,2,3,0,2,4],
                         [0,0,1,1,2,2,0,1,2,2,3,4,5,5,5,6,6,6])),shape=(7,7))
        spreadsheet_df = pd.DataFrame({'GS1':[1,0,0,0,0,0,0]})
        len_gene_names = 5
        run_parameters = {'number_of_iteriations_in_rwr': 500,'it_max': 10000,
                          'restart_tolerance': 0.0001,'restart_probability': 0.5}
        ret = perform_DRaWR(network_sparse, spreadsheet_df, len_gene_names, run_parameters)


        
if __name__ == '__main__':
    unittest.main()