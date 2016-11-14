import yaml
import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import pandas as pd
import numpy as np

class TestRun_fisher(TestCase):
    def setUp(self):
        self.run_parameters = {"gg_network_name_full_path":"../../data/networks/TEST_1_gene_gene.edge",
                               "pg_network_name_full_path":"../../data/networks/TEST_1_property_gene.edge",
                               "spreadsheet_name_full_path":"../../data/spreadsheets/TEST_1_spreadsheet.tsv",
                               "gene_names_map": "../../data/spreadsheets/TEST_1_spreadsheet_MAP.tsv",
                               "results_directory":"./tmp"}

    def tearDown(self):
        del self.run_parameters

    def test_run_fisher(self):
        ret = tl.run_fisher(self.run_parameters)
        ret['pval'] = ret['pval'].round(4)
        ret.index = np.arange(ret.shape[0])
        data1 = [['user_gene_set1', 'P7', 0.5229, 5, 2, 3, 2], ['user_gene_set2', 'P6', 0.2218, 5, 1, 3, 1], ['user_gene_set2', 'P7', 0.0000, 5, 1, 3, 0],['user_gene_set1', 'P6', 0.0000, 5, 2, 3, 0]]
        res1 = pd.DataFrame(data1)
        res1.columns = ["user_gene_set", "property_gene_set", "pval", "universe_count",
                       "user_count", "property_count", "overlap_count"]

        data2 = [['user_gene_set1', 'P7', 0.5229, 5, 2, 3, 2], ['user_gene_set2', 'P6', 0.2218, 5, 1, 3, 1], ['user_gene_set1', 'P6', 0.0000, 5, 2, 3, 0], ['user_gene_set2', 'P7', 0.0000, 5, 1, 3, 0]]
        res2 = pd.DataFrame(data2)
        res2.columns = ["user_gene_set", "property_gene_set", "pval", "universe_count",
                       "user_count", "property_count", "overlap_count"]
        
        comp = ret.equals(pd.DataFrame(res1)) or ret.equals(pd.DataFrame(res2))
        self.assertEqual(True, comp)

if __name__ == '__main__':
    unittest.main()