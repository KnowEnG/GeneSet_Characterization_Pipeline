import yaml
import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import pandas as pd
import numpy as np

class TestRun_fisher(TestCase):
    def setUp(self):
        self.run_parameters = {"gg_network_name_full_path":"../../data/networks/test1_gene_gene.edge.txt",
                               "pg_network_name_full_path":"../../data/networks/test1_property_gene.edge.txt",
                               "spreadsheet_name_full_path":"../../data/spreadsheets/test1_spreadsheet.tsv",
                               "results_directory":"../../data/verification"}
        # with open("unit_test_run_dir/fisher.yml", 'r') as file_handle:
        #     self.run_parameters = yaml.load(file_handle)

    def tearDown(self):
        del self.run_parameters

    def test_run_fisher(self):
        ret = tl.run_fisher(self.run_parameters)
        ret['pval'] = ret['pval'].round(4)
        ret.index = np.arange(ret.shape[0])
        data = [['user1', 'P7', 1.204, 5, 2, 3, 2], ['user1', 'P6', 0.0000, 5, 2, 3, 0]]
        res = pd.DataFrame(data)
        res.columns = ["user_gene_set", "property_gene_set", "pval", "universe_count",
                       "user_count", "property_count", "overlap_count"]
        self.assertEqual(True, ret.equals(pd.DataFrame(res)))

if __name__ == '__main__':
    unittest.main()