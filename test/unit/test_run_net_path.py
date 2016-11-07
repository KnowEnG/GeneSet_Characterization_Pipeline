import yaml
import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import pandas as pd
import numpy as np


class TestRun_net_path(TestCase):
    def setUp(self):
        self.run_parameters = {"gg_network_name_full_path":"../../data/networks/TEST_1_gene_gene.edge",
                            "pg_network_name_full_path":"../../data/networks/TEST_1_property_gene.edge",
                            "spreadsheet_name_full_path":"../../data/spreadsheets/TEST_1_spreadsheet.tsv",
                            "gene_names_map": "../../data/spreadsheets/TEST_1_spreadsheet_MAP.tsv",
                            "results_directory":"./tmp",
                            "rwr_max_iterations":500,
                            "rwr_convergence_tolerence":1.0e-4,
                            "rwr_restart_probability":0.5,        
                            "k_space":2}
        # with open("unit_test_run_dir/net_path.yml", 'r') as file_handle:
        #     self.run_parameters = yaml.load(file_handle)

    def tearDown(self):
        del self.run_parameters

    def test_net_path(self):
    
        ret = tl.run_net_path(self.run_parameters)
        ret.index = np.arange(ret.shape[0])
        res = pd.DataFrame([['P7', 'P6'],['P6', 'P7']], columns=['user_gene_set1', 'user_gene_set2'])
        comp = ret.equals(res)
        self.assertEqual(True, comp)

if __name__ == '__main__':
    unittest.main()

