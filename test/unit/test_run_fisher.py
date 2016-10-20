import yaml
import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import pandas as pd
import numpy as np

class TestRun_fisher(TestCase):
    def setUp(self):
        with open("unit_test_run_dir/fisher.yml", 'r') as file_handle:
            self.run_parameters = yaml.load(file_handle)

    def tearDown(self):
        del self.run_parameters

    def test_run_fisher(self):
        ret = tl.run_fisher(self.run_parameters)
        ret['pval'] = ret['pval'].round(4)
        ret.index = np.arange(ret.shape[0])
        data = [['drug1', 'P1', 3, 2, 1, 1, 0.6667], ['drug1', 'P2', 3, 2, 3, 2, 1.0000]]
        res = pd.DataFrame(data)
        res.columns = ['user gene', 'property', 'count', 'user count',
                       'gene count', 'overlap', 'pval']
        self.assertEqual(True, ret.equals(pd.DataFrame(res)))

if __name__ == '__main__':
    unittest.main()