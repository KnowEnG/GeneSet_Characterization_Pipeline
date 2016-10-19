import yaml
import unittest
from unittest import TestCase
import geneset_characterization_toolbox as tl
import pandas as pd
import numpy as np


class TestRun_DRaWR(TestCase):
    def setUp(self):
        with open("unit_test_run_dir/DRaWR_run_file.yml", 'r') as file_handle:
            self.run_parameters = yaml.load(file_handle)

    def tearDown(self):
        del self.run_parameters

    def test_run_DRaWR(self):
        ret = tl.run_DRaWR(self.run_parameters)
        ret.index = np.arange(ret.shape[0])
        res = pd.DataFrame(['P7', 'P6'], columns=['drug1'])
        comp = ret.equals(res)
        self.assertEqual(True, comp)

if __name__ == '__main__':
    unittest.main()