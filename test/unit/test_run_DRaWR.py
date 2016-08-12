import yaml
import unittest
from unittest import TestCase
import src.geneset_characterization_toolbox as tl
import pandas as pd


class TestRun_DRaWR(TestCase):
    def setUp(self):
        with open("tmp/DRaWR_run_file.yml", 'r') as file_handle:
            self.run_parameters = yaml.load(file_handle)

    def tearDown(self):
        del self.run_parameters

    def test_run_DRaWR(self):
        ret = tl.run_DRaWR(self.run_parameters)
        res = pd.DataFrame([['P7', 'P7'],['P6', 'P6']], columns=['drug1', 'base'])
        comp = ret.equals(res)
        self.assertEqual(True, comp)

if __name__ == '__main__':
    unittest.main()