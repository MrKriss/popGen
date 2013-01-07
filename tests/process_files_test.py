'''
Created on 21 Nov 2012

@author: musselle
'''
import os
import sys

import unittest

from utils import Cycler

class Test(unittest.TestCase):


    def setUp(self):
        os.chdir('/space/musselle/datasets/gazellesAndZebras')
        
        
        RecCycler = Cycler()


    def tearDown(self):
        pass


    def testOpenCloseFile(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()