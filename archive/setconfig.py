'''
Created on 3 Jul 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

from utils.general_utilities import get_path_prefix


import argparse
    
parser = argparse.ArgumentParser(description='Filter and clean up FastQ files.')
parser.add_argument('-r', '--rootdir', action='store', help='Root directory where all ')
parser.add_argument('-n', '--name', action='store', type=int, help='Name of the experiemnt, and subdirectory where data is stored.')
parser.add_argument('-i', '--input', action='store', type=int, help='Input file(s) to process. Will accept a glob')


parser.add_argument('minreads', action='store', type=int, help='Minimum value of threshold for total reads in clusters of any size')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
args = parser.parse_args()

class Config(object):
    pass
c = Config()
  
# Setup Configuration
prefix =  get_path_prefix()


# Default path locations



c.testing = testing
c.root_name = name
c.db_name = db_name
if testing:
    self.c.data_inpath =  os.path.join(prefix,name, 'testset')
else:
    self.c.data_inpath =  os.path.join(prefix, name, 'raw-data') 
c.barcode_inpath = os.path.join(prefix, name , 'barcodes')
c.filtered_outpath = os.path.join(prefix, name , 'processed-data')
c.tag_processed_outpath = os.path.join(prefix, name, 'processed-data')
c.tag_splitby_sample_outpath = os.path.join(prefix, name, 'processed-data', 'per-sample')
c.tag_splitby_subgroup_outpath = os.path.join(prefix, name, 'processed-data', 'per-subgroup')
c.clusters_outpath = os.path.join(prefix, name, 'clusters')
c.db_path = os.path.join(prefix,  name)
c.cdhit_path = os.path.expanduser("~/bin/cd-hit-v4.6.1")






















if __name__ == '__main__':
    pass