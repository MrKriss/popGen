'''
Created on 5 Jul 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

class Config(object):
    pass
c = Config()


c.testing = False
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