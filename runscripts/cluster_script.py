'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import glob

import socket

from preprocess import  Preprocessor, ConfigClass
from cluster import Clustering

#==============================================================================
''' Clusrter gzL6'''
#===============================================================================

""" FOR LANE 6 """

LANE =  '6'

#===============================================================================
# Setup Configuration
#===============================================================================
starting_dir = os.getcwd()

c = ConfigClass()

c.experiment_name = 'gzL6'

# Work out where data is stored
if socket.gethostname() == 'yildun':
    prefix = '/space/musselle/datasets'
elif socket.gethostname() == 'luca':
    prefix = '/home/musselle/san/data'
elif socket.gethostname() == 'gg-pc6':
    prefix = '/home/musselle/data'

# Set paths for IO
# Set paths 
c.tag_processed_outpath = os.path.join(prefix,'gazelles-zebras/lane%s' % (LANE), 'filtered_data')
c.clusters_outpath = os.path.join(prefix,'gazelles-zebras/clusters')

cluster_file_path = os.path.join(c.tag_processed_outpath, c.experiment_name + '_all_preprocessed.fasta')

#===============================================================================
# Cluster Data 
#===============================================================================

# Variations to run
batch_parameters = [ { 'c_thresh' : 0.95, 'allvall' : True},
                    { 'c_thresh' : 0.90, 'allvall' : True},
                   ]
                   
Clusterer = Clustering(c, cluster_file_path) 

Clusterer.run_batch_cdhit_clustering(batch_parameters, threads=20, postfix='_allvall')

## Display Summary
#summary(clustered_file)
