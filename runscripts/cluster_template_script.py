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
''' EXPERIMENT NAME'''
#===============================================================================

#===============================================================================
# Setup Configuration
#===============================================================================
starting_dir = os.getcwd()

c = ConfigClass()

c.experiment_name = 'my_experiment'

# Work out where data is stored
if socket.gethostname() == 'yildun':
    prefix = '/space/musselle/datasets'
elif socket.gethostname() == 'luca':
    prefix = '/home/musselle/san/data'
elif socket.gethostname() == 'gg-pc6':
    prefix = '/home/musselle/data'

# Set paths for IO
c.data_inpath =  os.path.join(prefix,'sticklebacks') 
c.barcode_inpath = os.path.join(prefix,'sticklebacks/barcodes')
c.filtered_outpath = os.path.join(prefix,'sticklebacks/filtered_data')
c.processed_outpath = os.path.join(prefix,'sticklebacks/filtered_data')
c.clusters_outpath = os.path.join(prefix,'sticklebacks/clusters')

# Setup input files glob
#os.chdir(c.data_inpath)
#raw_files = glob.glob('*[0-9].fastq.bgzf')
#raw_files.sort()
#c.raw_input_files = raw_files 

# Setup barcode files glob
#os.chdir(c.barcode_inpath)
#barcodes = glob.glob('*[0-9].txt')
#barcodes.sort()
#c.barcode_files = barcodes
#os.chdir(starting_dir)

# Set barcode file mode  
#c.barcode_files_setup = 'individual' # Each reads file has an associated barcode file 
#
# MIDtags
#c.cutsite = 'TGCAGG'
#c.max_edit_dist = 2
        
# FILTERING
# Whether to log reads that fail the filtering         
#c.log_fails = True
       
# Define Classes
#Preprocess = Preprocessor(c) 

try:
    cluster_file_path
except NameError:
    cluster_file_path = os.path.join(c.processed_outpath, c.experiment_name + '_all_preprocessed.fasta')

#===============================================================================
# Cluster Data 
#===============================================================================

# default Vars for clustering 
#default_vars = { 'c_thresh' : 0.90,
#                 'n_filter' : 8,
#                 'threads' : 1,
#                 'mem' : 0,
#                 'maskN' : False}

# Variations to run
batch_parameters = [ { 'c_thresh' : 0.95},
                    { 'c_thresh' : 0.95, 'maskN' : True},
                    { 'c_thresh' : 0.90},
                    { 'c_thresh' : 0.90, 'maskN' : True},
                    { 'c_thresh' : 0.85},
                    { 'c_thresh' : 0.85, 'maskN' : True},
                   ]
                   
Clusterer = Clustering(c, cluster_file_path) 

Clusterer.run_batch_cdhit_clustering(batch_parameters, threads=10)

## Display Summary
#summary(clustered_file)
