'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import glob

import socket

from utils import get_data_prefix
from preprocess import  Preprocessor, ConfigClass
from cluster import Clustering

#==============================================================================
''' Filter/Clean/Cluster SCRIPT FOR ALLL READS IN Gazelles-Zebras RAD data Lane 6'''
#===============================================================================

""" FOR LANE 6 """

LANE =  '6'
experiment_name = 'gzL%s' % LANE

#===============================================================================
# Setup Configuration
#===============================================================================
starting_dir = os.getcwd()

c = ConfigClass()

c.experiment_name = experiment_name

# Work out where data is stored on this machine
prefix = get_data_prefix()

# Set paths 
c.data_inpath =  os.path.join(prefix,'gazelles-zebras/lane%s' % (LANE)) 
c.barcode_inpath = os.path.join(prefix,'gazelles-zebras/barcodes')
c.filtered_outpath = os.path.join(prefix,'gazelles-zebras/lane%s' % (LANE), 'filtered_data')
c.tag_processed_outpath = os.path.join(prefix,'gazelles-zebras/lane%s' % (LANE), 'filtered_data')
c.tag_split_outpath = os.path.join(prefix,'gazelles-zebras/lane%s' % (LANE), 'filtered_data', 'individuals')
c.clusters_outpath = os.path.join(prefix,'gazelles-zebras/clusters')

# Setup input files and barcodes
os.chdir(c.data_inpath)
raw_files = glob.glob('*[0-9].fastq.bgzf')
raw_files.sort()
c.raw_input_files = raw_files 

os.chdir(c.barcode_inpath)
barcodes = glob.glob('*%s.txt' % (LANE))
barcodes.sort()
c.barcode_files = barcodes
os.chdir(starting_dir)

# Set barcode file mode  
c.barcode_files_setup = 'global' # one global file(s) for all barcodes used 

# Set interim file suffixes
c.filtered_files_postfix = '-pass'
c.tag_processed_files_postfix = '-clean'

# MIDtags
c.cutsite = 'TGCAGG'
c.max_edit_dist = 2
        
# FILTERING
# Whether to log reads that fail the filtering         
c.log_fails = True
       
# Define Class
Preprocess = Preprocessor(c) 

#===============================================================================
# Setup and run filter
#===============================================================================
Preprocess.filter_functions = [Preprocess.make_propN_filter(0.1),
                               Preprocess.make_phred_filter(25),
                               Preprocess.make_cutsite_filter(max_edit_dist=2),
                               Preprocess.make_overhang_filter('TCGAGG', 'GG', max_edit_dist=0)]

Preprocess.filter_reads_pipeline()

#===============================================================================
# Process and Correct MID tag 
#===============================================================================

Preprocess.process_MIDtag(max_edit_dist = 1, outfile_postfix='-clean')
Preprocess.cleanup('filtered') # Remove filtered intermediate files 

out = Preprocess.trim_reads(n = 1)
Preprocess.cleanup('tag_processed') # Remove MID tag processed intermediate files 

#===============================================================================
# Split Data into individuals based on MID tag 
#===============================================================================








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
                   
Clusterer = Clustering(c, infiles=out[0], inpath=[1]) 

Clusterer.run_batch_cdhit_clustering(batch_parameters, threads=10)

## Display Summary
#summary(clustered_file)