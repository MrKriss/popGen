'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import glob

import socket

from preprocess import  Workflow, ConfigClass

      
#==============================================================================
''' RUN SCRIPT FOR ALLL READS IN stickleback RAD data '''
#===============================================================================


#===============================================================================
# Setup Configuration
#===============================================================================
starting_dir = os.getcwd()

c = ConfigClass()

# Work out where data is stored
if socket.gethostname() == 'yildun':
    prefix = '/space/musselle/datasets'
elif socket.gethostname() == 'luca':
    prefix = '/home/musselle/san/data'

# Set paths 
c.data_inpath =  os.path.join(prefix,'sticklebacks') 
c.barcode_inpath = os.path.join(prefix,'sticklebacks/barcodes')
c.filtered_outpath = os.path.join(prefix,'sticklebacks/filtered_data')
c.processed_outpath = os.path.join(prefix,'sticklebacks/filtered_data')
c.clusters_outpath = os.path.join(prefix,'sticklebacks/clusters')

# Setup input files and barcodes
os.chdir(c.data_inpath)
#raw_files = glob.glob('*0.fastq.bgzf')
raw_files = glob.glob('sb_testdata.bgzf')
raw_files.sort()
c.raw_input_files = raw_files 

os.chdir(c.barcode_inpath)
barcodes = glob.glob('sb_testdata.txt')
barcodes.sort()
c.barcode_files = barcodes
os.chdir(starting_dir)

# Set barcode file mode  
c.barcode_files_setup = 'individual' # Each reads file has an associated barcode file 

# MIDtags
c.cutsite = 'TGCAGG'
c.max_edit_dist = 2
        
# FILTERING
# Whether to log reads that fail the filtering         
c.log_fails = True
       
# Define Class
Experiment = Workflow(c) 

#===============================================================================
# Setup and run filter
#===============================================================================
Experiment.filter_functions = [Experiment.make_propN_filter(0.1),
                               Experiment.make_phred_filter(25),
                               Experiment.make_cutsite_filter(max_edit_dist=2),
                               Experiment.make_overhang_filter('TCGAGG', 'GG', max_edit_dist=0)]

Experiment.filter_reads_pipeline()

#===============================================================================
# Process and Correct MID tag 
#===============================================================================
Experiment.process_MIDtag(max_edit_dist = 1, outfile_postfix='-clean')

#===============================================================================
# Cluster Data 
#===============================================================================
allreads_file = 'sb_' + 'allreads_preprocessed.fasta'
Experiment.trim_reads(outfile=allreads_file, n = 1)

# default Vars for clustering 
default_vars = { 'c_thresh' : 0.90,
                 'n_filter' : 8,
                 'threads' : 1,
                 'mem' : 0,
                 'mask' : False, 
                 'log_filename' : 'report.log'}
                
experiment_name = 'sb_clustered_reads'

# Variations to run
clustering_runs = [ { 'c_thresh' : 0.95},
                    { 'c_thresh' : 0.95, 'mask' : True},
                    { 'c_thresh' : 0.90},
                    { 'c_thresh' : 0.90, 'mask' : True},
                    { 'c_thresh' : 0.85},
                    { 'c_thresh' : 0.85, 'mask' : True},
                   ]
                   
for d in clustering_runs:
    
    inputs_dict = {}
    inputs_dict.update(default_vars)
    inputs_dict.update(d)
    
    dirname = experiment_name
    outfile = experiment_name
    if 'c_thresh' in d:
        dirname = dirname + '-c{:d}'.format(d['c_thresh']*100)
        outfile = outfile + '-c{:d}'.format(d['c_thresh']*100)
    if 'mask' in d:
        dirname = dirname + '-maskN'
        outfile = outfile + '-maskN'
    
    path = os.path.join(c.clusters_outpath, dirname)        
    if not os.path.exists(path):
        os.mkdir(path)
        
    path2outfile  = os.path.join(path, outfile)

    Experiment.run_cdhit_clustering(infile=allreads_file, outfile=path2outfile,
              **inputs_dict)


## Display Summary
#summary(clustered_file)
