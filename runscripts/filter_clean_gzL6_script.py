'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import glob

import socket

from preprocess import  Workflow, ConfigClass

      
#==============================================================================
''' RUN SCRIPT FOR ALLL READS IN stickleback RAD data '''
#===============================================================================

""" FOR LANE 6 """

LANE =  '6'

experiment_name = 'gzL6'

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
elif socket.gethostname() == 'gg-pc6':
    prefix = '/home/musselle/data'

# Set paths 
c.data_inpath =  os.path.join(prefix,'gazelles-zebras/lane%s' % (LANE)) 
c.barcode_inpath = os.path.join(prefix,'gazelles-zebras/barcodes')
c.filtered_outpath = os.path.join(prefix,'gazelles-zebras/filtered_data')
c.processed_outpath = os.path.join(prefix,'gazelles-zebras/filtered_data')
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
allreads_file = experiment_name + '_allreads_preprocessed.fasta'
Experiment.trim_reads(out_filename=allreads_file, n = 1)

# default Vars for clustering 
default_vars = { 'c_thresh' : 0.90,
                 'n_filter' : 8,
                 'threads' : 10,
                 'mem' : 0,
                 'maskN' : False}

# Variations to run
clustering_runs = [ { 'c_thresh' : 0.95},
                    { 'c_thresh' : 0.95, 'maskN' : True},
                    { 'c_thresh' : 0.90},
                    { 'c_thresh' : 0.90, 'maskN' : True},
                    { 'c_thresh' : 0.85},
                    { 'c_thresh' : 0.85, 'maskN' : True},
                   ]
                   
for d in clustering_runs:
    
    inputs_dict = {}
    inputs_dict.update(default_vars)
    inputs_dict.update(d)
    
    dirname = experiment_name + '_clustered_reads'
    outfile = experiment_name + '_clustered_reads'
    if 'c_thresh' in d:
        dirname = dirname + '-c{}'.format(int(d['c_thresh']*100))
        outfile = outfile + '-c{}'.format(int(d['c_thresh']*100))
    if 'maskN' in d:
        dirname = dirname + '-maskN'
        outfile = outfile + '-maskN'
    
    path = os.path.join(c.clusters_outpath, dirname)        
    if not os.path.exists(path):
        os.makedirs(path)
        
    path2outfile  = os.path.join(path, outfile)
    inputs_dict['log_filename'] = os.path.join(path, 'report.log')

    Experiment.run_cdhit_clustering(infile=allreads_file, outfile=path2outfile,
              **inputs_dict)


## Display Summary
#summary(clustered_file)
