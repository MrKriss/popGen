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

from preprocess import Workflow, ConfigClass

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
c.inpath =  os.path.join(prefix,'sticklebacks') 
c.barpath = os.path.join(prefix,'sticklebacks/barcodes')
c.filteredpath = os.path.join(prefix,'sticklebacks/filtered_data')
c.processedpath = os.path.join(prefix,'sticklebacks/filtered_data')

# Setup input files and barcodes
os.chdir(c.inpath)
raw_files = glob.glob('*[1-9].fastq.bgzf')
raw_files.sort()
c.raw_input_files = raw_files 

os.chdir(c.barpath)
barcodes = glob.glob('*[1-9].fastq.bgzf')
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

## Variables 
#c_thresh = 0.9
#n_filter = 8
#
#clustered_file = 'lane' + LANE + 'clustered_reads'
#cluster_cdhit(infile=allreads_file, outfile=clustered_file,
#              c_thresh=c_thresh, n_filter=n_filter)
#
## Display Summary
#summary(clustered_file)
