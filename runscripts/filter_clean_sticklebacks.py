'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import glob

from preprocess import filter_reads, process_MIDtag, trim_reads, setup_cutsite_filter, \
setup_overhang_filter, setup_phred_filter, setup_propN_filter, setup_illumina_filter, \
filter_reads_pipeline, Workflow, ConfigClass

from cluster import cluster_cdhit, summary

       
#==============================================================================
''' RUN SCRIPT FOR ALLL READS IN stickleback RAD data '''
#===============================================================================


#===============================================================================
# Setup Configuration
#===============================================================================
starting_dir = os.getcwd()

c = ConfigClass()

# Set paths 
c.inpath = '/home/musselle/san/data/stickleback', 
c.barpath = '/home/musselle/san/data/stickleback/barcodes'
c.filteredpath = '/home/musselle/san/data/stickleback/filtered_data/'
c.processedpath = '/home/musselle/san/data/stickleback/filtered_data/processed'

# Set input files and barcodes
os.chdir(c.paths.inpath)
raw_files = glob.glob('*[0-9].fastq.bgzf')
raw_files.sort()
c.raw_input_files = raw_files 

os.chdir(c.paths.barpath)
barcodes = glob.glob('*.txt')
barcodes.sort()
c.barcode_files = barcodes
os.chdir(starting_dir)

# Set barcode file mode  
c.barcode_files_setup = 'individual'

if c.barcode_files_setup == 'individual':
    # One barcode file per input file with matching names
    # Input check
    barnames = [b.split('.')[0] for b in c.barcodes]
    filenames = [f.split('.')[0] for f in c.filenames]
    for fname in filenames:
        if fname not in barnames:
            raise Exception('Set to individual barcode files, yet at least one input'
            'file name does not match the given barcode file names')

# MIDtags
c.cutsite = 'TGCAGG'
c.max_edit_dist = 1
        
# FILTERING
#----------
# Whether to log reads that fail the filtering         
c.log_fails = True       
# Output directory for filtered reads       
c.outdir = 'filtered_data'      
        
# Define Class
Experiment = Workflow(c) 

#===============================================================================
# Setup and run filter
#===============================================================================
Experiment.filter_functions = [setup_propN_filter(0.1),
                               setup_phred_filter(25),
                               setup_cutsite_filter('TCGAGG', 2),
                               setup_overhang_filter('TCGAGG', 'GG', 0)]

Experiment.filter_reads_pipeline()

#===============================================================================
# Process and Correct MID tag 
#===============================================================================
Experiment.process_MIDtag(max_edit_dist = 1, outfile_postfix='-clean')

#===============================================================================
# Cluster Data 
#===============================================================================
allreads_file = 'lane' + LANE + 'allreads-clean.fasta'
trim_reads(infiles=cleaned_files, inpath=cleaned_inpath, 
            outpath=cleaned_inpath, outfile=allreads_file, n=1)
# Variables 
c_thresh = 0.9
n_filter = 8

clustered_file = 'lane' + LANE + 'clustered_reads'
cluster_cdhit(infile=allreads_file, outfile=clustered_file,
              c_thresh=c_thresh, n_filter=n_filter)
#
## Display Summary
#summary(clustered_file)
