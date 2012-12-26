'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import glob

from filter import setup_filter, filter_reads
from file_conversions import process_MIDtag, reads2fasta
from cluster import cluster_cdhit, summary

#==============================================================================
''' RUNS SCRIPT FOR ALLL READS IN LANE 6 '''
#===============================================================================

LANE = '6'

starting_dir = os.getcwd()

# Set paths and file patterns 
datapath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
os.chdir(datapath)
raw_files = glob.glob('*[0-9].fastq.bgzf')
raw_files.sort()

#===============================================================================
# Setup and run filter
#===============================================================================
f = setup_filter({'phred': 20, 'propN': 0.10})
outdir = 'L6_phredprop_filtered'
filter_reads(infiles=raw_files, filepattern=True, 
             datapath=datapath, outdir=outdir, filterfunc=f)

# Update names and path
filtered_files = []
for name in raw_files:
    temp = name.split('.')
    temp[0] = temp[0] + '-pass'
    temp = '.'.join(temp) 
    filtered_files.append(temp)
filtered_datapath = datapath + '/' + outdir

#===============================================================================
# Process and Correct MID tag 
#===============================================================================
cleaned_file_postfix = '-clean' 
cleaned_outdir = 'cleaned_data'
barcode_pattern = '*[' + LANE + '].txt'
process_MIDtag(infiles=filtered_files, barcodes=barcode_pattern,
               barcode_pattern=True, datapath=filtered_datapath, 
               barcode_path=barpath, outfile_postfix=cleaned_file_postfix, 
               outdir=cleaned_outdir)
# Update names and path
cleaned_files = []
for name in filtered_files:
    temp = name.split('.')
    temp[0] = temp[0] + cleaned_file_postfix
    temp = '.'.join(temp) 
    cleaned_files.append(temp) 
cleaned_datapath = filtered_datapath + '/' + cleaned_outdir

#===============================================================================
# Cluster Data 
#===============================================================================
allreads_file = 'lane' + LANE + 'allreads-clean.fasta'
reads2fasta(infiles=cleaned_files, inpath=cleaned_datapath, outfile=allreads_file)

# Variables 
c_thresh = 0.9
n_filter = 8

clustered_file = 'lane' + LANE + 'clustered_reads'
#cluster_cdhit(infile=allreads_file, outfile=clustered_file,
#              c_thresh=c_thresh, n_filter=n_filter)
#
## Display Summary
#summary(clustered_file)



