'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import glob

from preprocess import filter_reads, process_MIDtag, reads2fasta, setup_cutsite_filter, \
setup_overhang_filter, setup_phred_filter, setup_propN_filter, setup_illumina_filter, \
filter_reads_pipeline

from cluster import cluster_cdhit, summary

#==============================================================================
''' RUNS SCRIPT FOR ALLL READS IN LANE 8 '''
#===============================================================================

LANE = '8'

starting_dir = os.getcwd()

# Set paths and file patterns 
#inpath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
#barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
inpath = '/home/musselle/san/data/lane' + LANE
barpath = '/home/musselle/san/data/barcodes'
os.chdir(inpath)
raw_files = glob.glob('*[0-9].fastq.bgzf')
raw_files.sort()

#===============================================================================
# Setup and run filter
#===============================================================================
filter_functions = [setup_illumina_filter(), 
                    setup_propN_filter(0.1),
                    setup_phred_filter(25),
                    setup_cutsite_filter('TCGAGG', 2),
                    setup_overhang_filter('TCGAGG', 'GG', 0)]

outdir = 'L' + LANE + '_filtered'

filter_reads_pipeline(infiles=raw_files, inpath=inpath, filterfuncs=filter_functions, 
                          outdir=outdir, log_fails=True)

# Update names and path
filtered_files = []
for name in raw_files:
    temp = name.split('.')
    temp[0] = temp[0] + '-pass'
    temp = '.'.join(temp) 
    filtered_files.append(temp)
filtered_inpath = os.path.join(inpath, outdir)

#===============================================================================
# Process and Correct MID tag 
#===============================================================================
cleaned_file_postfix = '-clean' 
cleaned_outdir = '' #'cleaned_data'
barcode_pattern = '*[' + LANE + '].txt'

process_MIDtag(infiles=filtered_files, barcodes =barcode_pattern,
                   barcode_pattern=True, 
                   inpath=filtered_inpath, barcode_path=barpath,
                   outfile_postfix=cleaned_file_postfix, outdir=cleaned_outdir, 
                   MIDtag_len = 6, max_edit_dist = 1, cutsite_len = 6)

# Update names and path
cleaned_files = []
for name in filtered_files:
    temp = name.split('.')
    temp[0] = temp[0] + cleaned_file_postfix
    temp = '.'.join(temp) 
    cleaned_files.append(temp) 
cleaned_inpath = os.path.join(filtered_inpath, cleaned_outdir)

#===============================================================================
# Cluster Data 
#===============================================================================
allreads_file = 'lane' + LANE + 'allreads-clean.fasta'
reads2fasta(infiles=cleaned_files, inpath=cleaned_inpath, 
            outpath=cleaned_inpath, outfile=allreads_file)
# Variables 
c_thresh = 0.9
n_filter = 8

clustered_file = 'lane' + LANE + 'clustered_reads'
#cluster_cdhit(infile=allreads_file, outfile=clustered_file,
#              c_thresh=c_thresh, n_filter=n_filter)
#
## Display Summary
#summary(clustered_file)
