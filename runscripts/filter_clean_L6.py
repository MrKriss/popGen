'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import glob

from utils.filter import setup_filter, filter_reads
from utils.file_conversions import process_MIDtag

#==============================================================================
''' RUNS SCRIPT FOR ALLL READS IN LANE 6 '''
#===============================================================================

#starting_dir = os.getcwd()
#
## Set paths and file patterns 
#datapath = '/space/musselle/datasets/gazellesAndZebras/lane6'
#barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
#os.chdir(datapath)
#raw_files = glob.glob('*[0-9].fastq.bgzf')
#raw_files.sort()
#
##===============================================================================
## Setup and run filter
##===============================================================================
#f = setup_filter({'phred': 20, 'propN': 0.10})
#outdir = 'L6_phredprop_filtered'
#filter_reads(infiles=raw_files, filepattern=True, 
#             datapath=datapath, outdir=outdir, filterfunc=f)

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
process_MIDtag(infiles=raw_files, barcodes='*[6].txt', barcode_pattern=True,
               datapath=filtered_datapath, barcode_path=barpath, 
               outfile_postfix=cleaned_file_postfix, outdir=cleaned_outdir)
# Update names and path
cleaned_files = []
for name in raw_files:
    temp = name.split('.')
    temp[0] = temp[0] + cleaned_file_postfix
    temp = '.'.join(temp) 
    cleaned_files.append(temp) 
cleaned_datapath = filtered_datapath + '/' + cleaned_outdir

#===============================================================================
# Cluster Data 
#===============================================================================









