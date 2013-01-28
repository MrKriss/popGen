'''
Created on 22 Jan 2013

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
''' RUNS SCRIPT FOR Lane 6 '''
#===============================================================================

LANE = '6'

starting_dir = os.getcwd()

# Set paths and file patterns 
inpath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
#barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
#inpath = '/home/musselle/san/data/lane' + LANE
#barpath = '/home/musselle/san/data/barcodes'

os.chdir(inpath)
raw_files = glob.glob('*[0-9].fastq.bgzf')
raw_files.sort()
os.chdir(starting_dir)

# Setup Filters

filter_functions = [setup_illumina_filter(), 
                    setup_propN_filter(0.1),
                    setup_phred_filter(25),
                    setup_cutsite_filter('TCGAGG', 2),
                    setup_overhang_filter('TCGAGG', 'GG', 0)]

filter_reads_pipeline(infiles=raw_files, inpath=inpath, filterfuncs=filter_functions, 
                          outdir='filtered_reads', log_fails=True)
