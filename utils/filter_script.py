'''
Created on 22 Jan 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import glob

from preprocess import setup_filter, filter_reads, process_MIDtag, reads2fasta
from cluster import cluster_cdhit, summary

#==============================================================================
''' RUNS SCRIPT FOR ALLL READS IN LANE 6 '''
#===============================================================================

LANE = '6'

starting_dir = os.getcwd()

# Set paths and file patterns 
#inpath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
#barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
inpath = '/home/musselle/san/data/lane' + LANE
barpath = '/home/musselle/san/data/barcodes'
os.chdir(inpath)
raw_files = glob.glob('*[0-9].fastq.bgzf')
raw_files.sort()




















if __name__ == '__main__':
    
    
    