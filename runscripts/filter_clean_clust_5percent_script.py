'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
from os.path import join as joinp
import sys 

import glob
import socket
import re

from utils import get_data_prefix
from preprocess import  Preprocessor, ConfigClass
from cluster import Clustering

from database import PopGen_DB

#==============================================================================
''' Filter/Clean/Cluster SCRIPT FOR ALLL READS IN Gazelles-Zebras RAD data Lane 6'''
#===============================================================================

experiment_name = 'gz' 

#===============================================================================
# Setup Configuration
#===============================================================================
starting_dir = os.getcwd()

c = ConfigClass()

c.experiment_name = experiment_name

# Work out where data is stored on this machine
prefix = get_data_prefix()

# Set paths 
#c.data_inpath =  joinp(prefix,'gazelles-zebras') 
c.data_inpath =  joinp(prefix,'gazelles-zebras', 'testset') 
c.barcode_inpath = joinp(prefix,'gazelles-zebras', 'barcodes')
c.filtered_outpath = joinp(prefix,'gazelles-zebras', 'processed-data')
c.tag_processed_outpath = joinp(prefix,'gazelles-zebras', 'processed-data')
c.tag_split_outpath = joinp(prefix,'gazelles-zebras', 'processed-data', 'per-individual')
c.clusters_outpath = joinp(prefix,'gazelles-zebras', 'clusters')

# Setup input files and barcodes
os.chdir(c.data_inpath)
#raw_files = glob.glob('*[0-9].fastq.bgzf')
raw_files = glob.glob('testset_1m.fastq.bgzf')
assert raw_files
raw_files.sort()
c.raw_input_files = raw_files 

os.chdir(c.barcode_inpath)
barcodes = glob.glob('*[68].txt')
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
c.log_fails = False
       
#===============================================================================
# Make Database for samples
#===============================================================================

db_path = joinp(prefix,'gazelles-zebras') 
db = PopGen_DB(joinp(db_path, 'gz_samples.db'), recbyname=True)
        
#L6_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[6].txt')) 
L8_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[8].txt')) 

#r1 = re.compile('lane6*.bgzf')
#r2 = re.compile('lane8*.bgzf')


#L6_datafiles = filter(r1.match, c.raw_input_files)
#L8_datafiles = filter(r2.match, c.raw_input_files)


# Associate subsets of the data files list to their respective barcode files.   
#db.add_barcodes_datafiles(L6_barcode_files, L6_datafiles)
#db.add_barcodes_datafiles(L8_barcode_files, L8_datafiles)

# Testing 
r3 = re.compile('testset_1m.fastq.bgzf')
datafiles = filter(r3.match, c.raw_input_files)
db.add_barcodes_datafiles(L8_barcode_files, datafiles)
#db.add_barcodes_datafiles(L6_barcode_files, datafiles)

# Define Preprocessing Class
Preprocess = Preprocessor(c, db) 

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
Preprocess.cleanup_files('filtered') # Remove filtered intermediate files 

#===============================================================================
# Split Data into individuals based on MID tag 
#===============================================================================
Preprocess.split_by_tags()
Preprocess.cleanup_files('tag_processed') # Remove MID tag processed intermediate files 

#===============================================================================
# Prepare for Clustering 
#===============================================================================
files2cluster, path = Preprocess.trim_reads(mode='separate', n=1)

#===============================================================================
# Cluster Data 
#===============================================================================

# default Vars for clustering 
#default_vars = { 'c_thresh' : 0.90,
#                 'n_filter' : 8,
#                 'threads' : 1,
#                 '
#                 'maskN' : False}

# Variations to run
batch_parameters = [ 
                    { 'c_thresh' : 1.0},
                    { 'c_thresh' : 0.90},
                   ]
                   
Clusterer = Clustering(c, infiles=files2cluster, inpath=path) 

Clusterer.run_batch_cdhit_clustering(batch_parameters, threads=1)

## Display Summary
#summary(clustered_file)