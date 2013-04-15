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
from cluster import ClusterClass

from database import PopGen_DB

starting_dir = os.getcwd()
c = ConfigClass()

#==============================================================================
''' Filter/Clean/Cluster SCRIPT FOR ALLL READS IN Gazelles-Zebras RAD data Lane 6'''
#===============================================================================

# These should be unique for each experiment, else results table is overwritten
c.experiment_name = 'gz-allg-allz'
c.experiment_description = '''Clustering all gazelles and all zebras separately'''

c.db_name = 'gz_allg-allz.db'

# Testing
testing = False 
if testing: 
    #testfile = 'testset_10m.fastq.bgzf'
    testfile = 'testset_500.fastq.bgzf'

#===============================================================================
# Setup Configuration
#===============================================================================

# Work out where data is stored on this machine
prefix = get_data_prefix()

# Set paths 
c.data_inpath =  joinp(prefix,'gazelles-zebras', 'raw-data') 
if testing:
    c.data_inpath =  joinp(prefix,'gazelles-zebras', 'testset')

c.barcode_inpath = joinp(prefix,'gazelles-zebras', 'barcodes')
c.filtered_outpath = joinp(prefix,'gazelles-zebras', 'processed-data')
c.tag_processed_outpath = joinp(prefix,'gazelles-zebras', 'processed-data')
c.tag_split_outpath = joinp(prefix,'gazelles-zebras', 'processed-data', 'per-species')
c.clusters_outpath = joinp(prefix,'gazelles-zebras', 'clusters')

# Setup input files and barcodes
os.chdir(c.data_inpath)
raw_files = glob.glob('*[0-9].fastq.bgzf')
if testing:
    raw_files = glob.glob(testfile)
assert raw_files
raw_files.sort()
c.raw_input_files = raw_files 

os.chdir(c.barcode_inpath)
barcodes = glob.glob('*[68].txt')
barcodes.sort()
c.barcode_files = barcodes
os.chdir(starting_dir)

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
# Set Parameters
#===============================================================================
# Define Preprocessing Class
Preprocess = Preprocessor(c) 

# For Filtering
#------------------------------------------------------------------------------ 
filter_params = {'propN': 0.1,
                 'phred': 25,
                 'cutsite_edit_dist' : 2,
                 'overhang_edit_dist' : 0,
                 'overhang_target': 'GG'}
c.filter_funtion_params = str(filter_params)

Preprocess.filter_functions = [
                Preprocess.make_propN_filter(filter_params['propN']),
                Preprocess.make_phred_filter(filter_params['phred']),
                Preprocess.make_cutsite_filter(max_edit_dist=filter_params['cutsite_edit_dist']),
                Preprocess.make_overhang_filter('TCGAGG', 'GG', filter_params['overhang_edit_dist'])
                ]

# For Clustering 
#------------------------------------------------------------------------------ 
# Default vars for clustering 
default_vars = { 'c_thresh' : 0.90,
                 'n_filter' : 8,
                 'maskN' : False}

# Varibles to change, 1 dictionary per run
run_parameters = [ 
                    { 'c_thresh' : 1.0},
                    { 'c_thresh' : 0.90},
                   ]

#===============================================================================
# Make/Update Database for Experiment
#===============================================================================
db_path = joinp(prefix,'gazelles-zebras') 
db = PopGen_DB(joinp(db_path, c.db_name), recbyname=True, new=True)
Preprocess.db = db # Pass database reference to Preprocessor Object
        
# Testing 
if testing:
    L8_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[8].txt')) 
    r3 = re.compile(testfile)
    datafiles = filter(r3.match, c.raw_input_files)
    db.add_barcodes_datafiles(L8_barcode_files, datafiles, datafile_type='raw_mixed')
else:
    L6_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[6].txt')) 
    L8_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[8].txt')) 
    
    r1 = re.compile('lane6*.bgzf')
    r2 = re.compile('lane8*.bgzf')
    
    L6_datafiles = filter(r1.match, c.raw_input_files)
    L8_datafiles = filter(r2.match, c.raw_input_files)
    
    # Associate subsets of the data files list to their respective barcode files.   
    db.add_barcodes_datafiles(L6_barcode_files, L6_datafiles)
    db.add_barcodes_datafiles(L8_barcode_files, L8_datafiles)


# Add Experimental details and config object in database
#------------------------------------------------------------------------------ 
c.exp_id = db.add_experiment(config=c, exp_type='clustering')

#===============================================================================
# Run Filtering
#===============================================================================
#Preprocess.filter_reads_pipeline()

#===============================================================================
# Process and Correct MID tag 
#===============================================================================
#Preprocess.process_MIDtag(max_edit_dist = 1, outfile_postfix='-clean')
#Preprocess.cleanup_files('filtered') # Remove filtered intermediate files 

#===============================================================================
# Splitting by species zebras and gazelles
#===============================================================================
Preprocess.set_input_files(data_inpath=c.tag_processed_outpath, file_pattern="*.fastq.bgzf")


subgroups = { 'zebra'  : '.*zebra.*',
            'gazelle' : '.*gazelle.*'}

Preprocess.split_by_subgroups(subgroups)
#Preprocess.cleanup_files('tag_processed') # Remove MID tag processed intermediate files 

#===============================================================================
# Prepare for Clustering 
#===============================================================================

# Separate all Barcodes into separate files for clustering individually
files2cluster, path = Preprocess.trim_reads(mode='separate', n=1)

#===============================================================================
# Cluster Data 
#===============================================================================
                   
Cluster = ClusterClass(infiles=files2cluster, inpath=path, defaults=default_vars) 
Cluster.c = c
Cluster.db = db

out_list = Cluster.run_batch_cdhit_clustering(run_parameters, threads=1)

## Display Summary
