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
import cPickle as pkl

from database2 import Popgen_db

starting_dir = os.getcwd()
c = ConfigClass()

# Work out where data is stored on this machine
prefix = get_data_prefix()

#==============================================================================
''' Filter and Clean SCRIPT FOR ALLL READS IN Gazelles-Zebras RAD-data'''
#===============================================================================

# These should be unique for each experiment, else results table is overwritten
#c.experiment_name = 'gz-allg-allz'
#c.experiment_description = '''Clustering all gazelles and all zebras separately'''

c.root = 'gazelles-zebras'

c.db_name = 'gz_allg-allz.db'

db_path = joinp(prefix, c.root) 

db = Popgen_db(joinp(db_path, c.db_name), recbyname=True, new=True)

# Testing
testing = False 
if testing: 
    #testfile = 'testset_10m.fastq.bgzf'
    testfile = 'testset_500.fastq.bgzf'

#===============================================================================
# Setup Configuration
#===============================================================================

# Set paths 
if testing:
    c.data_inpath =  joinp(prefix, c.root, 'testset')
else:
    c.data_inpath =  joinp(prefix, c.root, 'raw-data') 
c.barcode_inpath = joinp(prefix, c.root , 'barcodes')
c.filtered_outpath = joinp(prefix, c.root , 'processed-data')
c.tag_processed_outpath = joinp(prefix, c.root, 'processed-data')
c.tag_split_outpath = joinp(prefix, c.root, 'processed-data', 'per-species')
c.clusters_outpath = joinp(prefix, c.root, 'clusters')

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
# Update/input samples-datafiles info in Database 
#===============================================================================
        
# Testing 
if testing:
    L8_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[8].txt')) 
    r3 = re.compile(testfile)
    datafiles = filter(r3.match, c.raw_input_files)
    db.add_barcodes_datafiles(L8_barcode_files, datafiles, datafile_type='raw_mixed')
else:

    L6_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[6].txt')) 
    L8_barcode_files = glob.glob(joinp(c.barcode_inpath, '*[8].txt')) 
    L6_datafiles = glob.glob(joinp(c.data_inpath, 'lane6*bgzf'))
    L8_datafiles = glob.glob(joinp(c.data_inpath, 'lane8*bgzf'))

    # Associate subsets of the data files list to their respective barcode files.   
    db.add_barcodes_datafiles(L6_barcode_files, L6_datafiles, datafile_type='raw_mixed')
    db.add_barcodes_datafiles(L8_barcode_files, L8_datafiles, datafile_type='raw_mixed')

# Define Preprocessing Class
Preprocess = Preprocessor(c) 
Preprocess.db = db # Pass database reference to Preprocessor Object

## Add Experimental details and config object in database
##------------------------------------------------------------------------------ 
#c.exp_id = db.add_experiment(config=c, exp_type='clustering')

#===============================================================================
# Setup Filtering Parameters
#===============================================================================
p = {'filtering' : {'propN': 0.1,
                    'phred': 25,
                    'cutsite_edit_dist' : 2,
                    'overhang_edit_dist' : 0},
     'cleaning' : {'max_edit_dist' : 1 }}

# Insert into filter_parameters table
c.filterparam_id = db.insert_binary(p, col='params', table='filtering_parameters')

Preprocess.filter_functions = [
                Preprocess.make_propN_filter(p['filtering']['propN']),
                Preprocess.make_phred_filter(p['filtering']['phred']),
                Preprocess.make_cutsite_filter(max_edit_dist=p['filtering']['cutsite_edit_dist']),
                Preprocess.make_overhang_filter('TCGAGG', 'GG', p['filtering']['overhang_edit_dist'])
                ]

#===============================================================================
# Run Filtering
#===============================================================================
Preprocess.filter_reads_pipeline()

#===============================================================================
# Process and Correct MID tag 
#===============================================================================
Preprocess.process_MIDtag(max_edit_dist = p['cleaning']['max_edit_dist'])
Preprocess.cleanup_files('filtered') # Remove filtered intermediate files 

# Store or pass on config file to clustering section
# Pickle config 
pkl.dump(c, open('config.pkl', 'wb'))
