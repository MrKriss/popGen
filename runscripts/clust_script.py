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

#==============================================================================
''' Cluster SCRIPT FOR ALLL READS IN Gazelles-Zebras RAD data '''
#===============================================================================

testing = True

# Load Config and setup

# Work out where data is stored on this machine
prefix = get_data_prefix()
path2config = joinp(prefix, 'gazelles-zebras', 'config.pkl')

c = pkl.load(open(path2config, 'rb'))
Preprocess = Preprocessor(c)

# Set next files to process 
if testing:
    Preprocess.set_input_files(data_inpath=c.tag_processed_outpath, file_pattern='testset_500-pass-clean.fastq.bgzf')
else:  
    Preprocess.set_input_files(data_inpath=c.tag_processed_outpath, file_pattern='lane[68]*-clean.bgzf')

# Connect to database 
db_path = joinp(prefix, c.root) 
db = Popgen_db(joinp(db_path, c.db_name), recbyname=True)
Preprocess.db = db

#==============================================================================
# Add Experimental details and config object in database
#==============================================================================

# These should be unique for each experiment, else results table is overwritten
c.experiment_name = 'gz-allg-allz'
c.experiment_description = '''Clustering all gazelles and all zebras separately'''
c.exp_id = db.add_experiment(config=c, exp_type='clustering')

#===============================================================================
# Splitting by species zebras and gazelles
#===============================================================================
subgroups = { 'zebra'  : '.*zebra.*',
            'gazelle' : '.*gazelle.*'}

Preprocess.split_by_subgroups(subgroups)

# Preprocess.cleanup_files('tag_processed') # Remove MID tag processed intermediate files 

#===============================================================================
# Prepare for Clustering 
#===============================================================================

# Separate all reads into separate files for each group so they are clustered individually
files2cluster, path = Preprocess.trim_reads(mode='separate', n=1)


#===============================================================================
# Cluster Data 
#===============================================================================

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
                   
Cluster = ClusterClass(infiles=files2cluster, inpath=path, defaults=default_vars) 
Cluster.c = c
Cluster.db = db

out_list = Cluster.run_batch_cdhit_clustering(run_parameters, threads=1)

## Display Summary

























# Work out where data is stored on this machine
prefix = get_data_prefix()

# This should be unique for each experiment
c.experiment_name = 'gz-allg-allz'
c.experiment_description = '''Clustering all gazelles and all zebras separately'''

c.db_name = 'gz_allg-allz.db'

db_path = joinp(prefix,'gazelles-zebras') 
db = PopGen_DB(joinp(db_path, c.db_name), recbyname=True)


#===============================================================================
# Setup Configuration
#===============================================================================

# Update paths
c.data_inpath = joinp(prefix,'gazelles-zebras', 'raw-data')
c.barcode_inpath = joinp(prefix,'gazelles-zebras', 'barcodes')
c.filtered_outpath = joinp(prefix,'gazelles-zebras', 'processed-data')
c.tag_processed_outpath = joinp(prefix,'gazelles-zebras', 'processed-data')
c.tag_split_outpath = joinp(prefix,'gazelles-zebras', 'processed-data', 'per-species')
c.clusters_outpath = joinp(prefix,'gazelles-zebras', 'clusters')

# Set interim file suffixes
c.filtered_files_postfix = '-pass'
c.tag_processed_files_postfix = '-clean'

# MIDtags
c.cutsite = 'TGCAGG'
c.max_edit_dist = 2

# Get previous config file 
c2 =  db.get_binary('config', 'name', c.experiment_name, 'experiments')
c.filter_funtion_params = c2.filter_funtion_params


#===============================================================================
# Set Parameters
#===============================================================================
# Define Preprocessing Class
Preprocess = Preprocessor(c) 

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
Preprocess.db = db # Pass database reference to Preprocessor Object
        


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
