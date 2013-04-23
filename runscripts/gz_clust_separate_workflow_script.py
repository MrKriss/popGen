'''
Created on 10 Dec 2012

@author: musselle
'''

from utils.workflow import Workflow

'''----------------------------------------------------------------------------
               Workflow script for Gazelles-Zebras RAD-data 
-----------------------------------------------------------------------------'''

testing = False

experiment_name = 'gz_allg_allz'
description = 'Cluster all gazelles and zebras separately at 90% and 100%'

# Preprocessing 
#===============================================================================
W = Workflow()

W.create_new(name='zebras-gazelles-test', testing=testing)

# Make sure barcodes and datafiles are in their right folders

if testing:
    W.add_datafiles(data_files='testset_500.fastq.bgzf' , barcode_files='*[8].txt')
else:
    W.add_datafiles(data_files='lane6*bgzf' , barcode_files='*[6].txt')
    W.add_datafiles(data_files='lane8*bgzf' , barcode_files='*[8].txt')

# Set parameters
p = {'filtering' : {'propN': 0.1,
                    'phred': 25,
                    'cutsite_edit_dist' : 2,
                    'overhang_edit_dist' : 0},
     'cleaning' : {'max_edit_dist' : 1 }}

if testing:
    W.setup_preprocessing(infiles_pattern="testset_500.fastq.bgzf" , params=p)
else:
    W.setup_preprocessing(infiles_pattern="lane*.bgzf" , params=p)

W.run_preprocessing()

# Clustering 
#===============================================================================
if testing:
    W.add_experiment_name('gz-subgroups-test', 'Test for splitting files by subgroups')
else:
    W.add_experiment_name(experiment_name, description)

default_params = { 'c_thresh' : 0.90,
                   'n_filter' : 8,
                    'maskN' : False}

subgroups = { 'zebra'  : '.*zebra.*',
            'gazelle' : '.*gazelle.*'}

if testing:
    W.setup_clustering(mode='split_by_subgroups', infiles_pattern='test*-clean.fastq.bgzf',
                     default_params=default_params, subgroups=subgroups) 
#     W.setup_clustering(mode='split_by_tags', infiles_pattern='test*-clean.fastq.bgzf',
#                      default_params=default_params) 
else:
    W.setup_clustering(mode='split_by_subgroups', infiles_pattern='lane*-clean.fastq.bgzf',
                     default_params=default_params, subgroups=subgroups) 
#     W.setup_clustering(mode='split_by_tags', infiles_pattern='lane*-clean.fastq.bgzf',
#                      default_params=default_params) 

# Varibles to change, 1 dictionary per run
run_parameters = [ 
                    { 'c_thresh' : 1.0},
                    { 'c_thresh' : 0.90},
                   ]

W.run_clustering(run_parameters, threads=5)




