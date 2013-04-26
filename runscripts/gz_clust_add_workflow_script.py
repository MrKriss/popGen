'''
Created on 10 Dec 2012

@author: musselle
'''

from utils.workflow import Workflow

'''----------------------------------------------------------------------------
               Workflow script for Gazelles-Zebras RAD-data 
-----------------------------------------------------------------------------'''

testing = False

experiment_name = 'gz_byMIDs_90_100'
description = 'Cluster individual samples separately by MIDtag at 90% and 100% identity'

# Load previously processed data info
#===============================================================================
W = Workflow()

W.load(name='gazelles-zebras')
# W.create_new(name='gazelles-zebras', testing=testing)

# Clustering 
#===============================================================================
if testing:
    W.add_experiment_name('gz-subgroups-test', 'Test for splitting files by subgroups')
else:
    W.add_experiment_name(experiment_name, description)

default_params = { 'c_thresh' : 0.90,
                   'n_filter' : 8,
                    'maskN' : False}

# subgroups = { 'zebra'  : '.*zebra.*',
#             'gazelle' : '.*gazelle.*'}

if testing:
    W.setup_clustering(mode='split_by_subgroups', infiles_pattern='test*-clean.fastq.bgzf',
                     default_params=default_params, subgroups=subgroups) 
#     W.setup_clustering(mode='split_by_tags', infiles_pattern='test*-clean.fastq.bgzf',
#                      default_params=default_params) 
else:
#     W.setup_clustering(mode='split_by_subgroups', infiles_pattern='lane*-clean.fastq.bgzf',
#                      default_params=default_params, subgroups=subgroups) 
    W.setup_clustering(mode='split_by_tags', infiles_pattern='lane*-clean.fastq.bgzf',
                     default_params=default_params) 

# Varibles to change, 1 dictionary per run
run_parameters = [ 
                    { 'c_thresh' : 1.0},
                    { 'c_thresh' : 0.90},
                   ]

W.run_clustering(run_parameters, threads=5)

W.cleanup_files('fasta')


