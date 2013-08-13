'''
Created on 10 Dec 2012

@author: musselle
'''

from utils.workflow import Workflow

'''----------------------------------------------------------------------------
   Workflow script Template for adding experiment to Gazelles-Zebras RAD-data 
-----------------------------------------------------------------------------'''

experiment_name = ''
description = ''

dataset_name = ''

# Load previously processed data info
#===============================================================================
W = Workflow()
W.load(name=dataset_name)

# Clustering 
#===============================================================================

W.add_experiment_name(experiment_name, description)

default_params = { 'c_thresh' : 0.90,
                   'n_filter' : 8,
                    'maskN' : False}

#===============================================================================
# Choose method to split files 
#===============================================================================

subgroups = { 'zebra'  : '.*zebra.*',
            'gazelle' : '.*gazelle.*'}

splitby = None

if splitby == 'subgroups':
    W.setup_clustering(mode='split_by_subgroups', infiles_pattern='lane*-clean.fastq.bgzf',
                     default_params=default_params, subgroups=subgroups) 
elif splitby == 'tags':
    W.setup_clustering(mode='split_by_tags', infiles_pattern='lane*-clean.fastq.bgzf',
                     default_params=default_params) 
#===============================================================================

# Varibles to change, 1 dictionary per run
run_parameters = [ 
                    { 'c_thresh' : 1.0},
                    { 'c_thresh' : 0.90},
                   ]

W.run_clustering(run_parameters, threads=5)

W.cleanup_files('fasta')


