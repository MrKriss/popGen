'''
Created on 10 Dec 2012

@author: musselle
'''
import os 
from os.path import join as joinp
import sys 

import glob
import re

from utils import get_data_prefix
from utils.preprocess import  Preprocessor, ConfigClass
from utils.cluster import ClusterClass
import cPickle as pkl

from file_conversions import makeSQLindex

from utils.database import Popgen_db

class Workflow(object):
    ''' Container for all preprocessing, filtering, and exploritry analysis with 
    a particular dataset
    '''
    
    def __init__(self):
        self.c = ConfigClass()
    
    def create_new(self, name=None, db_name=None, testing=False):
        ''' Setup directory structure and initialise config file and database 
        for the given dataset name.'''
        
        if (name is None) or (type(name) is not str):
            raise Exception('Must specify a valid name for the dataset.')
        if db_name is None:
            db_name = name + '.db'
        
        # Setup Configuration
        prefix =  get_data_prefix()
        
        # Default path locations
        self.c.testing = testing
        self.c.root_name = name
        self.c.db_name = db_name
        if testing:
            self.c.data_inpath =  joinp(prefix,name, 'testset')
        else:
            self.c.data_inpath =  joinp(prefix, name, 'raw-data') 
        self.c.barcode_inpath = joinp(prefix, name , 'barcodes')
        self.c.filtered_outpath = joinp(prefix, name , 'processed-data')
        self.c.tag_processed_outpath = joinp(prefix, name, 'processed-data')
        self.c.tag_splitby_sample_outpath = joinp(prefix, name, 'processed-data', 'per-sample')
        self.c.tag_splitby_subgroup_outpath = joinp(prefix, name, 'processed-data', 'per-subgroup')
        self.c.clusters_outpath = joinp(prefix, name, 'clusters')
        self.c.db_path = joinp(prefix,  name)
        self.c.cdhit_path = os.path.expanduser("~/bin/cd-hit-v4.6.1")


        # Create directories of they dont exist
        for attr in dir(self.c): 
            if 'path' in attr:
                path = getattr(self.c, attr)
                if not os.path.exists(path):
                    os.makedirs(path)
        
        # Set interim file suffixes
        self.c.filtered_files_postfix = '-pass'
        self.c.tag_processed_files_postfix = '-clean'

        # MIDtags
        self.c.cutsite = 'TGCAGG'
        self.c.max_edit_dist = 2
        
        # FILTERING
        # Whether to log reads that fail the filtering         
        self.c.log_fails = False
        
        # Create new Database 
        self.db = Popgen_db(joinp(self.c.db_path, db_name), recbyname=True, new=True)
        
        # Save config
        f = open(joinp(self.c.db_path, '.' + name + '-config.pkl'), 'w')
        pkl.dump(self.c, f)
        
    def load(self, name=None, db_name=None, recbyname=True):
        ''' Load a pre-existing directory structure, config file and database 
        with the given dataset name.'''
        
        if (name is None) or (type(name) is not str):
            raise Exception('Must specify a valid name for the dataset.')
        if db_name is None:
            db_name = name + '.db'
        
        # Load config
        prefix = get_data_prefix()
        path2config = joinp(prefix, name, '.' + name + '-config.pkl')
        path2db = joinp(prefix, name, db_name)
        
        self.c = pkl.load(open(path2config))
        self.db = Popgen_db(path2db, recbyname=recbyname)
        
    def add_datafiles(self, data_files=None , barcode_files=None ):
        ''' Add datafiles and barcodes in pairs to the database.
         
        Each pair defines the samples present in the datafiles listed. 
        
        If 'files' or 'barcodes' is a str, it is interpreted as a glob to the 
        data_path / barcode_path respecively.
        '''
        
        if type(data_files) is str:
            data_files = glob.glob(joinp(self.c.data_inpath, data_files))
        if type(barcode_files) is str:
            barcode_files = glob.glob(joinp(self.c.barcode_inpath, barcode_files))
            
        # Input samples-datafiles info in Database 
        self.db.add_barcodes_datafiles(barcode_files, data_files, datafile_type='raw_mixed') 
                
    def setup_preprocessing(self, infiles_pattern, params=None,  param_id=None):
        ''' Setup the preprocessing function for the workflow '''
        
        # Get params if id given 
        if params is None and param_id is None:
            raise Exception("No parameters and no parameter id to lookup.")
        if param_id:
            params = self.db.get_binary('params', 'filtering_parameterID', param_id, table='filtering_parameters')
            assert params, "No data returned from database for param_id: %s" % param_id
            self.c.filterparam_id = param_id
        else:
            # Insert parameters dictionary into filter_parameters table
            self.c.filterparam_id = self.db.insert_binary(params, col='params', table='filtering_parameters')
            
        # Define Preprocessing Class and set inputs
        self.Preprocessor = Preprocessor(self.c)
        self.Preprocessor.db = self.db # Pass database reference to Preprocessor 
        
        self.Preprocessor.set_input_files(data_files=infiles_pattern, data_inpath=self.c.data_inpath)
         
        self.Preprocessor.filter_functions = [
                self.Preprocessor.make_propN_filter(params['filtering']['propN']),
                self.Preprocessor.make_phred_filter(params['filtering']['phred']),
                self.Preprocessor.make_cutsite_filter(max_edit_dist=params['filtering']['cutsite_edit_dist']),
                self.Preprocessor.make_overhang_filter('TCGAGG', 'GG', params['filtering']['overhang_edit_dist'])
                ]
        
        # Save addition to config file
        path = joinp(self.c.db_path, '.' + self.c.root_name + '-config.pkl')
        if os.path.exists(path):
            os.remove(path)
            pkl.dump(self.c, open(path, 'w'))
            
    def run_preprocessing(self):
        ''' Call the Preprocessing functions for the workflow '''

        params = self.db.get_binary('params', 'filtering_parameterID', self.c.filterparam_id, table='filtering_parameters')
            
        # Process and Correct MID tag 
        self.Preprocessor.filter_reads_pipeline()
        self.Preprocessor.process_MIDtag(max_edit_dist = params['cleaning']['max_edit_dist'])

        # Remove filtered intermediate files 
        self.Preprocessor.cleanup_files('filtered') 

    def setup_clustering(self, mode, infiles_pattern, default_params=None, subgroups=None):
        ''' Setup files for the Clustering function of the workflow. 
        
        Does the necessary splitting and trimming of files if specified in mode.
        '''
        
        # Input Checks
        if not hasattr(self, 'c'):
            raise Exception('Must first load a config file')
        if not hasattr(self, 'Preprocessor'):
            self.Preprocessor = Preprocessor(self.c)
            self.Preprocessor.db = self.db

        # Set files to process for clustering         
        self.Preprocessor.set_input_files(data_files=infiles_pattern, data_inpath=self.c.tag_processed_outpath)
        
        if mode == 'split_by_tags':
            (outfiles, outpath) = self.Preprocessor.split_by_tags()
            
            # Create index for files clustered
            makeSQLindex(outfiles, outpath)
            
            files2cluster, path = self.Preprocessor.trim_reads(mode='separate', 
                                outpath=self.c.tag_splitby_sample_outpath, n=1)
        elif mode == 'split_by_subgroups':
            if subgroups is None:
                raise Exception("No subgroups specified")
            (outfiles, outpath) = self.Preprocessor.split_by_subgroups(subgroups)
            
            # Create index for files clustered
            makeSQLindex(outfiles, outpath)
            
            files2cluster, path = self.Preprocessor.trim_reads(mode='separate',
                             outpath=self.c.tag_splitby_subgroup_outpath,  n=1)
        elif mode == 'no_split':
            # Create index for files clustered
            makeSQLindex(filepattern=infiles_pattern, data_inpath=self.c.tag_processed_outpath)
            
            files2cluster, path = self.Preprocessor.trim_reads(mode='grouped', n=1)
        else:
            raise Exception(' No valid mode specified. ')

        # Setup Clusterer
        self.Clusterer = ClusterClass(infiles=files2cluster, inpath=path, 
                                      config=self.c, db=self.db, defaults=default_params)

    def run_clustering(self, run_parameters, **kwargs):
        ''' Run CDHIT using the specified parameters. passes on any kwargs'''
        
        outputs_list = self.Clusterer.run_batch_cdhit_clustering(run_parameters, **kwargs)
        outnames_list, out_path, counters_list = zip(*outputs_list)
        
        return outnames_list, out_path, counters_list
        
        
    def add_experiment_name(self, name, description):
        ''' Add Experimental details and config object in database'''

        # These should be unique for each experiment, else results table is overwritten
        self.c.experiment_name = name
        self.c.experiment_description = description
        self.c.exp_id = self.db.add_experiment(config=self.c, exp_type='clustering')
        
    def cleanup_files(self, file_type):
        ''' Remove all intermediate files specified '''
        self.Preprocessor.cleanup_files(file_type)
        
        
        

