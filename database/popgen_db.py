'''
Created on 27 Jun 2013

@author: musselle
'''
import os
import cPickle as pkl
import glob

from utils import get_path_prefix
from general_utilities import set_trace

from release.lib.core_db import SQLdatabase

class Popgen_db(SQLdatabase):  
   
    def __init__(self, db_file="sample.db", recbyname=True, new=False):
       
        create_tables = 0
        if not os.path.exists(db_file) or new:
            create_tables = 1

        Database.__init__(self, db_file, recbyname) 
    
        if create_tables:
            self.create_tables()
       
    def create_tables(self):

        with self.con as con:
    
            # Not sure if i need this 
            con.execute('pragma foreign_keys=OFF')

            curs = con.cursor()        
    
            # Create Tables
            curs.execute('DROP TABLE IF EXISTS samples')
            curs.execute(''' CREATE TABLE samples (
            sampleId INTEGER PRIMARY KEY, 
            MIDtag TEXT, 
            description TEXT UNIQUE, 
            read_count INTEGER,
            read_file TEXT) ''')
            self.tables.append('samples')
            
            curs.execute('DROP TABLE IF EXISTS datafiles')
            curs.execute(''' CREATE TABLE datafiles (
            datafileId INTEGER PRIMARY KEY, 
            filename TEXT UNIQUE,
            type TEXT, 
            filtering_parameterId BLOB )''')
            self.tables.append('datafiles')
            
            curs.execute('DROP TABLE IF EXISTS samples_datafiles')
            curs.execute(''' CREATE TABLE samples_datafiles (
            linkId INTEGER PRIMARY KEY, 
            sampleId INTEGER, 
            datafileId INTEGER, 
            UNIQUE(sampleId, datafileId) ON CONFLICT IGNORE,
            FOREIGN KEY(sampleId) REFERENCES samples(sampleId),  
            FOREIGN KEY(datafileId) REFERENCES datafiles(datafileId) )''')
            self.tables.append('samples_datafiles')
            
            curs.execute('DROP TABLE IF EXISTS experiments ')
            curs.execute(''' CREATE TABLE experiments (
            experimentId INTEGER PRIMARY KEY, 
            name TEXT UNIQUE,
            type TEXT,
            description TEXT,
            config BLOB )''')
            self.tables.append('experiments')
            
            curs.execute('DROP TABLE IF EXISTS filtering_parameters ')
            curs.execute(''' CREATE TABLE filtering_parameters (
            filtering_parameterId INTEGER PRIMARY KEY, 
            params BLOB, 
            filtering_summary TEXT,
            cleaning_summary TEXT )''')
            self.tables.append('filtering_parameters')

            curs.execute('DROP TABLE IF EXISTS CDHIT_parameters ')
            curs.execute(''' CREATE TABLE CDHIT_parameters (
            CDHIT_parameterId INTEGER PRIMARY KEY, 
            params TEXT UNIQUE ON CONFLICT IGNORE )''')
            self.tables.append('CDHIT_parameters')
            
#            curs.execute('DROP TABLE IF EXISTS parameters ')
#            curs.execute(''' CREATE TABLE parameters (
#            parameterId INTEGER PRIMARY KEY, 
#            CDHIT_parameters TEXT,
#            filtering_parameters TEXT, 
#            UNIQUE(CDHIT_parameters, filtering_parameters) ON CONFLICT IGNORE )''')
#            self.tables.append('parameters')
            
#            curs.execute('DROP TABLE IF EXISTS experiments_parameters')
#            curs.execute(''' CREATE TABLE experiments_parameters (
#            linkId INTEGER PRIMARY KEY, 
#            experimentId INTEGER, 
#            parameterId INTEGER, 
#            FOREIGN KEY(experimentId) REFERENCES experiments(experimentId)  
#            FOREIGN KEY(parameterId) REFERENCES parameters(parameterId) )''')
#            self.tables.append('experiments_parameters')
            
            curs.execute('DROP TABLE IF EXISTS clust_results')
            curs.execute(''' CREATE TABLE clust_results (
            clust_resultId INTEGER PRIMARY KEY, 
            experimentId INTEGER,
            datafileId INTEGER, 
            CDHIT_parameterID TEXT,
            cluster_counter BLOB )''')
            self.tables.append('clust_results')
            
    def add_barcodes_datafiles(self, barcodefiles, datafiles, datafile_type=None):
        ''' Add entries to samples, datafiles and mappings table given a list of barcodefiles 
        and a list of datafiles.'''
           
        with self.con as con:
            curs = con.cursor()

            for datafile in datafiles:
                
                # Do not include path tail if its present
                datafile = os.path.split(datafile)[1]
                
                curs.execute('''INSERT OR IGNORE INTO datafiles(filename, type)
                                        VALUES(?,?)''', (datafile, datafile_type))
        
                for barcode in barcodefiles:
                    with open(barcode, 'r') as f: 
                        for line in f:
                            line = line.strip().split()
                            MIDtag = line[0] 
                            description = line[1]
            
                            curs.execute('''INSERT OR IGNORE INTO samples(MIDtag, description)
                                            VALUES(?,?)''', (MIDtag, description))
                    
                            # Find ID of last scanned Barcode                             
                            curs.execute('''SELECT sampleId FROM samples WHERE description=?''', (description,))
                            sample_id = curs.fetchone()['sampleId']
                            
                            # Find ID of last scanned datafile                             
                            curs.execute('''SELECT datafileId FROM datafiles WHERE filename=?''', (datafile,))
                            datafile_id = curs.fetchone()['datafileId']
                    
                            curs.execute('''INSERT INTO samples_datafiles(sampleId, datafileId)
                                            VALUES(?,?)''', (sample_id, datafile_id))
            
    def add_datafile(self, datafile, sample_desc_list, datafile_type=None ):
        ''' Add ONE new file to the datafile table along with appropriate entries into 
        samples_datafiles mappings table. '''
        
        with self.con as con:
            curs = con.cursor()
            
            curs.execute('''INSERT OR IGNORE INTO datafiles(filename, type)
                                        VALUES(?,?)''', (datafile,datafile_type))
            # Find ID of datafile                             
            curs.execute('''SELECT datafileId FROM datafiles WHERE filename=?''', (datafile,))
            datafile_id = curs.fetchone()['datafileId']
        
            for sample_desc in sample_desc_list:
                
                # Find ID of last scanned Barcode                             
                curs.execute('''SELECT sampleId FROM samples WHERE description=?''', (sample_desc,))
                sample_id = curs.fetchone()['sampleId']
                
                # Update Mapping table
                curs.execute('''INSERT INTO samples_datafiles(sampleId, datafileId)
                                            VALUES(?,?)''', (sample_id, datafile_id))
            return datafile_id 

    def add_experiment(self, config, exp_type='', name=None, description=None):
        ''' Add appropriate entries into experiment table. '''

        if name == None:
            name = config.experiment_name
        if description == None:
            description = config.experiment_description

        with self.con as con:
            curs = con.cursor()

            # Update experiments table            
            curs.execute('''INSERT INTO experiments(name, type, description) VALUES (?,?,?)'''.format(),
                            (config.experiment_name, exp_type, config.experiment_description) )
            exp_id = curs.lastrowid
            self.update_binary(config, col='config', target='experimentId', value=exp_id, table='experiments')
            
            return exp_id
    
    def add_results_parameters_datafiles(self, in_filename, out_filename, counter, config, CDHIT_params):
        ''' Add appropriate entries into results table for samplesID, ExperimetID, paramaterID,  
        datafilesID and cluster_counter.
        
        Also Updates related tables: parameters, experiments_parameters, datafiles and samples_datafiles
        '''    
        
        with self.con as con:
            curs = con.cursor()
            
            # Update Parameters table
            curs.execute('''INSERT OR IGNORE INTO CDHIT_parameters(params)
                                VALUES(?)''', (CDHIT_params,))
            # Find parameter ID, (entry may already exist so not necessarily lastID added)                           
            curs.execute('''SELECT CDHIT_parameterId FROM CDHIT_parameters WHERE params=?''', (CDHIT_params,))
            param_id = curs.fetchone()['CDHIT_parameterId']
            
            # Update Datafiles
            # Get samples description list for clustering input file  
            if in_filename.endswith('.fasta'):
                fname = in_filename.split('.')[:-1] + ['.bgzf']
                fname = ''.join(fname)
            recs = self.get_samples4datafile(fname, fields=['description'])
            sample_desc_list = [rec['description'] for rec in recs]
            # Add new entry in datafiles table and samples_datafiles mapping table 
            datafile_id = self.add_datafile(out_filename, sample_desc_list, datafile_type='Clusters')
            
            # Update Results table
            curs.execute('''INSERT INTO clust_results(experimentId, CDHIT_parameterId, datafileId)
                                VALUES(?,?,?)''', (config.exp_id, param_id, datafile_id))
            res_id = curs.lastrowid
            self.update_binary(counter, col='cluster_counter', target='clust_resultId', value=res_id, table='clust_results')
            
    def get_samples4datafile(self, filename, fields=['MIDtag']):
        ''' Return samples that are present for a given filename.
        
        Fields is a list f columns to return for each row '''
        
        querry1 = '''SELECT {0} 
                    FROM datafiles NATURAL JOIN samples_datafiles NATURAL JOIN samples 
                    WHERE filename GLOB ? '''.format(', '.join(fields))
    
#        # Equivilant querries
#        querry2 = '''SELECT {0} FROM datafiles d, samples_datafiles sd, samples s 
#                    WHERE d.datafileId=sd.datafileId 
#                    AND s.sampleId=sd.sampleId 
#                    AND filename GLOB ? '''.format(', '.join(fields))
#
#        querry3 = '''SELECT {0} 
#                    FROM datafiles AS d INNER JOIN samples_datafiles AS sd ON d.datafileId=sd.datafileId
#                    INNER JOIN samples AS s ON sd.sampleId=s.sampleId
#                    WHERE filename GLOB ? '''.format(', '.join(fields))
        
        with self.con as con:
        
            curs = con.cursor()
            curs.execute(querry1, (filename,))
            records = curs.fetchall()

        return records
    
    def get_clust_results4sample(self, fields, col, target):
        ''' Return record from '''
        
        querry1 = '''SELECT {0} 
                    FROM samples NATURAL JOIN samples_datafiles 
                    NATURAL JOIN datafiles NATURAL JOIN clust_results 
                    WHERE {1} GLOB {2}'''.format(fields, col, target)
        
        with self.con as con:
        
            curs = con.cursor()
            curs.execute(querry1, (fields, col, target))
            records = curs.fetchall()

        return records
    
    def get_cluster_counters4sample(self, sample_description="*", clust_parameters="*"): # fields='cluster_counter', col=, target):
        ''' Return list of cluster_counter dictionaries for a a given sample description '''
        
        table_name = 'samples2results'
        
#        querry1 = '''SELECT cluster_counter  
#                    FROM samples NATURAL JOIN samples_datafiles 
#                    NATURAL JOIN datafiles NATURAL JOIN clust_results 
#                    WHERE {1} GLOB {2}'''.format(fields, col, target)
        
        querry2 = '''SELECT {0} FROM {1} WHERE {2} GLOB ? AND {3} GLOB ?'''.format('cluster_counter, CDHIT_parameters', 
                        table_name, 'description', 'CDHIT_parameters')
        
        with self.con as con:
            curs = con.cursor()
        
            curs.execute('''CREATE VIEW IF NOT EXISTS {0} AS SELECT * 
                    FROM samples NATURAL JOIN samples_datafiles 
                    NATURAL JOIN datafiles NATURAL JOIN clust_results 
                    NATURAL JOIN CDHIT_parameters'''.format(table_name))
        
            curs.execute(querry2, (sample_description, clust_parameters))
            
            records = curs.fetchall()
            data = []
            for row in records:
                pickled_data = str(row[0])
                data.append((pkl.loads(pickled_data), row[1]))

        return data
    
    def get_cluster_counters4experiment(self, exp_name):
        ''' Return list of cluster_counter dictionaries for a a given sample description '''
        
        table_name = 'experiment2results'
        
#        querry1 = '''SELECT cluster_counter  
#                    FROM samples NATURAL JOIN samples_datafiles 
#                    NATURAL JOIN datafiles NATURAL JOIN clust_results 
#                    WHERE {1} GLOB {2}'''.format(fields, col, target)
        
        querry2 = '''SELECT {0} FROM {1} WHERE {2} GLOB ?'''.format('cluster_counter', 
                        table_name, 'name')
        
        with self.con as con:
            curs = con.cursor()
        
            curs.execute('''CREATE VIEW IF NOT EXISTS {0} AS SELECT * 
                    FROM experiments NATURAL JOIN clust_results NATURAL JOIN datafiles'''.format(table_name))
        
            curs.execute(querry2, (exp_name,))
            
            records = curs.fetchall()
            data = []
            for row in records:
                pickled_data = str(row[0])
                data.append(pkl.loads(pickled_data))

        return data
        
    def get_cluster_counters(self, fields="cluster_counter", target1='', value1='', target2=None, value2=None): # fields='cluster_counter', col=, target):
        ''' Return list of cluster_counter dictionaries for a a given sample description. 
        
        Also returns a list of datafiles for those cluster counters, and/or parameters if required
        '''
    
        table_name = 'all'
        
        querry = '''SELECT {0} FROM {1} WHERE {2} GLOB ?'''.format(
                    fields, table_name, target1)
        
        if target2 and value2:
            querry += " AND {3} GLOB ?".format(target2)
        
        with self.con as con:
            curs = con.cursor()
        
            curs.execute('''CREATE VIEW IF NOT EXISTS {0} AS SELECT * 
                    FROM samples NATURAL JOIN samples_datafiles 
                    NATURAL JOIN datafiles NATURAL JOIN clust_results 
                    NATURAL JOIN CDHIT_parameters NATURAL JOIN filtering_parameters
                    NATURAL JOIN experiments'''.format(table_name))
        
            if value2:
                curs.execute(querry, (value1, value2))
            else:
                curs.execute(querry, (value1,))
                
            records = curs.fetchall()
            data = []
            other_fields = []
            for row in records:
                pickled_data = str(row[0])
                data.append((pkl.loads(pickled_data), row[1]))
                other_fields.append(row[1:])
    
        return data, other_fields
    
                               
if __name__ == '__main__':
    
        # Create Database for Gazelles and Zebras

        # Config setup
        class Config(object):
            pass
        c = Config()
        prefix = get_path_prefix()
#        c.barcode_inpath = '/Users/chris/Dropbox/work/code/popGen/testData/barcodes'
        c.barcode_inpath = os.path.join(prefix,'gazelles-zebras', 'barcodes')
        
        db = Popgen_db('gz_samples.db', recbyname=True)
        
#        db.overwrite_table('samples', ('MIDtag TEXT, description TEXT, raw_datafile TEXT, counts INTEGER, readsfile TEXT, counter BLOB'))
        
#        db.overwrite_table('samples', ('''MIDtag TEXT, description TEXT, raw_datafile TEXT,
# counts INTEGER, readsfile TEXT, counter BLOB, clusters BLOB'''))

        
        L6_barcode_files = glob.glob(os.path.join(c.barcode_inpath, '*[6].txt')) 
        
        L8_barcode_files = glob.glob(os.path.join(c.barcode_inpath, '*[8].txt')) 
        
        L6_datafiles = ['L6file1', 'L6file2', 'L6file3' ]
        L8_datafiles = ['L8file1', 'L8file2', 'L8file3' ]
        
        db.add_barcodes_datafiles(L6_barcode_files, L6_datafiles)
        db.add_barcodes_datafiles(L8_barcode_files, L8_datafiles)
        
        data = db.select('* FROM samples') 


''' Useful for when I come to load in images'''

#import matplotlib.pyplot as plt
#import StringIO
#from matplotlib import numpy as np
#
#x = np.arange(0,np.pi*3,.1)
#y = np.sin(x)
#
#fig = plt.figure()
#plt.plot(x,y)
#
#imgdata = StringIO.StringIO()
#fig.savefig(imgdata, format='svg')
#imgdata.seek(0)  # rewind the data
#
#svg_dta = imgdata.buf  # this is svg data
#
#file('test.htm', 'w').write(svg_dta)  # test it
