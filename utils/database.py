'''
Created on 19 Mar 2013

class to interface with a sqlite database

@author: musselle
'''
import os
import sys
import cPickle as pkl
from subprocess import PIPE, Popen
import glob

import gzip
import numpy as np
from Bio import SeqIO

import sqlite3
from utils import get_data_prefix

class Database(object):
    """ Class to handle all python communication with a sqlite database file """
    def __init__(self, dbfile="default.db", recbyname=True):
        """ 
        recbyname - sets returned records by select to be callable by column names 
        """        
        if os.path.exists(dbfile):
            print 'Database found with matching file name.'
            print 'Connecting to database {0}'.format(dbfile)
        else:
            print 'Creating new Database file: {0}'.format(dbfile) 
                
        # Stored Vars
        self.con = sqlite3.connect(dbfile)
        self.dbfile = os.path.abspath(dbfile) 
        self.recbyname = recbyname
        self.tables = []

        if recbyname: 
            self.con.row_factory = sqlite3.Row
            print 'Setting Row_factory to named Rows' 

    def connect(self):
        ''' Open connection to the database '''
        self.con = sqlite3.connect(self.dbfile)
        if self.recbyname: 
            self.con.row_factory = sqlite3.Row
            print 'Setting Row_factory to named Rows' 
    
    def close(self):
        ''' Close connection to database'''
        self.con.close()
        
    def new_table(self, name, headers):
        """ Create a table with specified headers """       
        self.tables.append(name)
        
        with self.con as con:        
            cur = con.cursor()
            cur.execute("CREATE TABLE {0} (id INTEGER PRIMARY KEY, {1})".format(name, headers))        
        
    def overwrite_table(self, name, headers):
        """ Create or overwrite a table if it already exists with specified headers """       
        if name not in self.tables: 
            self.tables.append(name)
        
        with self.con as con:        
            cur = con.cursor()
            cur.execute("DROP TABLE IF EXISTS {0}".format(name))
            cur.execute("CREATE TABLE {0} (id INTEGER PRIMARY KEY, {1})".format(name, headers))     
        
    def select(self,cmd, *args):
        """ Select records from the database returned as tuple of Row Objects. """
        
        with self.con as con:        
            cur = con.cursor()
            # SELECT column_name(s) FROM table_name
            cur.execute("SELECT " + cmd, *args)
            records = cur.fetchall()
        
        return records

    def display(self, cmd, num=0, *args):
        """ Print out the records returned from the database querry 
        
        num = max number to display. default 0 is all records returned.
        """    
        with self.con as con:     
            cur = con.cursor()
            # SELECT column_name(s) FROM table_name
            cur.execute("SELECT " + cmd, *args)
            if num: 
                records = cur.fetchmany(num)
            else:  
                records = cur.fetchall()      
        
        if self.con.row_factory == sqlite3.Row:
            row_names = records[0].keys()
            print '  '.join(row_names)
            for row in records:
                print row
        else: 
            for row in records:
                print row
        
    def select_as_tuples(self, cmd, *args):
        ''' Select records from the database returned as tuple of tuples. '''
        tuples = []
        with self.con as con:     
            cur = con.cursor()
            # SELECT column_name(s) FROM table_name
            cur.execute("SELECT " + cmd, *args)
            records = cur.fetchall()
            for row in records:
                tuples.append(row)
        return tuples
    
    def insert(self,cmd):
        """ insert a new record to database and return the new primary key """
        newID = 0
        with self.con as con:        
            cur = con.cursor()
            cur.execute("INSERT " + cmd)
            newID = cur.lastrowid
        
        return newID
        
    def update(self, cmd, *args):
        """ Set values in a table """        

        with self.con as con:        
            cur = con.cursor()
            # UPDATE table_name SET column1=value, column2=value,... WHERE some_column=some_value
            cur.execute("UPDATE " + cmd, *args)
        
    def add_binary(self, obj, col, id, table=None):
        """ Store a binary object to the database """

        if table is None: 
            table=self.tables[0]       

        with self.con as con:        
            cur = con.cursor()
            b = sqlite3.Binary(pkl.dumps(obj))             
            # UPDATE table_name SET column1=value, column2=value,... WHERE some_column=some_value
            cur.execute('UPDATE {0} SET {1}=? WHERE id=?'.format(table, col), (b,id))
        
    def get_binary(self, col, id, table):
        """ Retrieve binary object in the database """
        
        if table is None: 
            table=self.tables[0]

        with self.con as con:        
            cur = con.cursor()
            cur.execute('SELECT {0} FROM {1} WHERE id=?'.format(col, table), (id,))        
            pickled_data = str(cur.fetchone()[0])
            data = pkl.loads(pickled_data)
        
        return data   
        
    def execute(self,sql, *args):
        """ execute any SQL statement but no return value given """
        
        with self.con as con:        
            cur = con.cursor()
            cur.execute(sql, *args)
            records = cur.fetchall()
        return records

    def executescript(self,sql, *args):
        """ execute any SQL statement but no return value given """
        
        with self.con as con:        
            cur = con.cursor()
            cur.executescript(sql, *args)
        
class PopGen_DB(Database):  
   
    def __init__(self, db_file="sample.db", recbyname=True):
       
        Database.__init__(self, db_file, recbyname) 
    
        if not os.path.exists(db_file):
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
            
            curs.execute('DROP TABLE IF EXISTS datafiles ')
            curs.execute(''' CREATE TABLE datafiles (
            datafileId INTEGER PRIMARY KEY, 
            filename TEXT UNIQUE,
            type TEXT )''')
            
            curs.execute('DROP TABLE IF EXISTS samples_datafiles')
            curs.execute(''' CREATE TABLE samples_datafiles (
            linkId INTEGER PRIMARY KEY, 
            sampleId INTEGER, 
            datafileId INTEGER, 
            FOREIGN KEY(sampleId) REFERENCES samples(sampleId)  
            FOREIGN KEY(datafileId) REFERENCES datafiles(datafileId) )''')
            
            curs.execute('DROP TABLE IF EXISTS experiments ')
            curs.execute(''' CREATE TABLE experiments (
            experimentId INTEGER PRIMARY KEY, 
            name TEXT UNIQUE,
            type TEXT,
            description TEXT, 
            config BLOB )''')
            
            curs.execute('DROP TABLE IF EXISTS parameters ')
            curs.execute(''' CREATE TABLE parameters (
            parameterId INTEGER PRIMARY KEY, 
            CD-HIT_parameters TEXT,
            filtering_parameters INTEGER, 
            experimentId INTEGER, 
            FOREIGN KEY(experimentId) REFERENCES experiments(experimentId) )''')
            
            curs.execute('DROP TABLE IF EXISTS results ')
            curs.execute(''' CREATE TABLE results (
            resultId INTEGER PRIMARY KEY, 
            sampleId INTEGER,
            experimentId INTEGER )''')
            
            curs.execute('DROP TABLE IF EXISTS results_datafiles ')
            curs.execute(''' CREATE TABLE results_datafiles (
            resultId INTEGER PRIMARY KEY, 
            datafileId INTEGER )''')
            
    def add_barcodes_datafiles(self, barcodefiles, datafiles):
        ''' Add entries to samples, datafiles and mappings table given a list of barcodefiles 
        and a list of datafiles.'''
           
        with self.con as con:
            curs = con.cursor()

            for datafile in datafiles:
                curs.execute('''INSERT OR IGNORE INTO datafiles(filename)
                                        VALUES(?)''', (datafile,))
        
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
            
    def add_datafile(self, datafile, sample_desc_list):
        ''' Add ONE new file to the datafile table along with appropriate entries into 
        samples_datafiles mappings table. '''
        
        with self.con as con:
            curs = con.cursor()
            
            curs.execute('''INSERT OR IGNORE INTO datafiles(filename)
                                        VALUES(?)''', (datafile,))
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

    def add_experiment_and_samples(self, c, barcodefiles):
        ''' Add appropriate entries into mappings table for samples_experiments. '''

        with self.con as con:
            curs = con.cursor()

            # Update experiments table            
            curs.execute('''INSERT OR REPLACE INTO experiments(name, type, description) VALUES (?,?,?)''',
                            (c.experiment_name, 'clustering', c.experiment_description) )
            expid = curs.lastid
            self.add_binary(c, 'config', id=expid, table='experiments')
            
            for barcode in barcodefiles:
                with open(barcode, 'r') as f: 
                    for line in f:
                        line = line.strip().split()
                        MIDtag = line[0] 
                        description = line[1]

                        # Find ID of last scanned Barcode                             
                        curs.execute('''SELECT sampleId FROM samples WHERE description=?''', (description,))
                        sample_id = curs.fetchone()['sampleId']
                        
                        curs.execute('''INSERT INTO results(sampleId, datafileId)
                                        VALUES(?,?)''', (sample_id, expid))



            
#    def add_experiment(self, experiment_details_dict, type, ):
#        ''' Add an experiment to the '''
#            

    def update_repalce_table(self, tablename, id = , column_values_dict):
        ''' Updates values in tablename with values in the passed dictionary.
        
        If id is not given, a new row is inserted using INSERT OR REPLACE '''
        
            
            
            
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
                               
if __name__ == '__main__':
    
        # Create Database for Gazelles and Zebras
        import cPickle as pickle
        
        # Config setup
        class Config(object):
            pass
        c = Config()
        prefix = get_data_prefix()
#        c.barcode_inpath = '/Users/chris/Dropbox/work/code/popGen/testData/barcodes'
        c.barcode_inpath = os.path.join(prefix,'gazelles-zebras', 'barcodes')
        
        db = PopGen_DB('gz_samples.db', recbyname=True)
        
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

        