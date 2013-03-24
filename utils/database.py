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
        
        database_already_exists = os.path.exists(dbfile)    
        
        # Stored Vars
        self.con = sqlite3.connect(dbfile)
        self.tables = []
        
        if database_already_exists:
            print 'Database found with matching file name.'
            print 'Connected to database {0}'.format(dbfile)
        else:
            print 'Creating new Database file: {0}'.format(dbfile) 
                
        if recbyname: 
            self.con.row_factory = sqlite3.Row
            print 'Setting Row_factory to named Rows' 

        self.dbname = dbfile
        self.con.close()
        
    def new_table(self, name, headers):
        """ Create a table with specified headers """
        
        self.tables.append(name)
        
        con = sqlite3.connect(self.dbfile) 
        
        with con:        
            cur = con.cursor()
            cur.execute("CREATE TABLE {0} (id INTEGER PRIMARY KEY, {1})".format(name, headers))        
        
    def overwrite_table(self, name, headers):
        """ Create or overwrite a table if it already exists with specified headers """
        
        if name not in self.tables: 
            self.tables.append(name)
        
        cur = self.con.cursor()
        try:
            cur.execute("DROP TABLE IF EXISTS {0}".format(name))
            cur.execute("CREATE TABLE {0} (id INTEGER PRIMARY KEY, {1})".format(name, headers))       
        except sqlite3.Error as e:                        
            print "Error %s:" % e.args[0]
        finally:
            self.con.commit()
            cur.close()    
        
    def select(self,cmd):
        """ select records from the database returned as tuple of tuples. """
        cur = self.con.cursor()  
        # SELECT column_name(s) FROM table_name
        cur.execute("SELECT " + cmd)
        records = cur.fetchall()
        cur.close()
        return records
    
    def display(self, cmd, num = 0):
        """ Print out the records returned from the database querry 
        
        num = max number to display. default 0 is all records returned.
        """    
        cur = self.con.cursor()  
        cur.execute("SELECT " + cmd)       
        if num: 
            records = cur.fetchmany(num)
        else:  
            records = cur.fetchall()      
        cur.close()   
        if self.con.row_factory == sqlite3.Row:
            row_names = records[0].keys()
            print '  '.join(row_names)
            for row in records:
                print row
        else: 
            for row in records:
                print row
        
    def insert(self,cmd):
        """ insert a new record to database and return the new primary key """
        newID = 0
        cur = self.con.cursor()            
        cur.execute("INSERT " + cmd)
        newID = cur.lastrowid
        self.con.commit()
        cur.close()
        return newID
        
    def update(self, cmd,*args):
        """ Set values in a table """        
        cur = self.con.cursor()            
        # UPDATE table_name SET column1=value, column2=value,... WHERE some_column=some_value
        cur.execute("UPDATE " + cmd, *args)
        self.con.commit()
        cur.close()    
    
    def add_binary(self, obj, col, id, table=None):
        """ Store a binary object to the database """
        
        if table is None: 
            table=self.tables[0]       
        cur = self.con.cursor()
        b = sqlite3.Binary(pkl.dumps(obj))             
        # UPDATE table_name SET column1=value, column2=value,... WHERE some_column=some_value
        cur.execute('UPDATE {0} SET {1}=? WHERE id=?'.format(table, col), (b,id))
        self.con.commit()
        cur.close()  
        
    def get_binary(self, col, id, table):
        """ Retrieve binary object in the database """
        
        if table is None: 
            table=self.tables[0]
        cur = self.con.cursor()
        cur.execute('SELECT {0} FROM {1} WHERE id=?'.format(col, table), (id,))        
        pickled_data = str(cur.fetchone()[0])
        data = pkl.loads(pickled_data)
        cur.close()
        return data   
        
    def execute(self,sql, *args):
        """ execute any SQL statement but no return value given """
        cur = self.con.cursor()  
        cur.execute(sql, *args)
        records = cur.fetchall()
        self.con.commit()
        cur.close()
        return records

    def executescript(self,sql, *args):
        """ execute any SQL statement but no return value given """
        cur = self.con.cursor()  
        cur.executescript(sql, *args)
        self.con.commit()
        cur.close()
        
        
#pdata = cPickle.dumps(data, cPickle.HIGHEST_PROTOCOL)
#curr.execute("insert into table (data) values (:data)", sqlite3.Binary(pdata))
#        
#        
#curr.execute("select data from table limit 1")
#for row in curr:
#  data = cPickle.loads(str(row['data']))     
#        
            
               
#    def register_pickler(self, obj):
#        """ Registers an adapter function and converter function that used 
#        pickle.dumps to store the passed object """
#        
#        def adapter_func(obj):
#            """Convert from in-memory to storage representation.
#            """
#            print 'adapter_func(%s)\n' % obj
#            return pickle.dumps(obj)
#        
#        def converter_func(data):
#            """Convert from storage to in-memory representation.
#            """
#            print 'converter_func(%r)\n' % data
#            return pickle.loads(data)
#
#        # Register the functions for manipulating the type.
#        sqlite3.register_adapter(obj, adapter_func)
#        sqlite3.register_converter(str(obj), converter_func)
                
  
class PopGen_DB(Database):  
   
    def __init__(self, db_file="sample.db", table_name=None, column_headers=None, recbyname=True):
       
        table_name = 'samples'
        table_headers = 'MIDtag TEXT, description TEXT, raw_datafile TEXT, count INTEGER, readsfile TEXT' 
       
        Database.__init__(self, db_file, table_name, table_headers, recbyname) 
    
    def add_barcodes(self, barcode_file, raw_datafile):   
        curs = self.con.cursor()
    
        with open(barcode_file, 'r') as f: 
            for line in f:
                line = line.strip().split()
                MIDtag = line[0] 
                description = line[1]

                curs.execute('''INSERT INTO samples (MIDtag, description, raw_datafile)
                                VALUES(?,?,?)''', (MIDtag, description, raw_datafile)) 
        self.con.commit()
        curs.close()   

    def add_barcodes_datafiles(self, barcodefiles, datafiles):
        ''' Creates samples table, datafile table and mappings table given a list of barcodefiles 
        and a list of datafiles.'''
           
        # Not sure if i need this 
        self.conn.execute('pragma foreign_keys=ON')
        
        curs = self.con.cursor()
        
        # Create Tables
        curs.execute(''' CREATE TABLE samples (
        Id INTEGER PRIMARY KEY AUTOINCREMENT, 
        MIDtag TEXT, 
        description TEXT, 
        raw_datafiles INTEGER, 
        FOREIGN KEY(raw_datafiles) REFERENCES datafiles(Id) )''')
        
        curs.execute(''' CREATE TABLE datafiles (
        Id INTEGER PRIMARY KEY AUTOINCREMENT, 
        filename TEXT  )''')
        
        curs.execute(''' CREATE TABLE samples_datafiles (
        Id INTEGER PRIMARY KEY AUTOINCREMENT, 
        sampleId INTEGER, 
        datafileId INTEGER, 
        FOREIGN KEY(sampleId) REFERENCES samples(Id)  
        FOREIGN KEY(datafileId) REFERENCES datafiles(Id) )''')
        
        for datafile in datafiles:

            curs.execute('''INSERT INTO datafiles(filename)
                                    VALUES(?)''', (datafile,))
                
            last_datafile = curs.lastrowid
    
            for barcode in barcode_files:
        
                with open(barcode, 'r') as f: 
                    for line in f:
                        line = line.strip().split()
                        MIDtag = line[0] 
                        description = line[1]
        
                        curs.execute('''INSERT INTO samples(MIDtag, description)
                                        VALUES(?,?)''', (MIDtag, description))
                        last_sample = curs.lastrowid
                
                        curs.execute('''INSERT INTO samples_datafiles(sampleId, datafileId)
                                        VALUES(?,?)''', (last_sample, last_datafile))
                
                self.con.commit()
                curs.close()   
            
    def add_link_table(self):
        
        self.db.execute(''' CREATE TABLE samples_datafiles
        
        ''')                
        
        
        INSERT INTO child VALUES(NULL, 'bobby');
SELECT last_insert_rowid(); -- gives the id of bobby, assume 2 for this example
INSERT INTO dog VALUES(NULL, 'spot');
SELECT last_insert_rowid(); -- gives the id of spot, assume 4 for this example
INSERT INTO child_dog VALUES(2, 4);
        
                             
                             
# Querry would be something like 
c.execute(''.join([
            'SELECT MIDtags ',
            'FROM samples ',
            'WHERE variable.id=var_func.variable_id ', 
            'AND function.id=var_func.function_id ',
            'AND variable.name="wert" ',
            'AND var_func.type="input" ',
            ]))                             
                             
                             
                               
if __name__ == '__main__':
    
        # Create Database for Gazelles and Zebras
        import cPickle as pickle
        
        # Config setup
        class Config(object):
            pass
        c = Config()
        prefix = get_data_prefix()
        c.barcode_inpath = os.path.join(prefix,'gazelles-zebras/barcodes')
        
        db = PopGen_DB('gz_samples.db', recbyname=True)
        
        db.overwrite_table('samples', ('MIDtag TEXT, description TEXT, raw_datafile TEXT, counts INTEGER, readsfile TEXT, counter BLOB'))
        
        L6_barcode_files = glob.glob(os.path.join(c.barcode_inpath, '*[6].txt')) 
        L8_barcode_files = glob.glob(os.path.join(c.barcode_inpath, '*[8].txt')) 
        
        for f in L6_barcode_files:
            db.add_barcodes(f, 'lane6*.bgzf') 
        for f in L8_barcode_files:
            db.add_barcodes(f, 'lane8*.bgzf')
            
        myarray = np.ones([5,5]) * np.random.rand(5,5)
        
        pikled_str_array = pickle.dumps(myarray) 
        
        bin_array = sqlite3.Binary(pikled_str_array)
              
        # UPDATE table_name SET column1=value, column2=value,... WHERE some_column=some_value 
        db.execute('UPDATE samples SET counter=? WHERE id=1', (bin_array,))
            
        data = db.select('* FROM samples') 

        # Test getting and setting binary 
        path = '/space/musselle/data/RAD-seq/gazelles-zebras/stats/L6PrimerTagsCounter.pkl'
        obj = pkl.load(open(path, 'r'))
        db.add_binary(obj, 'counter', id=1)
        got_obj = db.get_binary('counter', id=1, table='samples')
        
        assert obj == got_obj
        

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

        