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
    def __init__(self, db_file="default.db", table_name=None, column_headers=None):
       
        
        database_already_exists = os.path.exists(db_file)    
        
        self.db = sqlite3.connect(db_file)
        
        if database_already_exists:
            print 'Database found with matching file name.'
            print 'Connected to database {0}'.format(db_file)
        else:
            print 'Creating new Database file: {0}'.format(db_file) 
            if column_headers:
                self.create_table(table_name, column_headers)
        
        
    def new_table(self, name, headers):
        """ Create a table with specified headers """
        
        cur = self.db.cursor()
        try:
            cur.execute("CREATE TABLE {0} (id INTEGER PRIMARY KEY, {1})".format(name, headers))        
        except sqlite3.Error as e:                        
            print "Error %s:" % e.args[0]
        finally:
            self.db.commit()
            cur.close()    
        
    def overwrite_table(self, name, headers):
        """ Create or overwrite a table if it already exists with specified headers """
        
        cur = self.db.cursor()
        try:
            cur.execute("DROP TABLE IF EXISTS {0}".format(name))
            cur.execute("CREATE TABLE {0} (id INTEGER PRIMARY KEY, {1})".format(name, headers))       
        except sqlite3.Error as e:                        
            print "Error %s:" % e.args[0]
        finally:
            self.db.commit()
            cur.close()    
        
    def select(self,cmd):
        """ select records from the database """
        cur = self.db.cursor()  
        cur.execute("SELECT " + cmd)
        records = cur.fetchall()
        cur.close()
        return records   
        
    def insert(self,cmd):
        """ insert a new record to database and return the new primary key """
        newID = 0
        cur = self.db.cursor()            
        cur.execute("INSERT " + cmd)
        newID = cur.lastrowid
        self.db.commit()
        cur.close()
        return newID
        
    def execute(self,sql):
        """ execute any SQL statement but no return value given """
        cur = self.db.cursor()  
        cur.execute(sql)
        records = cur.fetchall()
        self.db.commit()
        cur.close()
        return records

    def executescript(self,sql):
        """ execute any SQL statement but no return value given """
        cur = self.db.cursor()  
        cur.executescript(sql)
        self.db.commit()
        cur.close()
  
class Sample_DB(Database):  
   
    def __init__(self, db_file="sample.db", table_name=None, column_headers=None):
       
        table_name = 'samples'
        table_headers = 'MIDtag TEXT, description TEXT, raw_datafile TEXT, count INTEGER, readsfile TEXT' 
       
        Database.__init__(self, db_file, table_name, table_headers) 
    
    def add_barcodes(self, barcode_file, raw_datafile):
    
        curs = self.db.cursor()
    
        with open(barcode_file, 'r') as f: 
            for line in f:
                line = line.strip().split()
                MIDtag = line[0] 
                description = line[1]

                curs.execute('''INSERT INTO samples (MIDtag, description, raw_datafile)
                                VALUES(?,?,?)''', (MIDtag, description, raw_datafile)) 
        self.db.commit()
        curs.close()
                               
if __name__ == '__main__':
    
        # Create Database for Gazelles and Zebras
        
        # Config setup
        class Config(object):
            pass
        c = Config()
        prefix = get_data_prefix()
        c.barcode_inpath = os.path.join(prefix,'gazelles-zebras/barcodes')
        
        db = Sample_DB('gz_samples.db')
        
        db.overwrite_table('samples', ('MIDtag TEXT, description TEXT, raw_datafile TEXT, count INTEGER, readsfile TEXT'))
        
        L6_barcode_files = glob.glob(os.path.join(c.barcode_inpath, '*[6].txt')) 
        L8_barcode_files = glob.glob(os.path.join(c.barcode_inpath, '*[8].txt')) 
        
        for f in L6_barcode_files:
            db.add_barcodes(f, 'lane6*.bgzf') 
        for f in L8_barcode_files:
            db.add_barcodes(f, 'lane8*.bgzf')
            
        data = db.select('* FROM samples') 
    
    

        