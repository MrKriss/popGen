'''
Created on 27 Jun 2013

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
from general_utilities import set_trace

from database.core import SQLdatabase



class Reads_db(SQLdatabase):
    ''' Database to hold all information on fastq sequence reads for an experiment
    
    '''

    def __init__(self, db_file="sample.db", recbyname=True, new=False):
       
        create_tables = 0
        if not os.path.exists(db_file) or new:
            create_tables = 1

        SQLdatabase.__init__(self, db_file, recbyname) 
    
        if create_tables:
            self.create_tables()
       
    def create_tables(self, overwrite=True):

        with self.con as con:
    
            curs = con.cursor()        
    
            # TODO: Write function to infer description format and use this to creste specific tables.
            # Currently only does Casava 1.8 format

            # Create Tables
            if overwrite:
                curs.execute('DROP TABLE IF EXISTS seqs')
            curs.execute(''' CREATE TABLE seqs (
            seqId INTEGER PRIMARY KEY,
            
            seq TEXT,
            phred TEXT, 
            length INTEGER,
            
            instrument TEXT, 
            runId INTEGER,
            flowcellId TEXT,
            laneId INTEGER,
            tileNo INTEGER,
            xcoord INTEGER,
            ycoord INTEGER,
            paired_end INTEGER,
            illuminaFilter TEXT,
            controlBits INTEGER,
            indexSeq TEXT) ''')
            self.tables.append('seqs')




if __name__ == '__main__':
    pass