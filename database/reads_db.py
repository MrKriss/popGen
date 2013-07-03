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
from utils.fileIO import SeqRecCycler


class Reads_db(SQLdatabase):
    ''' Database to hold all information on fastq sequence reads for an experiment
    
    '''

    def __init__(self, db_file="test.db", recbyname=True, new=False):
       
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
            seqId INTEGER PRIMARY KEY NOT NULL,
            
            seq TEXT NOT NULL,
            phred TEXT NOT NULL, 
            length INTEGER NOT NULL,
            individualId TEXT,
            meanPhred NUMBER, 
            clusterId INTEGER,
            
            instrument TEXT, 
            runId INTEGER,
            flowcellId TEXT,
            laneId INTEGER,
            tileNo INTEGER,
            xcoord INTEGER,
            ycoord INTEGER,
            pairedEnd INTEGER,
            illuminaFilter TEXT,
            controlBits INTEGER,
            indexSeq TEXT) ''')
            self.tables.append('seqs')

    def load_seqs(self, data_files=None, data_inpath=None, barcodes=None, barcode_inpath):
        ''' Load in all sequences in the specified files '''
        
        # Setup Files generator
        SeqRecGen = SeqRecCycler(data_files=data_files, data_inpath=data_inpath)
        
        with self.con as con:
            curs = con.cursor()
        
            for rec in SeqRecGen.recgen:
                
                # TODO: Write function to infer description format and use this to creste specific tables.
                # Currently only does Casava 1.8 format
                
                data = rec.description.split()
                data = [i.split(':') for i in data]
                # Format is 
                # [['@HWI-ST0747', '233', 'C0RH3ACXX', '6', '2116', '17762', '96407'], ['1', 'N', '0', '']]
                
                seq = rec.seq.tostring()
                length = len(seq)
                meanPhred = sum(rec.letter_annotations['phred_quality']) / float(length)
                phred = [chr(i+33) for i in rec.letter_annotations['phred_quality']]
                phred = ''.join(phred)
                
                # TODO: Add lookup to barcode individuals and store here
                
                
                instrument = data[0][0]
                runId = data[0][1]
                flowcellId = data[0][2]
                laneId = data[0][3]
                tileNo = data[0][4]
                xcoord = data[0][5]
                ycoord = data[0][6]
                pairedEnd = data[1][0]
                illuminaFilter = data[1][1]
                controlBits = data[1][2]
                indexSeq = data[1][3]
                
                print seq
                print phred
                
                # Insert record into table
                curs.execute('''INSERT INTO seqs
                 (seq, phred, meanPhred, length, instrument, runId, flowcellId, laneId, tileNo, 
                 xcoord, ycoord, pairedEnd, illuminaFilter, controlBits, indexSeq) 
                 VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);''', 
                 (seq, phred, meanPhred, length, instrument, runId, flowcellId, laneId, tileNo, 
                 xcoord, ycoord, pairedEnd, illuminaFilter, controlBits, indexSeq));
                

if __name__ == '__main__':
    pass