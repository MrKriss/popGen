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
from utils import get_path_prefix
from general_utilities import set_trace

from database.core import SQLdatabase
from utils.fileIO import SeqRecCycler


class Reads_db(SQLdatabase):
    ''' Database to hold all information on fastq sequence reads for an experiment
    '''

    def __init__(self, db_file="test.db", recbyname=True, new=False):
       
        SQLdatabase.__init__(self, db_file, recbyname) 
       
    def create_tables(self, table_name='seqs', overwrite=False):

        with self.con as con:
    
            curs = con.cursor()        
    
            # TODO: Write function to infer description format and use this to creste specific tables.
            # Currently only does Casava 1.8 format

            # Create Tables
#             curs.execute(''' CREATE TABLE seqs (
#             seqId INTEGER PRIMARY KEY NOT NULL,
#             
#             MIDseq TEXT NOT NULL,
#             MIDphred TEXT NOT NULL,
#             seq TEXT NOT NULL,
#             phred TEXT NOT NULL, 
#             length INTEGER NOT NULL,
#             individualId TEXT,
#             medianPhred INTEGER, 
#             clusterId INTEGER,
#             
#             instrument TEXT, 
#             runId INTEGER,
#             flowcellId TEXT,
#             laneId INTEGER,
#             tileNo INTEGER,
#             xcoord INTEGER,
#             ycoord INTEGER,
#             pairedEnd INTEGER,
#             illuminaFilter TEXT,
#             controlBits INTEGER,
#             indexSeq TEXT) ''')
#             self.tables.append('seqs')
            
            if overwrite:
                curs.execute('DROP TABLE IF EXISTS {0}'.format(table_name))
            curs.execute(''' CREATE TABLE IF NOT EXISTS {0} (
            seqId INTEGER PRIMARY KEY NOT NULL,
            
            MIDseq TEXT NOT NULL,
            MIDphred TEXT NOT NULL,
            seq TEXT NOT NULL,
            phred TEXT NOT NULL, 
            length INTEGER NOT NULL,
            individualId TEXT,
            meanPhred INTEGER, 
            clusterId INTEGER,
            
            description TEXT,
            
            pairedEnd INTEGER,
            illuminaFilter TEXT,
            controlBits INTEGER,
            indexSeq TEXT) '''.format(table_name))
            self.tables.append('{0}'.format(table_name))
            
#             if overwrite:
#                 curs.execute('DROP TABLE IF EXISTS meta')
#             curs.execute(''' CREATE TABLE meta (
#             Id INTEGER PRIMARY KEY NOT NULL,
#             tableId TEXT,
#             RADseqType TEXT,
#             MIDlength INTEGER NOT NULL,
#             Cutsite1 TEXT, 
#             Cutsite2 TEXT)''')
#             self.tables.append('meta')
#             

    def load_seqs(self, data_files=None, barcode_files=None, table_name='seqs'):
        ''' Load in all sequences in the specified files to the database 
        
        Barcodes for the MIDtag samples are added to the database if given.
         
        The list of Barcodes in the files given must be unique. 
        Therefore if there are duplicates, the data and barcodes must be added in sets of
        unique barcodes for the data files added. 
         
        If 'data_files' or 'barcode_files' is a str, it is interpreted as a glob to the 
        data_path / barcode_path respecively.
        '''
        
        start_dir = os.getcwd()
        
        if barcode_files:
            # input checks
            if type(barcode_files) is str:
                
                # Check for path head included in data_files
                if not os.path.split(barcode_files)[0]:
                    # Append the abs path to current directory 
                    current_dir = os.getcwd()
                    barcode_files = os.path.join(current_dir, barcode_files)
                
                # Glob the file pattern 
                barcode_files = glob.glob(barcode_files)
                assert barcode_files, "No barcode files found at destination"
                barcode_files.sort()

            elif type(barcode_files) is list or type(barcode_files) is tuple:
                if not os.path.split(barcode_files[0])[0]:
                    current_dir = os.getcwd()
                    barcode_files = [os.path.join(current_dir, x) for x in barcode_files]
            else:
                raise Exception('Invalid entry for barcode_files.')

            # Store barcode dictionary            
            MIDs = []
            individuals = []
            for barcode_file in barcode_files:
                with open(barcode_file, 'rb') as f: 
                    for line in f:
                        parts = line.strip().split()
                        MIDs.append(parts[0])
                        individuals.append(parts[1])
                        
                    diff = len(MIDs) - len(set(MIDs))
                    if diff > 0:
                        raise Exception('MID tags in barcode files are not unique.\n{0} duplicates found'.format(diff))
                    else:
                        barcode_dict = dict(zip(MIDs, individuals))
                        MIDlength = len(MIDs[0])
               
        # Setup Files generator
        print type(data_files) is file
        if type(data_files) is file:
            recgen = SeqIO.parse(data_files, 'fastq')
        else: 
            SeqRecGen = SeqRecCycler(data_files=data_files)
            recgen = SeqRecGen.recgen
        
        with self.con as con:
            curs = con.cursor()
            for rec in recgen:
                
                # TODO: Write function to infer description format and use this to create specific tables.
                # Currently only does Casava 1.8 format

                # Store sequence MIDtag and phred info 
                fullseq = rec.seq.tostring()
                MIDseq = fullseq[:MIDlength]
                if barcode_files:
                    individualId = barcode_dict[MIDseq]

                # Sequence stored does not include MID tag, this is stored separately.
                seq = fullseq[MIDlength:]
                length = len(seq)
                
                fullphred = [chr(i+33) for i in rec.letter_annotations['phred_quality']]
                MIDphred = fullphred[:MIDlength]
                MIDphred = ''.join(MIDphred)
                phred = fullphred[MIDlength:]
                phred = ''.join(phred)
                meanPhred = round(sum(rec.letter_annotations['phred_quality']) / float(len(fullseq)))
                
                # Store data in Sequence ID string
                data = rec.description.split()
                data = [i.split(':') for i in data]
                # Format is 
                # [['@HWI-ST0747', '233', 'C0RH3ACXX', '6', '2116', '17762', '96407'], ['1', 'N', '0', '']]
                
                description = ':'.join(data[0])
                
#                 instrument = data[0][0]
#                 runId = data[0][1]
#                 flowcellId = data[0][2]
#                 laneId = data[0][3]
#                 tileNo = data[0][4]
#                 xcoord = data[0][5]
#                 ycoord = data[0][6]
                pairedEnd = data[1][0]
                illuminaFilter = data[1][1]
                controlBits = data[1][2]
                indexSeq = data[1][3]
                
                # Insert record into table
#                 curs.execute('''INSERT INTO seqs
#                  (seq, phred, MIDseq, MIDphred, individualId, meanPhred, length, instrument, runId, flowcellId, laneId, tileNo, 
#                  xcoord, ycoord, pairedEnd, illuminaFilter, controlBits, indexSeq) 
#                  VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);''', 
#                  (seq, phred, MIDseq, MIDphred, individualId, medianPhred, length, instrument, runId, flowcellId, laneId, tileNo, 
#                  xcoord, ycoord, pairedEnd, illuminaFilter, controlBits, indexSeq));
                
                curs.execute('''INSERT INTO {0}
                 (seq, phred, MIDseq, MIDphred, individualId, meanPhred, length, description,
                  pairedEnd, illuminaFilter, controlBits, indexSeq) 
                 VALUES (?,?,?,?,?,?,?,?,?,?,?,?);'''.format(table_name), 
                 (seq, phred, MIDseq, MIDphred, individualId, meanPhred, length, description,
                  pairedEnd, illuminaFilter, controlBits, indexSeq));
                

if __name__ == '__main__':
    pass