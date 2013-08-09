'''
Created on 27 Jun 2013

@author: musselle
'''
import os, sys, time, glob, gzip 
import cPickle as pkl
from subprocess import PIPE, Popen

import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sqlite3
from utils.general_utilities import set_trace, get_path_prefix

from database.core import SQLdatabase
from utils.fileIO import SeqRecCycler, inputfile_check, outputfile_check
from utils.clusterIO import ClusterObj, parse, sortby

class Reads_db(SQLdatabase):
    ''' Database to hold all information on fastq sequence reads for an experiment
    '''

    def __init__(self, db_file="test.db", recbyname=True, new=False):
       
        SQLdatabase.__init__(self, db_file, recbyname) 
        
    def create_seqs_table(self, table_name='seqs', overwrite=False):

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
            self.seqs_table_name = table_name
            
    def create_meta_table(self, overwrite=False):
        
        with self.con as con:   
        
            con.execute(''' ''')
        
        
            if overwrite:
                con.execute('DROP TABLE IF EXISTS meta')
            con.execute(''' CREATE TABLE meta (
            Id INTEGER PRIMARY KEY NOT NULL,
            lastclusterid INTEGER )''')
    
            self.tables.append('meta')
            
            # Initialise
            con.execute(''' INSERT INTO meta (lastclusterid) Values (0); ''')
    
    #             tableId TEXT,
    #             RADseqType TEXT,
    #             MIDlength INTEGER NOT NULL,
    #             Cutsite1 TEXT, 
    #             Cutsite2 TEXT)''')

    def create_cluster_table(self, table_name='clusters', overwrite=False):
        ''' Make cluster table in database '''
        
        with self.con as con:
    
            curs = con.cursor()

            if overwrite:
                curs.execute('DROP TABLE IF EXISTS {0}'.format(table_name))
            
            curs.execute(''' CREATE TABLE IF NOT EXISTS {0} (
            clusterId INTEGER PRIMARY KEY NOT NULL,
            repseqid INTEGER NOT NULL,
            size INTEGER NOT NULL) '''.format(table_name))
            self.tables.append('{0}'.format(table_name))


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
    
    def write_reads(self, sql_query, file_handle, format='fasta', ignoreup2=0):
        """ Write records returned by the querry to one large fasta or fastq 
        
         file_handle -- A file object or string specifying a filename. 
        
         ignoreup2 -- starting index of sequence that are written, used to miss out 
         cutsite if necessary. 
         
         """
        
        # Output check
        file_handle = outputfile_check(file_handle)
        
        with self.con as con:
            
            tic = time.time()
            print >> sys.stderr, 'Executing sql querry....', 
            record_curs = con.execute(sql_query)
            print >> sys.stderr, ' Done!'
            print >> sys.stderr, 'Records returned in {0}'.format(
                time.strftime('%H:%M:%S', time.gmtime(time.time() - tic))) 
            
            toc = time.time()
            print >> sys.stderr, 'Writing records to {0} format....'.format(format), 
            
            rec_count = 0
            for rec in record_curs:
                seq_rec = SeqRecord(Seq(rec['seq'][ignoreup2:]), id=str(rec['seqid']), description='')
                SeqIO.write(seq_rec, file_handle, format=format)
                rec_count += 1 
            
            print >> sys.stderr, ' Done!'
            print >> sys.stderr, '\n{0} records written successfully to {1}\nin {2}'.format(rec_count,
                file_handle.name, time.strftime('%H:%M:%S', time.gmtime(time.time() - toc))) 
            
            if file_handle.name not in ['<stdout>', '<stderr>']:
                file_handle.close()
                
        return file_handle 
    
    def load_cluster_file(self, cluster_file_handle, cluster_table_name='clusters', 
                          overwrite=False, fmin=2, fmax=None, skipsort=False):
        ''' Load in a clustering file into the database 
        
        By default singletons are not loaded as cutoff = 2
        
        can also set a fmin and fmax threshold if only clusters of a certain 
        size are to be added.
        
        '''
        if type(cluster_file_handle) == str:
            if not cluster_file_handle.endswith('.clstr'):
                cluster_file_handle = cluster_file_handle + '.clstr'
        
        cluster_file_handle = inputfile_check(cluster_file_handle)
        
        if not skipsort:
            # Filter out singletons and sort clusters in accending order
            print >> sys.stderr, 'Sorting cluster file %s ...' % (cluster_file_handle.name)
            sorted_cluster_file = sortby(cluster_file_handle, reverse=True, 
                                         mode='cluster_size', outfile_postfix='-subset', cutoff=fmin)
            cluster_file_handle = inputfile_check(sorted_cluster_file)

        
        print >> sys.stderr, 'Importing cluster file %s  to database...' % (cluster_file_handle.name)

        cluster_gen = parse(cluster_file_handle)
        
        if overwrite:
            with self.con as con:
                con.execute(''' DROP TABLE IF EXISTS {0} '''.format(cluster_table_name))
                con.execute(''' UPDATE seqs SET clusterid = NULL''')
        
        # Create clusters table if necessary
        self.create_cluster_table(cluster_table_name)
        
        if fmax:
            for cluster in cluster_gen:
                if cluster.size <= fmax and cluster.size >= fmin:
                    self.load_single_clusterobj(cluster, cluster_table_name, overwrite)
        else:    
            for cluster in cluster_gen:
                self.load_single_clusterobj(cluster, cluster_table_name, overwrite)
        
    def load_cluster_file2(self, cluster_file_handle, cluster_table_name='clusters', 
                          overwrite=False, fmin=2, fmax=None, skipsort=False):
        ''' Load in a clustering file into the database 
        
        By default singletons are not loaded as cutoff = 2
        
        can also set a fmin and fmax threshold if only clusters of a certain 
        size are to be added.
        
        '''
        if type(cluster_file_handle) == str:
            if not cluster_file_handle.endswith('.clstr'):
                cluster_file_handle = cluster_file_handle + '.clstr'
        
        cluster_file_handle = inputfile_check(cluster_file_handle)
        
        if not skipsort:
            # Filter out singletons and sort clusters in accending order
            print >> sys.stderr, 'Sorting cluster file %s ...' % (cluster_file_handle.name)
            sorted_cluster_file = sortby(cluster_file_handle, reverse=True, 
                                         mode='cluster_size', outfile_postfix='-subset', cutoff=fmin)
            cluster_file_handle = inputfile_check(sorted_cluster_file)

        
        print >> sys.stderr, 'Importing cluster file %s  to database...' % (cluster_file_handle.name)

        cluster_gen = parse(cluster_file_handle)
        
        if overwrite:
            with self.con as con:
                con.execute(''' DROP TABLE IF EXISTS {0} '''.format(cluster_table_name))
                con.execute(''' UPDATE seqs SET clusterid = NULL''')
        
        # Create clusters table if necessary
        self.create_cluster_table(cluster_table_name)
        
        if fmax:
            for cluster in cluster_gen:
                if cluster.size <= fmax and cluster.size >= fmin:
                    self.load_single_clusterobj2(cluster, cluster_table_name, overwrite)
        else:    
            for cluster in cluster_gen:
                self.load_single_clusterobj2(cluster, cluster_table_name, overwrite)
        
        
    def load_single_clusterobj(self, clusterobj, cluster_table_name, overwrite=False):
        ''' Load in a single cluster object to the database. '''
        
        seq_table_name = 'seqs'
        
        with self.con as con:
            
            # Just insert cluster into clusters table and incriment id 
            curs = con.cursor()
            
            curs.execute(''' INSERT INTO {0}
                    (repseqid, size) VALUES (?,?)'''.format(cluster_table_name), 
                    (int(clusterobj.rep_seq_id) , clusterobj.size) )
                
            last_clusterid = curs.lastrowid
                
            # Process members
            
            # Rewrite this with execdb =ute many
            
            
            
            for seqid in clusterobj.members_id:
                con.execute(''' UPDATE {0} SET
                    clusterId = ? WHERE seqId = ?'''.format(seq_table_name), 
                    (last_clusterid, int(seqid)))

#             # Just insert cluster into clusters table and incriment id 
#             curs = con.cursor()
#             curs.execute(''' INSERT INTO {0}
#                     (clusterid, repseqid, size) VALUES (?,?,?)'''.format(cluster_table_name), 
#                     ( clusterobj.id, int(clusterobj.rep_seq_id) , clusterobj.size) )
#                 
#             last_clusterid = curs.lastrowid
#
#             # Process members
#             for seqid in clusterobj.members_id :
#                 con.execute(''' UPDATE {0} SET
#                     clusterId = ? WHERE seqId = ?'''.format(seq_table_name), 
#                     (clusterobj.id, int(seqid)))
                
    def load_single_clusterobj2(self, clusterobj, cluster_table_name, overwrite=False):
        ''' Load in a single cluster object to the database. '''
        
        seq_table_name = 'seqs'
        
        with self.con as con:
            
            # Just insert cluster into clusters table and incriment id 
            curs = con.cursor()
            
            curs.execute(''' INSERT INTO {0}
                    (repseqid, size) VALUES (?,?)'''.format(cluster_table_name), 
                    (int(clusterobj.rep_seq_id) , clusterobj.size) )
                
            last_clusterid = curs.lastrowid
                
            # Process members
            
            # Rewrite this with execute many
            
            member_tuples = [( last_clusterid, int(seqid) ) for seqid in clusterobj.members_id]    
            
            con.executemany(''' UPDATE {0} SET
                    clusterId = ? WHERE seqId = ?'''.format(seq_table_name), member_tuples)

#             # Just insert cluster into clusters table and incriment id 
#             curs = con.cursor()
#             curs.execute(''' INSERT INTO {0}
#                     (clusterid, repseqid, size) VALUES (?,?,?)'''.format(cluster_table_name), 
#                     ( clusterobj.id, int(clusterobj.rep_seq_id) , clusterobj.size) )
#                 
#             last_clusterid = curs.lastrowid
#
#             # Process members
#             for seqid in clusterobj.members_id :
#                 con.execute(''' UPDATE {0} SET
#                     clusterId = ? WHERE seqId = ?'''.format(seq_table_name), 
#                     (clusterobj.id, int(seqid)))
                

if __name__ == '__main__':
    pass