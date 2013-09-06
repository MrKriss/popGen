'''
Created on 27 Jun 2013

@author: musselle
'''
import os, sys, time, glob, gzip, csv
import cPickle as pkl
from subprocess import PIPE, Popen

from cStringIO import StringIO

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
    
    4 TABLES:
        SEQS - all info on the reads themselves. Either derived from sequence ID header.
                or calculated as the file is processed.
                
            seqId -     INTEGER PRIMARY KEY
            sampleId -  Maps to a sampled individual in the SAMPLES table (INTEGER), 
            MIDphred  - Phred score for the MIDtag portion of the read (TEXT NOT NULL),
            seq  -      Read portion of the sequence, with the MIDtag already removed  (TEXT NOT NULL),
            phred -     The corresponding Phred score for the sequence read (TEXT NOT NULL), 
            length -    Length of the read (INTEGER NOT NULL),
            meanPhred - Mean phred score for the sequence read, rounded to nearest int (INTEGER), 
            description The text in the seqid header of the read, minus the information extracted bellow TEXT,
            pairedEnd - Whether the sequence is part of a paired end experiment and the numebr (1 or 2) INTEGER,
            illuminaFilter - Whether the read was flagged up by the illumina filter as being too low quality (Y/N) TEXT,
            controlBits - INTEGER,
            indexSeq  - Index sequence used for multiplexing. TEXT
            
        SAMPLES - info on individuals sampled
        
            
        
        
        CLUSTERS - info on clusters 
        
        
        
        MEMBERS - Link table for cluster-seqid membership
        
        
    
    '''

    def __init__(self, db_file="test.db", recbyname=True, new=False):
       
        SQLdatabase.__init__(self, db_file, recbyname) 
        
    def create_seqs_table(self, table_name='seqs', overwrite=False):

        with self.con as con:
    
            curs = con.cursor()        
    
            # TODO: Write function to infer description format and use this to creste specific tables.
            # Currently only does Casava 1.8 format

            
            if overwrite:
                curs.execute('DROP TABLE IF EXISTS {0}'.format(table_name))
            curs.execute(''' CREATE TABLE IF NOT EXISTS {0} (
            seqId INTEGER PRIMARY KEY NOT NULL,
            sampleId INTEGER,
            
            MIDphred TEXT NOT NULL,
            seq TEXT NOT NULL,
            phred TEXT NOT NULL, 
            length INTEGER NOT NULL,
            
            meanPhred INTEGER, 
            description TEXT,
            
            pairedEnd INTEGER,
            illuminaFilter TEXT,
            controlBits INTEGER,
            indexSeq TEXT) '''.format(table_name))
            self.tables.append('{0}'.format(table_name))
            self.seqs_table_name = table_name
            
#     def create_meta_table(self, overwrite=False):
#         
#         with self.con as con:   
#         
#             if overwrite:
#                 con.execute('DROP TABLE IF EXISTS meta')
#             con.execute(''' CREATE TABLE meta (
#             Id INTEGER PRIMARY KEY NOT NULL,
#             exp_name TEXT )''')
#             self.tables.append('meta')
            
    def create_cluster_table(self, table_name=None, overwrite=False):
        ''' Create a cluster table in database with the specified name. 
        
        overwrite - if true will drop any table that already exists with the specified name.
        '''
        
        if table_name is None:
            table_name = 'clusters' 
        
        with self.con as con:
    
            curs = con.cursor()

            if overwrite:
                curs.execute('DROP TABLE IF EXISTS {0}'.format(table_name))
            
            curs.execute(''' CREATE TABLE IF NOT EXISTS {0} (
            clusterId INTEGER PRIMARY KEY NOT NULL,
            selfsimilarity REAL, 
            repseqid INTEGER NOT NULL,
            size INTEGER NOT NULL) '''.format(table_name))
            self.tables.append(table_name)
            
    def create_samples_table(self, overwrite=False):
        
        with self.con as con:   
        
            if overwrite:
                con.execute('DROP TABLE IF EXISTS samples')
            con.execute(''' CREATE TABLE samples (
            sampleId INTEGER PRIMARY KEY NOT NULL,
            MIDtag TEXT,
            description TEXT UNIQUE,
            type TEXT,
            read_count INTEGER )''')
    
            self.tables.append('samples')
            
    def create_members_table(self, table_name=None, overwrite=False):
        ''' Create a members table in database with the specified name. 
        
        overwrite - if true will drop any table that already exists with the specified name.
        '''
        if table_name is None:
            table_name = 'members' 
        
        with self.con as con:   
        
            if overwrite:
                con.execute('DROP TABLE IF EXISTS {0}'.format(table_name))
            con.execute(''' CREATE TABLE {0} (
            linkId INTEGER PRIMARY KEY NOT NULL,
            clusterid INTEGER,
            seqid INTEGER, 
            UNIQUE (clusterid, seqid))'''.format(table_name))
    
            self.tables.append(table_name)
            

    def load_seqs(self, data_files=None, barcode_files=None, table_name='seqs', buffer_max=100000):
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
            descriptions = []
            for barcode_file in barcode_files:
                with open(barcode_file, 'rb') as f: 
                    for line in f:
                        parts = line.strip().split()
                        MIDs.append(parts[0])
                        descriptions.append(parts[1])
                        
                    diff = len(MIDs) - len(set(MIDs))
                    if diff > 0:
                        raise Exception('MID tags in barcode files are not unique.\n{0} duplicates found'.format(diff))
                    else:
                        # Store barcodes in database 
                        ids = []
                        with self.con as con:
                            curs = con.cursor()
                            for i in range(len(MIDs)):
                                curs.execute('''INSERT OR IGNORE INTO samples(MIDtag, description)
                                            VALUES(?,?)''', (MIDs[i], descriptions[i]))
                                # Find ID of last scanned Barcode                             
                                curs.execute('''SELECT sampleId FROM samples WHERE description=?''', (descriptions[i],))
                                ids.append(curs.fetchone()['sampleId'])
                                
                        barcode_dict = dict(zip(MIDs, ids))
                        MIDlength = len(MIDs[0])
                                          
        # Setup Files generator
        if type(data_files) is file:
            recgen = SeqIO.parse(data_files, 'fastq')
        else: 
            SeqRecGen = SeqRecCycler(data_files=data_files)
            recgen = SeqRecGen.recgen
        
        # Setup buffer
        data_buffer = []
        data_buffer_count = 0
        index_name = 'sampleidIndex'
        
        with self.con as con:
            curs = con.cursor()
            
            # Drop index if it exists
            con.execute(''' DROP INDEX IF EXISTS {0} '''.format(index_name))
            
            for rec in recgen:
                
                # TODO: Write function to infer description format and use this to create specific tables.
                # Currently only does Casava 1.8 format

                # Store sequence MIDtag and phred info 
                fullseq = rec.seq.tostring()
                MIDseq = fullseq[:MIDlength]
                if barcode_files:
                    sampleId = barcode_dict[MIDseq]

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
                
                pairedEnd = data[1][0]
                illuminaFilter = data[1][1]
                controlBits = data[1][2]
                indexSeq = data[1][3]
                
                # Write data to memory file 
                data_buffer.append((seq, phred, MIDphred, str(sampleId), str(meanPhred), str(length), description,
                                                pairedEnd, illuminaFilter, controlBits, indexSeq))
                data_buffer_count += 1
                
                if data_buffer_count > buffer_max:
                    # Insert data and reset buffer

                    con.executemany('''INSERT INTO {0}
                     (seq, phred, MIDphred, sampleId, meanPhred, length, 
                     description,pairedEnd, illuminaFilter, controlBits, indexSeq) 
                     VALUES (?,?,?,?,?,?,?,?,?,?,?);'''.format(table_name), data_buffer)
                    
                    data_buffer = []
                    data_buffer_count = 0
            
            # End of generator. Flush remaining data buffer
                    
            con.executemany('''INSERT INTO {0}
             (seq, phred, MIDphred, sampleId, meanPhred, length, description,
              pairedEnd, illuminaFilter, controlBits, indexSeq) 
             VALUES (?,?,?,?,?,?,?,?,?,?,?);'''.format(table_name), data_buffer)
            
#                 curs.execute('''INSERT INTO {0}
#                  (seq, phred, MIDphred, sampleId, meanPhred, length, description,
#                   pairedEnd, illuminaFilter, controlBits, indexSeq) 
#                  VALUES (?,?,?,?,?,?,?,?,?,?,?);'''.format(table_name), 
#                  (seq, phred, MIDphred, sampleId, meanPhred, length, description,
#                   pairedEnd, illuminaFilter, controlBits, indexSeq));
            
            # Rebuild index on sampleID
            con.execute('''CREATE INDEX {indexname} ON {tablename}(sampleId)'''.format(
                    indexname=index_name, tablename=table_name))
    
    def update_type(self, pattern, short_description):
        " Assiciate a symbol to a particular description pattern "
        
        with self.con as con:
            con.execute(''' UPDATE samples SET type = ? WHERE description GLOB ? ''', (short_description,pattern))
            
    
    def get_cluster_by_size(self, size_min, size_max, items=['seqId', 'seq', 'phred', 'sampleId'], table_prefix=None):
        ''' Return all clusterobjs within a certain size range '''
        
        if table_prefix is None:
            clusters_table = 'clusters'
        else:
            clusters_table = table_prefix+ '_clusters'
        
        sql_query = ''' SELECT * FROM {clusters} WHERE size BETWEEN ? AND ?'''.format(clusters=clusters_table)
        
        with self.con as con:
            
            c = con.execute(sql_query, (size_min, size_max))
            clusterobjs = []
            
            # get list of clusterId within the size range
            for row in c:
                clusterobjs.append(self.get_cluster_by_id(row['clusterid'], items=items, table_prefix=table_prefix))
                
        return clusterobjs
    
    def get_cluster_by_id(self, cluster_id, items=['seqId', 'seq', 'phred', 'sampleId'], table_prefix=None ):
        ''' Return the cluster object for the given id '''
    
        if table_prefix is None:
            members_table = 'members'
            clusters_table = 'clusters'
        else:
            members_table = table_prefix + '_members'
            clusters_table = table_prefix+ '_clusters'
    
        # Set which parameters to get 
        get_seq = 0
        get_phred = 0
        get_sampleid = 0
        if 'seq' in items:
            get_seq =  True
        if 'phred' in items:
            get_phred = True
        if 'sampleId' in items:
            get_sampleid = True
    
        with self.con as con:
            
            repseq_sql_query = ''' SELECT repseqid, size, clusterid FROM {clusters} 
                WHERE clusterId = ? '''.format(clusters=clusters_table)
             
            members_sql_query = ''' SELECT {items}  FROM seqs 
                JOIN {members} USING (seqId) JOIN {clusters} USING (clusterId)
                WHERE clusterId = ?'''.format(items = ','.join(items), 
                        members=members_table, clusters=clusters_table) 
    
            clusterobj = ClusterObj()
            
            #===================================================================
            # Get Repseq Info   
            #===================================================================
            t1 = time.time()
            c = con.execute(repseq_sql_query, (cluster_id,))
            print 'Time for repseq_sql_query: ', time.strftime('%H:%M:%S', time.gmtime(time.time() - t1))
            
            cluster_row = c.fetchone()
            
            clusterobj.rep_seq_id = cluster_row['repseqid']
            clusterobj.size = cluster_row['size']
            clusterobj.id = cluster_row['clusterId']
            
            # get repseq_seq
            t1 = time.time()
            c = con.execute(''' SELECT {items} FROM seqs WHERE seqId = ?'''.format(items = ','.join(items)), (clusterobj.rep_seq_id,))
            print 'Time for seq_sql_query: ', time.strftime('%H:%M:%S', time.gmtime(time.time() - t1))

            row = c.fetchone()
            
            if get_seq: 
                clusterobj.rep_seq = row['seq']
            if get_phred:
                phred_ascii = row['phred']
                phred_list = [ord(c) - 33 for c in phred_ascii]
                clusterobj.rep_phred = np.array(phred_list)
            if get_sampleid:
                    clusterobj.rep_sample_id = row['sampleId']
            
            #===================================================================
            # Get Members Info
            #===================================================================
            t1 = time.time()
            curs = con.execute(members_sql_query, (cluster_id,))
            print 'Time for members_sql_query: ', time.strftime('%H:%M:%S', time.gmtime(time.time() - t1))

            for row in curs:
            
                if get_seq: 
                    clusterobj.members_id.append(row['seqId'])
                    clusterobj.members_seq.append(row['seq'])
                if get_phred:
                    phred_ascii = row['phred']
                    phred_list = [ord(c) - 33 for c in phred_ascii]
                    clusterobj.members_phred.append(np.array(phred_list))
                if get_sampleid:
                    clusterobj.members_sample_id.append(row['sampleId'])
                    
        return clusterobj          
    
    def write_reads(self, out_file_handle, output_format='fasta', filter_expression=None, 
                     startidx=0, rowbuffer=100000, overwrite=False):
        """ Write records returned by the querry to one large fasta or fastq 
        
        Defaults is to search by GLOBing the individual descriptions with the search_query.
        
            If sql_query = True, search_query is passed as a full sql statment.
            If use_type_column=True, search is done by GLOBing the individual type column instead.
        
        
        out_file_handle -- A file object or string specifying a filename. 
        
        startidx -- starting base index of DNA sequence that is written, used to miss out 
        cutsite if necessary.        
        """
        
        # Output check
        out_file_handle = outputfile_check(out_file_handle, mode='a', overwrite=overwrite)
         
        query = '''SELECT seqid, seq, phred  
                    FROM seqs INNER JOIN samples ON seqs.sampleId=samples.sampleId'''
            
        if filter_expression:
            query +=  ' WHERE {0}'.format(filter_expression)
        
        with self.con as con:
            
            toc = time.time()
            print >> sys.stderr, 'Writing records to {0} format....'.format(output_format), 

            c = con.execute(query)
            returned_records = c.fetchmany(rowbuffer)
            rec_count = 0
            while returned_records:
                
                for rec in returned_records:
                    seq_rec = SeqRecord(Seq(rec['seq'][startidx:]), id=str(rec['seqid']), description='')
                    
                    if output_format=='fastq':
                        seq_rec.letter_annotations['phred_quality'] = [ ord(x)-33 for x in rec['phred']]
                    
                    SeqIO.write(seq_rec, out_file_handle, format=output_format)
                    rec_count += 1 
                
                # Fetch next batch of records from cursor          
                returned_records = c.fetchmany(rowbuffer)
            
            print >> sys.stderr, ' Done!'
            print >> sys.stderr, '\n{0} records written successfully to {1}\nin {2}'.format(rec_count,
                out_file_handle.name, time.strftime('%H:%M:%S', time.gmtime(time.time() - toc))) 
            
            if out_file_handle.name not in ['<stdout>', '<stderr>']:
                out_file_handle.close()
                
        return out_file_handle 
    
    def calculate_reads_per_individual(self, individuals=None):
        """ Fill in individuals table with total reads of each 

        """
        
        # Get list of individuals
        with self.con as con :
        
            if individuals is None:    
                c = con.execute(''' SELECT sampleId from samples ''')
                rows = c.fetchall()
            
                num = len(rows)
                print >> sys.stderr, 'Updating read counts for {0} samples'.format(num),
                
                for r in rows:
                    # for each, find the total number of reads
                    c = con.execute(''' SELECT count(*) FROM seqs WHERE sampleId = ?''',  (r['sampleId'],))
                    
                    total = c.fetchone()['count(*)'] 
                    
                    # Fill in table  
                    con.execute(''' UPDATE samples SET read_count = ? WHERE sampleId = ? ''', (total , r['sampleId'] ))
            else:
                assert type(individuals) is list, 'Invalid parameter. Expected a list of sample ids'
                rows = individuals

                for r in rows:
                    # for each, find the total number of reads
                    c = con.execute(''' SELECT count(*) FROM seqs WHERE sampleId = ?''',  (r,))
                    total = c.fetchone()['count(*)'] 
                    
                    # Fill in table  
                    con.execute(''' UPDATE samples SET read_count = ? WHERE sampleId = ? ''', (total , r))
                    
            # Build index on Cluster size 
            con.execute('CREATE INDEX IF NOT EXISTS read_countIndex ON samples(read_count)')
    
    def load_cluster_file(self, cluster_file_handle, table_prefix=None, 
                          overwrite=False, fmin=2, fmax=None, skipsort=False, buffer_max=1000000):
        ''' Load in a clustering file into the database 
        
        By default singletons are not loaded as cutoff = 2
        
        can also set a fmin and fmax threshold if only clusters of a certain 
        size are to be added.
        
        '''
        # define names
        if table_prefix is None:
            members_table_name = 'members'
            cluster_table_name = 'clusters'
            index_name = 'clustersizeIndex'
        else:
            members_table_name = table_prefix + '_members'
            cluster_table_name = table_prefix + '_clusters'
            index_name = table_prefix + '_clustersizeIndex'
        
        # input checks
        if type(cluster_file_handle) == str:
            if not cluster_file_handle.endswith('.clstr'):
                cluster_file_handle = cluster_file_handle + '.clstr'
        
        cluster_file_handle = inputfile_check(cluster_file_handle)
        
        # Sort file
        if not skipsort:
            # Filter out singletons and sort clusters in accending order
            print >> sys.stderr, 'Sorting cluster file %s ...' % (cluster_file_handle.name)
            sorted_cluster_file = sortby(cluster_file_handle, reverse=True, 
                                         mode='cluster_size', outfile_postfix='-subset', cutoff=fmin)
            cluster_file_handle = inputfile_check(sorted_cluster_file)

        print >> sys.stderr, 'Importing cluster file %s  to database...' % (cluster_file_handle.name)

        # Overwrite/ make tables if necessary
        if overwrite:
            with self.con as con:
                con.execute(''' DROP TABLE IF EXISTS {0} '''.format(cluster_table_name))
                con.execute(''' DROP TABLE IF EXISTS {0} '''.format(members_table_name))
        
        self.create_cluster_table(cluster_table_name)
        self.create_members_table(members_table_name)
        
        # Make cluster generator. Returns all cluster info
        cluster_gen = parse(cluster_file_handle, similarity_count=True, edit_dist=True)
        
        # Buffer to hold clusters in memory then write all at once
        cluster_info_list = []
        cumulative_cluster_size = 0 
        
        # Find starting cluster id 
        c = self.con.execute(''' SELECT COUNT(*) FROM {0}'''.format(cluster_table_name))
        clusterid = c.fetchone()['count(*)'] + 1
        
        # Dorp any index on Cluster size 
        with self.con as con:
            con.execute(''' DROP INDEX IF EXISTS {0} '''.format(index_name))
        
        if fmax:
            for cluster in cluster_gen:
                if cluster.size <= fmax and cluster.size >= fmin:
                            
                    cluster_info_list.append( ( clusterid, cluster.rep_seq_id, cluster.size, cluster.members_id)  )
                 
                    clusterid += 1 
                    cumulative_cluster_size += cluster.size
                    
                    if cumulative_cluster_size > buffer_max:
                        self.load_batch_clusterdata(cluster_info_list, table_prefix)
                        cluster_info_list = []
        else:    
            for cluster in cluster_gen:
                cluster_info_list.append( ( clusterid, cluster.rep_seq_id, cluster.size, cluster.members_id)  )
                clusterid += 1 
                cumulative_cluster_size += cluster.size
                
                if cumulative_cluster_size > buffer_max:
                    self.load_batch_clusterdata(cluster_info_list, table_prefix)
                    cluster_info_list = []
        
        # Final flush of data
        if cluster_info_list:
            self.load_batch_clusterdata(cluster_info_list, table_prefix)
            
        # Rebuild index on Cluster size 
        with self.con as con:
            con.execute('''CREATE INDEX {indexname} ON {tablename}(size)'''.format(
                    indexname=index_name, tablename=cluster_table_name))
            
    def load_batch_clusterdata(self, data_structure, exp_name=None):
        ''' Load in a single cluster object to the database. '''
        
        if exp_name is None:
            members_table_name = 'members'
            cluster_table_name = 'clusters'
        else:
            members_table_name = exp_name + '_members'
            cluster_table_name = exp_name + '_clusters'
              
        with self.con as con:
            
            # Just insert cluster into clusters table and incriment id 
            curs = con.cursor()
            
            g = ((x[0], x[1], x[2]) for x in data_structure)
            
            curs.executemany(''' INSERT INTO {0}
                    (clusterid, repseqid, size) VALUES (?,?,?) '''.format(cluster_table_name), g )
                
            for x in data_structure:
                g = ( ( x[0], int(seqid) ) for seqid in x[3])
                
                # Processing members using executemany is much faster for many inserts
                con.executemany(''' INSERT INTO {0}
                     (clusterId, seqId) VALUES (?,?) '''.format(members_table_name), g)

if __name__ == '__main__':
    pass