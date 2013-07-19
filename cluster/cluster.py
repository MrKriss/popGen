'''
Created on 17 Jul 2013

@author: musselle
'''
"""
Conduct fast clustering on a database of reads
"""
import os 
import sys 
import socket

import numpy as np 
import matplotlib.pyplot as plt

from editdist import distance

import gzip
import argparse


from Bio import SeqIO, bgzf 
from collections import Counter

# Work out where data is stored
if socket.gethostname() == 'yildun':
    data_prefix = '/space/musselle/data'
    work_prefix = '/home/pgrad/musselle/ubuntu/workspace/'
        
elif socket.gethostname() == 'luca':
    data_prefix = '/home/musselle/san/data'
    work_prefix = '/home/musselle/'
       
elif socket.gethostname() == 'gg-pc6':
    data_prefix = '/home/musselle/data'
    work_prefix = '/home/musselle/'

sys.path.append(os.path.join(work_prefix, 'popGen'))
sys.path.append(os.path.join(work_prefix, 'popGen/utils'))
sys.path.append(os.path.join(work_prefix, 'popGen/plotscripts'))
sys.path.append(os.path.join(work_prefix, 'popGen/runscripts'))
sys.path.append(os.path.join(work_prefix, 'popGen/clf'))

del work_prefix
del data_prefix

from database.reads_db import Reads_db


class ClusterReadsObj():
    
    
    def __init__(self, args, db, records):
        self.args = args
        
        self.edit_threshold = args.similarity
        
    
    def run(self, record_curs, args, db):
        """" Run the clustering process. """

        # Setup data structures 
        repseq_dtype = np.dtype([('repid', 'u4'), 
                                 ('repseq', 'a88')]) 
        members_dtype = np.dtype(np.uint32)

        cluster_dict = {}

        max_cluster_id = 0 
        
        
        returned_records = True
        
        while returned_records is not None:
            
            returned_records = record_curs.fetchmany()
        
            for read in returned_records:
                
                if not cluster_dict:
                    # First iteration: set first read as a seed and initialise cluster structure
                    
                    cluster_dict[max_cluster_id] = [np.array(1, dtype=repseq_dtype), 
                                                    np.array(1, dtype=members_dtype)]
                    
                    max_cluster_id += 1
                else:
                    # go through all clusters and work out distance between read and seeds 
                    dist_vec = np.zeros(max_cluster_id, dtype=int)
                    
                    for i in range(max_cluster_id):
                        dist_vec[i] = distance(cluster_dict[i][0][1] , read['seq'])
                    
                    first_min_idx = dist_vec.argmin()
                    
                    
                    len_seq = len(read['seq'])
                    
                    
                    if ((len_seq - dist_vec[first_min_idx]) / len_seq) >= args.similarity:
                        # Add read to the membership of that cluster 
                        pass
                    else:
                        pass
                        # Start a new seed
        
        
        
        
class ClusterDictionary(dict):
    def __init__(self, *args, **kw):
        
        if 'seqlen' in kw:
            seqlen = 'a' + str(kw['seqlen'])
            self.repseq_dtype = np.dtype([('repid', 'u4'), 
                                          ('repseq', seqlen)]) 
        else:
            self.repseq_dtype = np.dtype([('repid', 'u4'), 
                                          ('repseq', 'a88')]) 
            
        self.members_dtype = np.dtype(np.uint32)
        self.aprox_memsize = 0
        
        super(ClusterDictionary,self).__init__(*args, **kw)
        
        self.max_cluster_id = len(super(ClusterDictionary,self).keys())
        
#     def __setitem__(self, key, value):
#         super(ClusterDictionary,self).__setitem__(key, value)
#     
#     def __iter__(self):
#         return super(ClusterDictionary,self).__iter__()
    
    def add_cluster(self, repid, repseq):
        ''' Add a new cluster with the given sequence and id as the representative sequence '''
        self[self.max_cluster_id] = [np.array([( repid , repseq )], dtype=self.repseq_dtype), 
                                     np.zeros(0, dtype = self.members_dtype)]
        self.aprox_memsize += self[self.max_cluster_id][0].nbytes
        self.max_cluster_id += 1
        
        
    def add_member(self, clusterid, memberids):
        num_members = len(memberids)
        self[clusterid][1] = np.append(self[clusterid][1], np.array([memberids], dtype=self.members_dtype))
        self.aprox_memsize += (4 * num_members)

        
    def flush2db(self, db, flush_thresh = 1000,  seq_table_name='seqs',  clust_table_name='clusters'):
        ''' Update/insert all clusters into database with members. Then flush dictionary '''
        
        
        processing_msg = ['10%', '25%', '...', '....', '.....']
        
        # Make table if not present
        if clust_table_name not in db.tables:
            db.create_cluster_table(clust_table_name, overwrite=True)
            
        with db.con as con:
            # Update Clusterids
            for clustid in range(self.max_cluster_id):
                
                members_size = len(self[clustid][1])    
                
                if members_size < flush_thresh:
                    continue
                
                if not clustid % 1000:
                    print 'done 1000'

                repseqid = self[clustid][0][0]['repid']
                
                records = con.execute(''' SELECT size FROM {0} WHERE clusterId = ?'''.format(clust_table_name), (clustid,) )
                row = records.fetchone()
    
                # Process clusterid 
                if row is None:
                    # Cluster id does not yet exist, Insert new cluster
                    clust_size = members_size + 1

                    con.execute(''' INSERT INTO {0}
                    (clusterid, repseqid, size) VALUES (?,?,?)'''.format(clust_table_name), 
                    ( clustid , int(repseqid) , clust_size ) )
                else: 
                    # Cluster id exists, update cluster size.
                    con.execute(''' UPDATE {0} SET
                        size = ? WHERE clusterid = ?'''.format(clust_table_name), 
                        ( row['size'] + members_size, clustid))
                    
                # Process members
                for member in self[clustid][1]:
                    con.execute(''' UPDATE {0} SET
                        clusterId = ? WHERE seqId = ?'''.format(seq_table_name), 
                        ( clustid, int(member)))
                
                # Clear members
                self[clustid][1] = np.zeros(0, dtype = self.members_dtype)
                self.aprox_memsize -= (members_size * 4)
            
#     def keys(self):
#         return self.itemlist
#     def values(self):
#         return [self[key] for key in self]  
#     def itervalues(self):
#         return (self[key] for key in self)







if __name__ == '__main__':
    
    # Parse arguments

    parser = argparse.ArgumentParser(description='Run Simple clustering on reads in the database.')
    parser.add_argument('-i',  dest='input', required=True, nargs='+',
                        help='Database file where reads are stored (/path/filename)')
    parser.add_argument('-q',  dest='query', required=True, nargs='+',
                        help='Querry to fetch records with.')
    parser.add_argument('-s',  dest='similarity', required=True, nargs='+',
                        help='Threshold for percentage similarity between clusters.')
    
    args = parser.parse_args()
    
    
    # Fetch records 
    db = Reads_db(args.input, recbyname=True)
    
    records = db.con.execute(args.query)
    
    
    # Run clustering 
    
    
    
    # Store results 

    
    