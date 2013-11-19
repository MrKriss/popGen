#! /usr/bin/env python 

'''
Created on 17 Jul 2013

@author: musselle
'''

import sys
import time
import argparse
from editdist import distance

import numpy as np

from rapier.lib.reads_db import Reads_db


class ClusterReadsObj():
    
    def __init__(self, args):
        self.args = args
        
        self.edit_threshold = args.similarity

        self.cluster_dict = ClusterDictionary()
        
        self.record_count = 0
    
    def run(self, record_curs, db):
        """" Run the clustering process. """

        for read in record_curs:
            
            self.record_count += 1
            
            if not self.cluster_dict:
                # First iteration: set first read as a seed and initialise cluster structure
                self.cluster_dict.add_cluster(read['seqid'], read['seq'])

            else:
                # go through all clusters and work out distance between read and seeds 
                dist_vec = np.zeros(self.cluster_dict.max_cluster_id, dtype=int)
                
                for i in range(self.cluster_dict.max_cluster_id):
                    dist_vec[i] = distance(self.cluster_dict[i][0][0]['repseq'] , read['seq'])
                
                first_min_idx = dist_vec.argmin()
                
                len_seq = len(read['seq'])
                
                if ((len_seq - dist_vec[first_min_idx]) / len_seq) >= args.similarity:
                    # Add read to the membership of that cluster 
                    self.cluster_dict.add_member(first_min_idx, read['seqid'])
                else:
                    # Start a new seed
                    self.cluster_dict.add_cluster(read['seqid'], read['seq'])
        
            # Check at the end of every batch of records from fetchmany()
            if args.maxmemory != 0:
                if self.record_count % 100 == 0:
                    if self.cluster_dict.aprox_memsize > self.args.maxmemory:
                        self.cluster_dict.flush2db(db)
        
        print 'Finished Clustering!\nAbout to flush all remaining records to database'
        self.cluster_dict.flush2db(db, flush_thresh=-1)
            
        
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

        
    def flush2db(self, db, flush_thresh = 100,  seq_table_name='seqs',  clust_table_name='clusters'):
        ''' Update/insert all clusters into database with members. Then flush dictionary '''
        
        # Make table if not present
        if clust_table_name not in db.tables:
            db.create_cluster_table(clust_table_name, overwrite=True)
            
        with db.con as con:
            # Update Clusterids
            
            # progress msg
            increment = self.max_cluster_id * .25
            x = 0
            
            for clustid in range(self.max_cluster_id):
                
                members_size = len(self[clustid][1])    
                
                if members_size < flush_thresh:
                    continue
                
                # Update status msg
                if clustid >= increment * x:
                    print '\r{0}% complete'.format(str(25*x)),
                    x +=1

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

    toc = time.time()

    parser = argparse.ArgumentParser(description='Run CDHIT clustering on reads in the database.')
    parser.add_argument('-i',  dest='input', required=True,
                        help='Database file where reads are stored (/path/filename)')
    parser.add_argument('-q',  dest='query', default=None,
                        help='Querry to fetch records with. Default will cycle through all records in database')
    parser.add_argument('-s',  dest='similarity', required=True,
                        help='Threshold for percentage similarity between clusters.')
    parser.add_argument('-m',  dest='maxmemory', default=0,
                        help='Maximum memory to use for the clustering. Default = 0 (unlimited)')

    args = parser.parse_args()
    
    # Fetch records 
    db = Reads_db(args.input, recbyname=True)

    if args.query is None:
        args.query = '''SELECT * FROM seqs'''    

    record_curs = db.con.execute(args.query)
    
    # Run clustering 
    clustering =  ClusterReadsObj(args)
    clustering.run(record_curs, db)
    
    total_t = time.time() - toc    
    print >> sys.stderr, 'Clustered {0} records into {1} clusters in {2}'.format(
              clustering.record_count, len(clustering.cluster_dict),
              time.strftime('%H:%M:%S', time.gmtime(total_t)))
    
    