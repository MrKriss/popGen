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
                    
                    else:
                        # Start a new seed
        
        
        
        
        
        # Run 
        
        
    class ClusterDictionary(dict):


        
class odict(dict):
         def __init__(self, *args, **kw):
             super(odict,self).__init__(*args, **kw)
             self.itemlist = super(odict,self).keys()
         def __setitem__(self, key, value):
              # TODO: what should happen to the order if
              #       the key is already in the dict       
             self.itemlist.append(key)
             super(odict,self).__setitem__(key, value)
         def __iter__(self):
             return iter(self.itemlist)
         def keys(self):
             return self.itemlist
         def values(self):
             return [self[key] for key in self]  
         def itervalues(self):
             return (self[key] for key in self)







if __name__ == '__main__':
    pass