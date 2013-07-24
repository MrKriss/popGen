#! /space/musselle/EPD7-3/bin/python 

'''
Created on 6 Jun 2013

@author: musselle

Read a CDHIT cluster file from stdin and filter based on parsed parameters

write those clusters that  

'''
import sys
import os
import socket

# Work out where data is stored
if socket.gethostname() == 'yildun':
    data_prefix = '/space/musselle/data'
    work_prefix = '/home/pgrad/musselle/ubuntu/workspace/'
        
elif socket.gethostname() == 'luca':
    data_prefix = '/home/musselle/san/data'
    work_prefix = '/home/indexfilemusselle/'
       
elif socket.gethostname() == 'gg-pc6':
    data_prefix = '/home/musselle/data'
    work_prefix = '/home/musselle/'

sys.path.append(os.path.join(work_prefix, 'popGen'))
sys.path.append(os.path.join(work_prefix, 'popGen/utils'))
sys.path.append(os.path.join(work_prefix, 'popGen/plotscripts'))
sys.path.append(os.path.join(work_prefix, 'popGen/runscripts'))
sys.path.append(os.path.join(work_prefix, 'popGen/clf'))

import argparse
from collections import Counter

from utils.ClusterIO import sortby, parse, Cluster
from Bio import SeqIO

# Setup parser
parser = argparse.ArgumentParser(description='Procedure to get the genuine representative sequences for each clusters in a CDHIT.clstr file')
parser.add_argument('infile', action='store', help='Path to the input CD-HIT file for all reads')
parser.add_argument('indexfile', action='store', help='Path to the index file for all reads')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()

# Run through all Clusters
cluster_gen = parse(args.infile, idx_file_path=args.indexfile, edit_dist=False, similarity_count=True)

# Slow, Only do this once.
seqrec_lookup = SeqIO.index_db(args.indexfile)

rep_seq_correction_count = 0
clusters_count = 0

with open('clusters.temp', 'wb') as clust_file:

    for cluster in cluster_gen:
        # if fraction of reads 100% similar to representative seq is in majority,
        # skip and write to new file
        seq_similarity_ranking = cluster.similarity_counter.most_common()
        
        print seq_similarity_ranking[0][0]
        
        if seq_similarity_ranking[0][0] != '100.00':
            # Work out which sequence is most representative 
            cluster.get_unique_seq(ignoreup2=6, db=seqrec_lookup)
            most_common_seq = cluster.unique_seq.most_common()[0][0]
            
            # Find the id of first sequence with the most common sequence and update rep seq
            for idx, mem in enumerate(cluster.members_seq):
                
                if mem[6:] == most_common_seq:
                    
                    # store temps
                    old_rep_seq = cluster.rep_seq
                    old_rep_seq_desc = cluster.rep_seq_id
                    
                    # update new rep seq
                    cluster.rep_seq = mem
                    cluster.rep_seq_id = cluster.members_id[idx]

                    # Remove from members and add old rep seq
                    del cluster.members_id[idx]
                    del cluster.members_seq[idx]
                    cluster.members_id.append(old_rep_seq_desc)
                    cluster.members_seq.append(old_rep_seq)
                    
                    rep_seq_correction_count += 1
                    break
            
        # Write new .clstr file
        cluster.write2clstr(clust_file, db=seqrec_lookup)
        clusters_count += 1
        
    print 'Wrote {0} Clusters and corrected {1} ({2:.2f}%)'.format(clusters_count, rep_seq_correction_count, 
                                                            100 * float(rep_seq_correction_count)/clusters_count)    
                     
# # Sort file 
# sortby(args.filename, reverse=True, mode='reads_per_cluster', outfile_postfix='-filtered', cutoff=args.minreads, 
#        clustersize_min=args.min, clustersize_max=args.max)
# 
# # Rewrite file
# old_filename_parts = args.filename.name.split('.')
# old_filename_parts[0]  += '-sortedby_' + 'reads_per_cluster'
# new_filename = '.'.join(old_filename_parts)
