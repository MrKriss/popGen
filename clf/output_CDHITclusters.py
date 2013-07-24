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

import argparse
from collections import Counter

from utils.ClusterIO import sortby, parse
from Bio import SeqIO

# Setup parser
parser = argparse.ArgumentParser(description='Procedure to write FastQ output file for clusters in a CDHIT.clstr file')
parser.add_argument('indexfile', action='store', help='Path to the index file for all reads')
parser.add_argument('ouputdir', action='store', help='Directory where outputs are stored to the index file for all reads')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()

if not os.path.exists(args.ouputdir):
    os.makedirs(args.ouputdir)

db = SeqIO.index_db(args.indexfilepath)

idx_file_dir = os.path.split(args.indexfilepath)[0]

cluster_gen = parse(sys.stdin, idx_file_path=args.indexfile, edit_dist=False)

size_counter = Counter()
 
count = 0
 
for cluster in cluster_gen:

    count += 1 
    
    # Write to output file
    size_counter[str(cluster.size)] += 1
    
    # Get the sequence records for the cluster 
    seqs = []
    if os.getcwd() != idx_file_dir:
        os.chdir(idx_file_dir)
    seqs.append(db[cluster.rep_seq_id])
    for member in cluster.members_id:
        seqs.append(db[member])
    
    if os.getcwd() != self.output_dirpath:
        os.chdir(self.output_dirpath)
        
    # Write cluster to a file 
    fname = "clustersize{0}-No{1}.fastq".format(str(cluster.size), str(size_counter[str(cluster.size)]))
    output_handle = open(fname, "wb")
    SeqIO.write(seqs, output_handle, "fastq")
        
    elif self.reads_per_cluster_size_counter[str(cluster.size)] < self.filter_params['min_reads']:
        break













# Sort file 
sortby(args.filename, reverse=True, mode='reads_per_cluster', outfile_postfix='-filtered', cutoff=args.minreads, 
       clustersize_min=args.min, clustersize_max=args.max)

# Rewrite file
old_filename_parts = args.filename.name.split('.')
old_filename_parts[0]  += '-sortedby_' + 'reads_per_cluster'
new_filename = '.'.join(old_filename_parts)

