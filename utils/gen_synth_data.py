'''
Created on Feb 26, 2013

@author: chris
'''
import os 
import sys
import numpy as np

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

import random
import string 

# Random ID generator
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

# Random sequence generator 
def seq_generator(length, chars ='ATGC'):
    return ''.join(random.choice(chars) for x in range(length))

def mid_generator(mid_lib):
    pass
    
def cause_SNP(seq):
    pass

def cause_indel(seq):
    pass

def add_error():
    pass

def f():
        '''   '''
        pass

# Example
simple_seq = Seq("GATC")
simple_seq_r = SeqRecord(simple_seq, id="AC12345")

# Vars 
#(length of seq, Number of times it appears)
seq_len_and_num_per_cluster = [(100, 1000), (50, 2000), (110, 3000), (150, 500), (90, 1000), (90, 1000), (90, 1000)]

out_filename = 'test_seqs.fasta'
outpath = None

if not out_filename.endswith('.fasta'):
    out_filename = out_filename + '.fasta'

if outpath is None:
    outpath = os.getcwd()

target_file_path = os.path.join(outpath, out_filename)

# Check file does not already exist
if os.path.isfile(target_file_path):
    os.remove(target_file_path)


cluster_count = 0 

with open(target_file_path, 'a') as out_filehdl:

    for tup in seq_len_and_num_per_cluster:
        
        cluster_count += 1
        new_seq_cluster = True
        
        seq = Seq(seq_generator(tup[0], 'ATGC'))
        
        for num_seq in range(tup[1]):
            # append a copy of the random seq to the file with unique id 
            # ID also holds the cluster number it belongs too 
            id = '.'.join([id_generator(4), str(cluster_count), str(num_seq)])
            
            # Any errors to add, do so here.
            
            seq_rec = SeqRecord(seq, id=id)
            SeqIO.write(seq_rec, out_filehdl, 'fasta')
            
if __name__ == '__main__':
    pass
    