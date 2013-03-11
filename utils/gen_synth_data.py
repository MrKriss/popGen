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


# Simple test for indels and error rates on CD-HIT clustering

# Number of randomly generated reads (num cutsites x 2) (25,000 * 2)
# Poisson Depth of coverage



# Option to Diplodise and add alleles 

# Number of individuals 
# number of sites covered per in







# Stacks methodology
#-------------------
# Get reference sequence
# Extract reads in both directions from each cut site (CCTGCA^GG)
# Re-diplodise the genome by creating alleles (2 copies of each read)
# Introduce SNPs in these alleles at a rate of 0.5%
# 'Sequence' to a mean depth determined by drawing from a poisson distribution with means of (10, 20, 30)
# Simulated Sequencing error rates which increased linearly along the illumina read (means of 0.5%, 1%, and 3%)
# Each run on simulated data involved 10 replicates

filepath = '/space/musselle/reference-genomes/Sticklebacks/Gasterosteus_aculeatus.BROADS1.56.dna.toplevel.fa'
 
def synthdata_fromgenome(genome_file, cutsite = 'CCTGCAGG'):
    ''' Generate synthetic 'Reads' from a reference genome ''''
    
    # Sort reference genome in order 
    
    D = SeqIO.index(genome_file, 'fasta')
    .keys
    
    # Find all cut sites in the reference genome
    
    # Cycle through all scaffolds























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
    