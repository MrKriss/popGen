'''
Created on 11 Mar 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

from gen_synth_data import seq_generator
from cluster import cluster_cdhit

f = open('indel_test.fasta', 'w')

seq1 = Seq(seq_generator(100))

seq2 = seq1[1:] + 'T'
seq3 = seq1[:50] + seq1[51:] + 'G'
seq4 = seq1[:25] + 'T' + seq1[25:-1]
seq5 = seq1[:75] + 'A' + seq1[75:-1]

seqR1 = SeqRecord(seq1, id='original')
seqR2 = SeqRecord(seq2, id='deletion0')
seqR3 = SeqRecord(seq3, id='deletion50')
seqR4 = SeqRecord(seq4, id='insertion25')
seqR5 = SeqRecord(seq5, id='insertion75')

SeqIO.write([seqR1, seqR2, seqR3, seqR4, seqR5], f, 'fasta')
f.flush()
f.close()

cluster_cdhit('indel_test.fasta','indel_test_clust', 0.9, 8, maskN=0)

with open('indel_test_clust.clstr', 'r') as outf:
    for line in outf:
        print line