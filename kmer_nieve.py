'''
Created on 3 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt


from collections import Counter

seq = 'AGATAGATAGACACAGAAATGGGACCACAC'
kmers = {}
k = 4

for i in range(len(seq) - k + 1):
   kmer = seq[i:i+k]
   if kmers.has_key(kmer):
      kmers[kmer] += 1
   else:
      kmers[kmer] = 1

for kmer, count in kmers.items():
   print kmer + "\t" + str(count)


if __name__ == '__main__':
    pass