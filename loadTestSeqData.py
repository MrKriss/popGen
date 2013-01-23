'''
Created on 21 Nov 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
from Bio import SeqIO
import glob


from utils import Cycler

if __name__ == '__main__':
    
    dataLoc = '/space/musselle/datasets/gazellesAndZebras/testdata'
    RecCycl = Cycler(filepattern ='*[0-9].fastq.bgzf', inpath = dataLoc, maxnumseq = 5)
    
    for rec in RecCycl.recgen:
        print rec.seq
        a = np.array(rec.letter_annotations['phred_quality'])
        print rec.letter_annotations['phred_quality']
        print 'Mean Phred = %d' % a.mean(), ' Filter = %s' % rec.description.split()[1].split(':')[1] == 'Y'
        
        
        
    