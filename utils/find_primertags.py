'''
Created on 14 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

import glob
from utils import Cycler
from utils import pklsave
import time

def tags_counter(infiles=None, sl=(6,12), filepattern=False, datapath=''):
    '''Writes the reads (without the MID tag) to one large fasta file 
    for clustering'''
    
    from collections import Counter
    
    TagCounter = Counter()
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, datapath=datapath)
    
    c = 0
    
    toc = time.time()
    cum_t = 0
    for rec in RecCycler.recgen:
    
        primer_tag = rec.seq[sl[0]: sl[1]].tostring()
        TagCounter[primer_tag] += 1
        c+=1
        if not c % 10000000:
            loop_t = time.time() - toc - cum_t
            cum_t += loop_t
            print 'Processed 10 million reads after {1}'.format(
                time.strftime('%H:%M:%S', time.gmtime(loop_t)))
          
    total_t = time.time() - toc
    print '\nProcessed {0} files in {1}'.format(
            time.strftime('%H:%M:%S', time.gmtime(total_t))) 
        
    return TagCounter


if __name__ == '__main__':
     
    #===========================================================================
    ''' RUNS SCRIPT FOR ALLL READS IN LANE 6 '''
    #===========================================================================
    
    LANE = '8'
    
    # Set paths and file patterns 
#    datapath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
    datapath = '/home/musselle/data/lane' + LANE
    os.chdir(datapath)
    raw_files = glob.glob('*[0-9].fastq.bgzf')
    raw_files.sort()
    
    TagsCounter = tags_counter(infiles = raw_files, sl=(6,12))
    pklsave(TagsCounter, 'L{0}_TagsCount'.format(LANE))





    