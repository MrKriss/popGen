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


def primertags_counter(infiles=None, filepattern=False, datapath='', 
                outfile='outfile.fasta'):
    '''Writes the reads (without the MID tag) to one large fasta file 
    for clustering'''
    
    from collections import Counter
    
    TagCounter = Counter()
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, datapath=datapath)
    
    for rec in RecCycler.recgen:
    
        primer_tag = rec.seq[6:12].tostring()
        TagCounter[primer_tag] += 1
        
    return TagCounter


if __name__ == '__main__':
     
    #===========================================================================
    ''' RUNS SCRIPT FOR ALLL READS IN LANE 6 '''
    #===========================================================================
    
    LANE = '6'
    
    # Set paths and file patterns 
    datapath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    os.chdir(datapath)
    raw_files = glob.glob('*[0-9].fastq.bgzf')
    raw_files.sort()
    
    TagsCounter = primertags_counter(infiles = raw_files)





    