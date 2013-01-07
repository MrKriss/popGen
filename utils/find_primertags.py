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

def tags_counter(infiles=None, slice=(6,12), filepattern=False, datapath='', 
                outfile='outfile.fasta'):
    '''Writes the reads (without the MID tag) to one large fasta file 
    for clustering'''
    
    from collections import Counter
    
    TagCounter = Counter()
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, datapath=datapath)
    
    for rec in RecCycler.recgen:
    
        primer_tag = rec.seq[slice[0]: slice[1]].tostring()
        TagCounter[primer_tag] += 1
        
    return TagCounter


if __name__ == '__main__':
     
    #===========================================================================
    ''' RUNS SCRIPT FOR ALLL READS IN LANE 6 '''
    #===========================================================================
    
    LANE = '8'
    
    # Set paths and file patterns 
    datapath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    os.chdir(datapath)
    raw_files = glob.glob('*[0-9].fastq.bgzf')
    raw_files.sort()
    
    TagsCounter = tags_counter(infiles = raw_files, slice=(6,12))
    pklsave(TagsCounter, 'L{0}_TagsCount'.format(LANE))





    