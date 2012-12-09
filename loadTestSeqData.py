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
from utils.utils import Cycler

if __name__ == '__main__':
    
    dataLoc = '/space/musselle/datasets/gazellesAndZebras/lane6'
    RecCycl = Cycler(filepattern ='*[0-9].fastq.bgzf', datapath = dataLoc, maxNumSeq = 2)
    
    for rec in RecCycl.recGen:
        print rec
    