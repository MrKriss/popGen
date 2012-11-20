'''
Created on 19 Nov 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

from basicAnalysis import Cycler


def filterReads(inFiles = None, fileType = '', dataPath = ''):
    ''' Filter reads based on criteria 
    
    Default is to use Machine Specific read filter 
    
    '''
    
    RecCycler = Cycler(inFiles = None, fileType = '', dataPath = '')
    
    for s in RecCycler.seqGen:
        
        if s.description.split()[1].split(':')[1] == 'N': 
    
    
    



if __name__ == '__main__':
    pass