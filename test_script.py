'''
Created on 7 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
from utils import make_MIDdict
import editdist as ed

if __name__ == '__main__':
    
    barcode_path = '/space/musselle/datasets/gazellesAndZebras/barcodes/'
    tags = make_MIDdict('*[6].txt', filepattern=True, inpath=barcode_path)

    distmat = np.zeros((len(tags), len(tags)))
    MIDs = [key[:6] for key in tags.iterkeys()]
    MIDs.sort()

    for i in range(distmat.shape[0]):
        for j in range(distmat.shape[1]):
            distmat[i,j] = ed.distance(MIDs[i], MIDs[j])
        
        


