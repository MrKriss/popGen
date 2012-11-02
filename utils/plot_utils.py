'''
Created on 24 Oct 2012

@author: musselle
'''

import os 

import numpy as np 
import matplotlib.pyplot as plt 


def makeHist(files , dataPath, bins = 50):
    ''' Function to plot histograms of stored data files'''

    os.chdir(dataPath)
 
    data = [0] * len(files)
    
    for i, f in enumerate(files):
        
        name = 'f' + str(i)
        vars()[name] = np.load(f)
        data[i] = vars()[name] 

    plt.figure()
    plt.hist(data, bins = bins)
    plt.ylabel('Frequency')


if __name__ == '__main__' : 
    

    dataPath = '/space/musselle/datasets/gazellesAndZebras'
    fs1 = ['lane6_meanPhred.npy', 'lane8_meanPhred.npy']
    fs2 = ['lane6_propN.npy', 'lane8_propN.npy']
    
    makeHist(fs1, dataPath)
    plt.title('Mean Phred Score')
    makeHist(fs2, dataPath)
    plt.title('Proportion of Ns')
