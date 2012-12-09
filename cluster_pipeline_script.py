'''
Created on 30 Nov 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

from filter import setup_filter, filter_reads


# Set data path 
datapath = '/space/musselle/datasets/gazellesAndZebras'
files = 'testdata_5percent.bgzf'

# Process and Correct MID tag 



# Set filter
f = setup_filter({'phred': 15, 'propN': 0.05})
outdir = 'phredprop_testdata005'
filter_reads(infiles=files, datapath=datapath, outdir=outdir, filterfunc=f)




filtered_file = 'testdata_5percent-pass.bgzf'





# Get Stats








if __name__ == '__main__':
    pass