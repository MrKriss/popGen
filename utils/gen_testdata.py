'''
Created on 30 Nov 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import time 
import socket

from Bio import SeqIO, bgzf
from utils import Cycler

start_dir = os.getcwd()

# Work out where data is stored
if socket.gethostname() == 'yildun':
    prefix = '/space/musselle/datasets'
elif socket.gethostname() == 'luca':
    prefix = '/home/musselle/san/data'


path = os.path.join(prefix, 'sticklebacks/')

maxnumseq = 10000 # 5% of 16,000,000

RecCycler = Cycler(filepattern='*[0].fastq.bgzf', data_inpath = path, maxnumseq=maxnumseq)

output_filename = 'sb_testdata.bgzf'

print '\nGenerating a dataset of the first {0} reads from each file.'.format(maxnumseq)

# If file already exists, overwrite it.
os.chdir(path)
if os.path.isfile(output_filename):
    f = open(output_filename, 'w')
    f.close()

output_filehdl = bgzf.BgzfWriter(output_filename, mode='a')

toc = time.time()

numwritten = SeqIO.write(RecCycler.recgen, output_filehdl , 'fastq')
print '{0} records written'.format(numwritten)

output_filehdl.flush()
output_filehdl.close()

print 'Total of {0} Sequence reads written to file {1}'.format(numwritten, output_filename) 

total_t = time.time() - toc    
print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                        time.gmtime(total_t)))
os.chdir(start_dir)
