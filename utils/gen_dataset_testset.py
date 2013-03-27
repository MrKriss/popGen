'''
Created on 30 Nov 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import time 


from Bio import SeqIO, bgzf
from utils import Cycler


start_dir = os.getcwd()

#lane6path = '/space/musselle/data/RAD-seq/gazelles-zebras/raw-data'
#lane8path = '/space/musselle/data/RAD-seq/gazelles-zebras/raw-data'
path = '/space/musselle/data/RAD-seq/gazelles-zebras/raw-data'

maxnumseq = 1000

#RecCyclerL6 = Cycler(filepattern='*[0-9].fastq.bgzf', data_inpath = lane6path, maxnumseq=maxnumseq)
#RecCyclerL8 = Cycler(filepattern='*[0-9].fastq.bgzf', data_inpath = lane8path, maxnumseq=maxnumseq)

# Just from lane 8
RecCycler = Cycler(filepattern='lane8*[0-9].fastq.bgzf', data_inpath = path, maxnumseq=maxnumseq)

output_filename = 'testset_1000.fastq.bgzf'

print '\nGenerating a dataset of the first {0} reads from each file.'.format(maxnumseq)

# If file already exists, overwrite it.
os.chdir('/space/musselle/data/RAD-seq/gazelles-zebras/testset')
if os.path.isfile(output_filename):
    f = open(output_filename, 'w')
    f.close()

output_filehdl = bgzf.BgzfWriter(output_filename, mode='a')

total_numwritten = 0

toc = time.time()

numwritten = SeqIO.write(RecCycler.recgen, output_filehdl , 'fastq')
print '{0} records written'.format(numwritten)
total_numwritten += numwritten

#numwritten = SeqIO.write(RecCycler.recgen , output_filehdl , 'fastq')
#print '{0} records written'.format(numwritten)
#total_numwritten += numwritten

output_filehdl.close()
print 'Total of {0} Sequence reads written to file {1}'.format(total_numwritten, output_filename) 

total_t = time.time() - toc    
print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                        time.gmtime(total_t)))
os.chdir(start_dir)
