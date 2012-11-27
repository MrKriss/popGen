'''
Created on 19 Nov 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

from basicAnalysis import Cycler
from Bio import SeqIO, bgzf
import gzip
import time

def filterReads(inFiles = None, fileType = '', dataPath = ''):
    ''' Filter reads based on criteria 
    
    Default is to use Machine Specific read filter 
    
    Only specific to Casava 1.8 Illumina output format at current
    
    '''   
    
    # Have to create two separate generators as copying a generator is not possible
    RecCycler1 = Cycler(inFiles = inFiles, fileType = fileType, dataPath = dataPath)
    RecCycler2 = Cycler(inFiles = inFiles, fileType = fileType, dataPath = dataPath)
    
    toc = time.time()
    cum_t = 0
    
    for recordGen1 in RecCycler1.seqFileGen:
        # Both gens are identical at this point so iterating over one is the same length as the other
        
        # Load second record generator
        recordGen2 = RecCycler2.seqFileGen.next()

        print 'Processing {0}'.format(RecCycler1.curFileName) 
        
        # Construct file names
        name = RecCycler1.curFileName.split('.')  
        
        passFileName = [name[0] + '-pass'] + name[1:] 
        passFileName = '.'.join(passFileName)
        failFileName = [name[0] + '-fail']  + name[1:]
        failFileName = '.'.join(failFileName)
        name = '.'.join(name)
        
        # Setup Generators              
        passGen = (rec for rec in recordGen1 if rec.description.split()[1].split(':')[1] == 'N')
        failGen = (rec for rec in recordGen2 if rec.description.split()[1].split(':')[1] == 'Y')
        
        if name.endswith('.bgzf'):
            passFileHdl = bgzf.BgzfWriter(passFileName)
            failFileHdl = bgzf.BgzfWriter(failFileName)
        elif name.endswith('.fastq'):
            passFileHdl = open(passFileName, 'wb')
            failFileHdl = open(failFileName, 'wb')
        elif name.endswith('.gz'):
            passFileHdl = gzip.open(passFileName, 'wb')
            failFileHdl = gzip.open(failFileName, 'wb')
        else:
            print 'Input file format not supported'
            sys.exit()
        
        print 'Writing passes to file {0} ....'.format(passFileName)
        numWritten = SeqIO.write( passGen , passFileHdl , 'fastq')
        passFileHdl.close()
        print '{0} records written'.format(numWritten)
        
        print 'Writing fails to file  {0} ....'.format(failFileName)
        numWritten = SeqIO.write( failGen , failFileHdl , 'fastq')
        failFileHdl.close()
        print '{0} records written'.format(numWritten)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished file {0} after {1}'.format(RecCycler1.curFileNum, time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
        
    total_t = time.time() - toc    
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', time.gmtime(total_t)))

if __name__ == '__main__':
    
    dataLoc = '/space/musselle/datasets/gazellesAndZebras/lane6'
    filterReads(None, fileType='*[0-9].fastq.bgzf', dataPath = dataLoc)
    
    
    
    