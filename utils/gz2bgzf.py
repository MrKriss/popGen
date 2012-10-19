'''
Created on Oct 15, 2012

@author: chris
'''
import os 
from Bio.bgzf import BgzfWriter
import gzip
from Bio import SeqIO
import time


def gz2bgzf(inFiles, fileType = '', dataPath = '' SQLindex = True):
    ''' Convert the list of files from .gz to .bgzf, A
    And also produce an SQL index if needed. 
    
    TODO: 
    
    inFiles can either be a string of the file name, a list of file names
    or a regular expression of the files to 
    
    '''
    
  
    if dataPath:
        os.chdir(dataPath)
    
    if type(inFiles) == str:
        inFiles = [inFiles]
  
    start_time = time.time() 
    
    for fileName in inFiles: 
        toc = time.time()
        
        # Checks for type of input
        if fileName.split('.')[-1] == '.gz':
            f2read = gzip.open(fileName)
            bgzfFileName = '.'.join(fileName.split('.')[:-1]) + '.bgzf'
            
        if fileName.split('.')[-1] == '.fastq':
            f2read = open(fileName)
            bgzfFileName = fileName + '.bgzf'

        print "Producing BGZF output from {0}...".format(fileName)
        # Drop .gz and append .bgzf
        
        w = BgzfWriter(bgzfFileName, 'wb')
        while True:
            data = f2read.read(65536)
            w.write(data)
            if not data:
                break
        w.close()
        print '{0} written successfully'.format(bgzfFileName)

        conv_t = time.time() - toc 
        print 'Finished converting {0}\n after {1}\n'.format(fileName, time.strftime('%H:%M:%S', time.gmtime(conv_t)))
        
        if SQLindex == True:
            makeSQLindex(bgzfFileName)
      
    total_t = time.time() - start_time
    print 'Finished all processing {0} files in {1}'.format(len(inFiles), time.strftime('%H:%M:%S', time.gmtime(total_t)))
  
def makeSQLindex(inFiles):
    ''' Creates an SQL index out of either an uncompressed file or a compressed .bfzf file 
    
    if inFiles is list, goes through all file names in list
    
    '''
    if type(inFiles) == str:
        inFiles = [inFiles]

    for fileName in inFiles: 
        tak = time.time()
        print 'Writing SQL index file for {0} ...'.format(fileName)
        idxFileName = '.'.join(fileName.split('.')[:-2]) + '.idx'
        SeqIO.index_db(idxFileName, fileName , 'fastq')
        print '{0} written successfully'.format(idxFileName)
        idx_t = time.time() - tak
        print 'Finished Indexing to {0}\n after {1}\n'.format(idxFileName, time.strftime('%H:%M:%S', time.gmtime(idx_t)))
    
if __name__ == '__main__':
    
    dataPath = '/space/musselle/datasets/gazellesAndZebras/'
    os.chdir(dataPath)
    
    # Convert all gz files to bgzf files
    filesLane8_1 = ["lane8_NoIndex_L008_R1_00%i.fastq.bgzf" % (i+1) for i in range(9)]
    filesLane8_2 = ["lane8_NoIndex_L008_R1_0%i.fastq.bgzf" % (i+10) for i in range(1)]
    filesLane6_1 = ["lane6_NoIndex_L006_R1_00%i.fastq.bgzf" % (i+1) for i in range(9)]
    filesLane6_2 = ["lane6_NoIndex_L006_R1_0%i.fastq.bgzf" % (i+10) for i in range(2)]
    
#    files = filesLane8_1 + filesLane8_2 + filesLane6_1 + filesLane6_2 
#    files = filesLane8_1 + filesLane8_2 
    files = filesLane6_1 + filesLane6_2 
    
#    files = ['lane8_NoIndex_L008_R1_001.fastq.gz']
    
#    gz2bgzf(files)
    dataPath = '/space/musselle/datasets/gazellesAndZebras/'
    os.chdir(dataPath +'lane6/')
    makeSQLindex(files)
    
    
    
    
    
    
    
    
    