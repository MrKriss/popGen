'''
Created on Oct 15, 2012

@author: chris
'''
import os 
from Bio.bgzf import BgzfWriter
import gzip
from Bio import SeqIO
import time


def gz2bgzf(inFiles = None, fileType = '', dataPath = '', SQLindex = True):
    ''' Convert the list of files from .gz to .bgzf,
    And also produce an SQL index if needed. 
    
    inFiles accepts str of file name of list of str for filenames. 
    If not specified will look at file type and glob the result to inFiles. 
  
    '''
    if dataPath:
        os.chdir(dataPath)
  
    # Handle multiple types of input
    if not inFiles:
        # Fetch files by file types
        assert fileType, 'No files listed and No file type specified.'
        import glob
        inFiles = glob.glob(fileType)
    elif type(inFiles) == str:
        # Convert to list
        inFiles = [inFiles]
  
    start_time = time.time() 
    
    for fileName in inFiles: 
        toc = time.time()
        
        # Checks for type of input
        if fileName.split('.')[-1] == '.gz':
            # Drop .gz and append .bgzf
            f2read = gzip.open(fileName)
            bgzfFileName = '.'.join(fileName.split('.')[:-1]) + '.bgzf'  
        elif fileName.split('.')[-1] == '.fastq':
            f2read = open(fileName)
            # Append .bgzf
            bgzfFileName = fileName + '.bgzf'

        print "Producing BGZF output from {0}...".format(fileName)
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
  
def makeSQLindex(inFiles = None, fileType = '', dataPath = ''):
    ''' Creates an SQL index out of either an uncompressed file or a compressed .bfzf file 
    
    if inFiles is list, goes through all file names in list
    
    '''
    if dataPath:
        os.chdir(dataPath)
  
    # Handle multiple types of input
    if not inFiles:
        # Fetch files by file types
        assert fileType, 'No files listed and No file type specified.'
        import glob
        inFiles = glob.glob(fileType)
    elif type(inFiles) == str:
        # Convert to list
        inFiles = [inFiles]

    for fileName in inFiles: 
        tak = time.time()
        print 'Writing SQL index file for {0} ...'.format(fileName)
        idxFileName = fileName.split('.')[0] + '.idx'
        SeqIO.index_db(idxFileName, fileName , 'fastq')
        print '{0} written successfully'.format(idxFileName)
        idx_t = time.time() - tak
        print 'Finished Indexing to {0}\n after {1}\n'.format(idxFileName, time.strftime('%H:%M:%S', time.gmtime(idx_t)))
    
if __name__ == '__main__':
    
    dataPath = '/space/musselle/datasets/gazellesAndZebras/'
    
    gz2bgzf(None, '*.gz', dataPath = dataPath + 'lane6/')
    gz2bgzf(None, '*.gz', dataPath = dataPath + 'lane8/')
    
    
    
    
    
    
    
    
    
    