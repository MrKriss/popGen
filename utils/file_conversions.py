'''
Created on Oct 15, 2012

@author: chris
'''
import os 
from Bio.bgzf import BgzfWriter
import gzip
from Bio import SeqIO
import time


def gz2bgzf(infiles = None, filetype = '', datapath = '', SQLindex = True):
    ''' Convert the list of files from .gz to .bgzf,
    And also produce an SQL index if needed. 
    
    infiles accepts str of file name of list of str for filenames. 
    If not specified will look at file type and glob the result to infiles. 
  
    '''
    if datapath:
        os.chdir(datapath)
  
    # Handle multiple types of input
    if not infiles:
        # Fetch files by file types
        assert filetype, 'No files listed and No file type specified.'
        import glob
        infiles = glob.glob(filetype)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]
  
    start_time = time.time() 
    
    for filename in infiles: 
        toc = time.time()
        
        # Checks for type of input
        if filename.split('.')[-1] == '.gz':
            # Drop .gz and append .bgzf
            f2read = gzip.open(filename)
            bgzfFileName = '.'.join(filename.split('.')[:-1]) + '.bgzf'  
        elif filename.split('.')[-1] == '.fastq':
            f2read = open(filename)
            # Append .bgzf
            bgzfFileName = filename + '.bgzf'

        print "Producing BGZF output from {0}...".format(filename)
        w = BgzfWriter(bgzfFileName, 'wb')
        while True:
            data = f2read.read(65536)
            w.write(data)
            if not data:
                break
        w.close()
        print '{0} written successfully'.format(bgzfFileName)

        conv_t = time.time() - toc 
        print 'Finished converting {0}\n after {1}\n'.format(filename, time.strftime('%H:%M:%S', time.gmtime(conv_t)))
        
        if SQLindex == True:
            makeSQLindex(bgzfFileName)
      
    total_t = time.time() - start_time
    print 'Finished all processing {0} files in {1}'.format(len(infiles), time.strftime('%H:%M:%S', time.gmtime(total_t)))
  
def makeSQLindex(infiles = None, filetype = '', datapath = ''):
    ''' Creates an SQL index out of either an uncompressed file or a compressed .bgzf file 
    
    if infiles is list, goes through all file names in list
    
    '''
    if datapath:
        os.chdir(datapath)
  
    # Handle multiple types of input
    if not infiles:
        # Fetch files by file types
        assert filetype, 'No files listed and No file type specified.'
        import glob
        infiles = glob.glob(filetype)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]

    for filename in infiles: 
        tak = time.time()
        print 'Writing SQL index file for {0} ...'.format(filename)
        idx_filename = filename.split('.')[0] + '.idx'
        SeqIO.index_db(idx_filename, filename , 'fastq')
        print '{0} written successfully'.format(idx_filename)
        idx_t = time.time() - tak
        print 'Finished Indexing to {0}\n after {1}\n'.format(idx_filename, time.strftime('%H:%M:%S', time.gmtime(idx_t)))
    
if __name__ == '__main__':
    
    dataPath = '/space/musselle/datasets/gazellesAndZebras/'
    
    gz2bgzf(None, '*.gz', datapath = dataPath + 'lane6/')
    gz2bgzf(None, '*.gz', datapath = dataPath + 'lane8/')
    
    
    
    
    
    
    
    
    
    