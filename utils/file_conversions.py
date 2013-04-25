'''
Created on 19 Mar 2013

@author: musselle
'''
import os
import sys
import time
import gzip
from subprocess import call
from collections import Counter
import glob

import numpy as np
from Bio import SeqIO, bgzf
import editdist as ed

from utils import Cycler, smartopen

def file2bgzf(infiles=None, filepattern=False, data_inpath='', SQLindex=True):
    ''' Convert given list of files from .gz or .fastq to .bgzf,
    And also produce an SQL index if needed. 
    
    infiles accepts str of file name of list of str for filenames. 
    If not specified will look at file type and glob the result to infiles. 
  
    '''
    if data_inpath:
        os.chdir(data_inpath)
  
    # Handle multiple types of input for infiles
    assert infiles is not None, 'No files listed or file pattern specified.'         
    if filepattern:
        # Fetch files by file types using glob
        import glob 
        infiles = glob.glob(infiles)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]
  
    start_time = time.time() 
    
    for filename in infiles: 
        toc = time.time()
        
        f = smartopen(filename)
        
        # Checks for type of input
        if filename.endswith('.gz'):
            # Drop .gz and append .bgzf
            bgzfFileName = '.'.join(filename.split('.')[:-1]) + '.bgzf'  
        elif filename.endswith('.fastq'):
            # Append .bgzf
            bgzfFileName = filename + '.bgzf'

        print "Producing BGZF output from {0}...".format(filename)
        w = bgzf.BgzfWriter(bgzfFileName, 'wb')
        while True:
            data = f.read(65536)
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
 
def makeSQLindex(infiles=None, data_inpath='', mode='grouped', outname=None):
    ''' Creates an SQL index out of either an uncompressed file or a compressed .bgzf file 
    
    if infiles is a string it is interpreted as a glob
    
    if infiles is list, goes through all file names in list.
    
     - mode  - grouped: all files are indexed to a single index file, specified by outname 
    
    '''
    
    starting_dir = os.getcwd()
    
    if data_inpath:
        os.chdir(data_inpath)
        
    if outname is None:
        outname = 'reads.idx'
  
    if type(infiles) is str:
        # Fetch files by file types using glob
        import glob 
        infiles = glob.glob(infiles)
    elif type(infiles) is not list and type(infiles) is not tuple:
        raise Exception("Invalid input files specified.")

    assert infiles, 'No files found, or no files passed.'

    # Handle multiple types of input for infiles
    if mode == 'grouped':
        idx_filename = outname
        tak = time.time()
        print 'Writing {0} files to SQL index ...'.format(len(infiles))
        SeqIO.index_db(idx_filename, infiles , 'fastq')
        idx_t = time.time() - tak
        print 'Finished Indexing to {0}\n after {1}\n'.format(idx_filename, time.strftime('%H:%M:%S', time.gmtime(idx_t)))

    elif mode == 'separate':
    
        for filename in infiles: 
            tak = time.time()
            print 'Writing SQL index file for {0} ...'.format(filename)
            idx_filename = filename.split('.')[0] + '.idx'
            SeqIO.index_db(idx_filename, filename , 'fastq')
            print '{0} written successfully'.format(idx_filename)
            idx_t = time.time() - tak
            print 'Finished Indexing after {1}\n'.format(time.strftime('%H:%M:%S', time.gmtime(idx_t)))

    if os.getcwd() != starting_dir:
        os.chdir(starting_dir) 


def file2fasta(filename):
    ''' Convert fastq file to fasta file '''
    
    handle = smartopen(filename)
    out_filename = filename.split('.')[0] + '.fasta'
    
    count = SeqIO.convert(handle, 'fastq', out_filename, 'fasta')
    
    print 'Converted {0} records to file\n{1}'.format(count, out_filename)

    






if __name__ == '__main__':
    pass