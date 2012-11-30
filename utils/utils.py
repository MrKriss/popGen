'''
Created on 20 Nov 2012

@author: musselle
'''
import os
import sys
import cPickle as pkl
from subprocess import PIPE, Popen
import glob

import gzip
import numpy as np
from Bio import SeqIO


class Cycler(object):
    ''' Object to hold generators that yield Sequence record objects or
    generators for all sequences records from the given file list

    Takes care of file input checks. Each next call returns the next file
    containing the next set of Sequence record objects.

    INPUTS
    inFiles - Single file as string or list of files to process as strings with
                full extensions.
    fileType - if no inFiles specified, runs glob on this string pattern in the
                path specified by dataPath e.g. glob.glob('*.fastq') for all
                files ending in .fastq
    dataPath - directory data files are stored in, will change to this at start
               of script.

    METHODS
    __init_files_gen - initiate files generator
    __init_rec_gen - initiate record generator

    ATTRIBUTES
    .seqfilegen - a generator that returns a generator per file to iterate over
                its sequence records.
    .recgen - a generator that iterates over all records serially for all files
                in list.
    .curfilename - Current file name in list of given files (inFiles) being
                iterated over.
    .curfilenum - Current file number in list of given files (inFiles).
    .maxNoSeq - Number of sequence records to return from each file. Default is
                all. Useful for testing.

    '''

    def __init__(self, infiles=None, filetype='', datapath='', maxnumseq=None):
        ''' Constructor '''

        if datapath:
            os.chdir(datapath)
        self.datapath = datapath
         
        # Handle multiple types of input for infiles
        if not infiles:
            # Fetch files by file types
            assert filetype, 'No files listed and No file type specified.'
            infiles = glob.glob(filetype)
        elif type(infiles) == str:
            # Convert to list
            infiles = [infiles]

        print '\nGenerator initiated to process the following files...'
        for f in infiles:
            print f

        self.numfiles = len(infiles)
        self.maxnumseq = maxnumseq

        '''May need to process files individually or all at once, so two types
        of generators needed'''
        # Generator that yields a sequence record iterator (generator)
        # for the next file in the list
        self.seqfilegen = self.__init_files_gen(infiles)

        # Generator that yields the next an individual record by running though
        # all files given in the list.
        self.recgen = self.__init_rec_gen()        

    def __init_files_gen(self, infiles):
        ''' Constructs the next file generator object '''
        # Generic file handling code
        for filenum, filename in enumerate(infiles):

            if self.datapath and os.getcwd() != self.datapath:
                os.chdir(self.datapath)

            self.curfilenum = filenum
            self.curfilename = filename

            # Check file extensions
            if filename.endswith('.idx'):
                # Raise warning if first call to generator
                print '''Warning: Processing .idx files takes ~3 times longer
                than using .bgzf or .fastq files.'''
                print 'Processing .idx files currently unsuported.'
                print '''Accepted formats for cycling through all records = .gz
                .bgzf and .fastq'''
                sys.exit()

            elif (filename.endswith('.gz') or filename.endswith('.bgzf') or
                                            filename.endswith('.fastq')):
                try:
                    yield SeqIO.parse(smartopen(filename), format='fastq')
                except IOError as e:
                    print e
                    print 'Invalid file name, or {0} may not exist.'.format(filename)
            else:
                print 'File extension for {0} currently unsupported.'.format(filename) 
                print 'Accepted formats for cycling through all records = .gz .bgzf and .fastq'

    def __init_rec_gen(self):
        ''' Return next Sequence Record '''

        if self.maxnumseq:
            for rec_file in self.seqfilegen:
                count = 0
                for record in rec_file:
                    count += 1
                    if count <= self.maxnumseq: 
                        yield record
                    else:
                        break
        else:
            for rec_file in self.seqfilegen:
                for record in rec_file:
                    yield record
        

def pklsave(obj, fileName):
    ''' Pickle the given object '''
    with open(fileName + '.pkl', 'wb') as f:
        pkl.dump(obj, f)
        
def smartopen(filename,*args,**kwargs):
    '''opens with open unless file ends in .gz or .bgzf, then use gzip.open

    in theory should transparently allow reading of files regardless of compression
    
    gzip can read in both .gz and .bgzf files.
    '''
    if filename.endswith('.gz') or filename.endswith('.bgzf'):
        return gzip.open(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)

def find_numrec(filename):
    ''' Find number of records in a fastq file. 
    
    Makes assumption that all records come from same machine and have a fastq format
    file that shares the filename given. 
    
    Uses subprocess to call 'grep -c 'filename'
    '''
    
    if filename.endswith('.idx') or filename.endswith('.gz') \
                            or filename.endswith('.bgzf') or filename.endswith('.fastq'):
        # Check fastq file exists.
        if filename.endswith('.fastq'):
            fastq_filename = filename
        else:
            fastq_filename = filename.split('.')[0] + '.fastq' 
        try:
            with open(fastq_filename, 'rb') as f:
                # Load in machine name
                sequencerID = f.readline().strip().split(':')[0]
        except IOError as e:
                print e
                print 'Invalid file name, or {0} may not exist.'.format(fastq_filename)
        filename = fastq_filename
    else:
        print 'File extension not recognised.'
        sys.exit()
    
#    cmd = 'cat {0} | parallel --pipe --block 2M grep -c "{1}"'.format(filename, sequencerID)
    cmd = 'grep -c "{1}" {0}'.format(filename, sequencerID)
    
    process = Popen(cmd, stdout=PIPE, shell=True)

    return int(process.communicate()[0].strip())

if __name__ == '__main__':
    
    os.chdir('/space/musselle/datasets/gazellesAndZebras/lane6')
    out = find_numrec('lane6_NoIndex_L006_R1_001.fastq')