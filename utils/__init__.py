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

# TODO: Write function to sumarise fails.log files


class Cycler(object):
    ''' Object to hold generators that yield Sequence record objects or
    generators for all sequences records from the given file list

    Takes care of file input checks. Each next call returns the next file
    containing the next set of Sequence record objects.
    
    NOTE: list of filenames is returned in random order. could sort, but at 
    current no real need, as it only needs to cycle through all files.  
    

    INPUTS
    infiles - Single file as string or list of files to process as strings with
                full extensions.
    filepattern - if no infiles specified, runs glob on this string pattern in the
                path specified by data_inpath e.g. glob.glob('*.fastq') for all
                files ending in .fastq
    data_inpath - directory data files are stored in, will change to this at start
               of script.

    METHODS
    __init_files_gen - initiate files generator
    __init_rec_gen - initiate record generator

    ATTRIBUTES
    .seqfilegen - a generator that returns a generator per file to iterate over
                its sequence records.
    .recgen - a generator that iterates over all records serially for all files
                in list.
    .curfilename - Current file name in list of given files (infiles) being
                iterated over.
    .curfilenum - Current file number in list of given files (infiles).
    .maxNoSeq - Number of sequence records to return from each file. Default is
                all. Useful for testing.

    '''

    def __init__(self, infiles=None, filepattern='', data_inpath='', maxnumseq=None):
        ''' Constructor '''

#        if data_inpath:
#            if data_inpath not in sys.path:
#                sys.path.insert(0, data_inpath)
        
        # Remember data_inpath 
        self.data_inpath = data_inpath
         
        # Handle multiple types of input for infiles
        if not infiles:
            # Fetch files by file types
            assert filepattern, 'No files listed and No file type specified.'
            infiles = glob.glob(os.path.join(data_inpath,filepattern))
        elif type(infiles) == str:
            # Convert to list
            infiles = [infiles]

        print '\nGenerator initiated to process the following files...'
        for f in infiles:
            print f

        self.numfiles = len(infiles)
        self.infiles = infiles
        
        # Max number of record to run the generator for (default = all)
        self.maxnumseq = maxnumseq

        ''' May need to process files individually or all at once, so two types
        of generators needed '''
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

            self.curfilenum = filenum
            self.curfilename = filename

            next_data_file_loc = os.path.join(self.data_inpath, filename)

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
                    yield SeqIO.parse(smartopen(next_data_file_loc), format='fastq')
                except IOError as e:
                    print e
                    raise Exception(('Invalid file name {0},'
                     ' or path {1} may not exist.').format(filename, next_data_file_loc))
            else:
                print 'File extension for {0} currently unsupported.'.format(filename) 
                print 'Accepted formats for cycling through all records = .gz .bgzf and .fastq'
                raise Exception

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

import string
def search_file(filename, search_path, pathsep=os.pathsep):
    """ Given a search path, find file with requested name """
    for path in string.split(search_path, pathsep):
        candidate = os.path.join(path, filename)
        if os.path.exists(candidate): return os.path.abspath(candidate)
    return None

#if _ _name_ _ == '_ _ _main_ _':
#    search_path = '/bin' + os.pathsep + '/usr/bin'  # ; on Windows, : on Unix
#    find_file = search_file('ls',search_path)
#    if find_file:
#        print "File found at %s" % find_file
#    else:
#        print "File not found"

def pklsave(obj, filename):
    ''' Pickle the given object '''
    if not filename.endswith('pkl'):
        filename = filename + '.pkl'
    with open(filename, 'wb') as f:
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

def make_MIDdict(infiles=None, filepattern=False, data_inpath=''):
    ''' Function to load in MIDs from a list of files into a dictionary '''

    # Handle multiple types of input for infiles
    assert infiles is not None, 'No files listed or file pattern specified.'         
    if filepattern:
        # Fetch files by file types using glob
        import glob 
        st_dir = os.getcwd()
        os.chdir(data_inpath)
        infiles = glob.glob(infiles)
        os.chdir(st_dir)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]

    # Run through files and store barcodes in a Dictionary object.
    # keys are starting tags (MID (6 BP) + cutsite (6BP))
    tags = {}
    
    for filename in infiles:
        with open(os.path.join(data_inpath,filename), 'rb') as f:
            for line in f:
                elem = line.split()
                tags[elem[0]] = elem[1] 
    return tags

if __name__ == '__main__':
    
    os.chdir('/space/musselle/datasets/gazellesAndZebras/lane6')
    out = find_numrec('lane6_NoIndex_L006_R1_001.fastq')