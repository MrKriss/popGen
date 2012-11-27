'''
Created on 20 Nov 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import glob
import cPickle as pkl
import gzip
from subprocess import PIPE, Popen


from Bio import SeqIO

class Cycler(object):
    ''' Object to hold generators that yield Sequence record objects or 
    generators for all sequences records from the given file list. 
    
    Takes care of file input checks. Each next call returns the next file
    containing the next set of Sequence record objects. 
    
    INPUTS
    inFiles - Single file as string or list of files to process as strings with full extensions. 
    fileType - if no inFiles specified, runs glob on this string pattern in the path specified by 
                dataPath e.g. glob.glob('*.fastq') for all files ending in .fastq
    dataPath - directory data files are stored in, will change to this at start of script.
    
    
    METHODS
    __init_files_gen - initiate files generator
    __init_rec_gen - initiate record generator
     
    ATTRIBUTES
    .seqFileGen - a generator that returns a generator per file to iterate over its sequence records. 
    .recGen - a generator that iterates over all records serially for all files in list. 
    .curFileName - Current file name in list of given files (inFiles) being iterated over. 
    .curFileNum - Current file number in list of given files (inFiles).
    .maxNoSeq - Number of sequence records to return from each file. Default is all. Useful for testing. 
    
    '''     

    def __init__(self, inFiles = None, fileType = '', dataPath = '', maxNumSeq = None):
        ''' Constructor '''
        
        if dataPath:
            os.chdir(dataPath)
      
        # Handle multiple types of input for inFiles
        if not inFiles:
            # Fetch files by file types
            assert fileType, 'No files listed and No file type specified.'
            import glob
            inFiles = glob.glob(fileType)
        elif type(inFiles) == str:
            # Convert to list
            inFiles = [inFiles]
    
        print '\nProcessing the following files...'
        for f in inFiles:
            print f
            
        self.numFiles = len(inFiles)
        self.maxNumSeq = maxNumSeq
        
        '''May need to process files individually or all at once, so two types of
        generators needed''' 
        # Generator that yields a sequence record iterator (generator)
        # for the next file in the list 
        self.seqFileGen = self.__init_files_gen(inFiles)
        
        # Generator that yields the next an individual record by running though all files
        # Given in the list  
        self.recGen = self.__init_rec_gen()
        
    def __init_files_gen(self, inFiles):
        ''' Constructs the next file generator object '''
        # Generic file handling code
        for fileNum, fileName in enumerate(inFiles):   
            
            self.curFileNum = fileNum     
            self.curFileName = fileName
                 
            # Check file extensions
            if fileName.endswith('.idx'):
                # Raise warning if first call to generator
                print 'Warning: Processing .idx files takes ~3 times longer than using .bgzf or .fastq files.'
                print 'Processing .idx files currently unsuported.'
                print 'Accepted formats for cycling through all records = .gz .bgzf and .fastq'
                sys.exit()
               
            elif fileName.endswith('.gz') or fileName.endswith('.bgzf') or fileName.endswith('.fastq'): 
                try:
                    yield SeqIO.parse(smartopen(fileName), format = 'fastq')
                except IOError as e:
                    print e
                    print 'Invalid file name, or {0} may not exist.'.format(fileName)
            else:
                print 'File extension for {0} currently unsupported.'.format(fileName) 
                print 'Accepted formats for cycling through all records = .gz .bgzf and .fastq'
    
    def __init_rec_gen(self):
        ''' Return next Sequence Record '''
        
        if self.maxNumSeq:
            for rec_file in self.seqFileGen:
                count = 0
                for record in rec_file:
                    count += 1
                    if count <= self.maxNumSeq:  
                        yield record
                    else:
                        break
        else:
            for rec_file in self.seqFileGen:
                for record in rec_file:
                    yield record
        

def pklSave(obj, fileName):
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

def findNumRec(fileName):
    ''' Find number of records in a fastq file. 
    
    Makes assumption that all records come from same machine and are in fastq format
    
    Uses subprocess properly to call grep in parallel for a list of files.
    '''
    
    if fileName.endswith('.idx') or fileName.endswith('.gz') \
                            or fileName.endswith('.bgzf') or fileName.endswith('.fastq'):
        # Check fastq file exists.
        if fileName.endswith('.fastq'):
            fastqFileName = fileName
        else:
            fastqFileName = fileName.split('.')[0] + '.fastq' 
        try:
            with open(fastqFileName, 'rb') as f:
                # Load in machine name
                sequencerID = f.readline().strip().split(':')[0]
        except IOError as e:
                print e
                print 'Invalid file name, or {0} may not exist.'.format(fastqFileName)
        fileName = fastqFileName
    else:
        print 'File extension not recognised.'
        sys.exit()
    
#    cmd = 'cat {0} | parallel --pipe --block 2M grep -c "{1}"'.format(fileName, sequencerID)
    cmd = 'grep -c "{1}" {0}'.format(fileName, sequencerID)
    
    process = Popen(cmd, stdout=PIPE, shell=True)

    return int(process.communicate()[0].strip())

if __name__ == '__main__':
    
    os.chdir('/space/musselle/datasets/gazellesAndZebras/lane6')
    out = findNumRec('lane6_NoIndex_L006_R1_001.fastq')