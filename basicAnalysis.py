'''
Created on Oct 9, 2012

@author: chris
'''

from Bio import SeqIO
import numpy as np
import time
import os
import sys
import matplotlib.pyplot as plt
import pickle as pkl
import subprocess


class Cycler(object):
    ''' Object to hold generators that yield Sequence record objects or 
    individual sequences from the given file list. 
    
    Takes care of file input checks. Each next call returns the next file
    containing the next set of Sequence record objects. 
    
    INPUTS
    inFiles - Single file as string or list of files to process as strings with full extensions. 
    fileType - if no inFiles specified, runs glob on this string pattern in the path specified by 
                dataPath e.g. glob.glob('*.fastq') for all files ending in .fastq
    dataPath - directory data files are stored in, will change to this at start of script.
    
    
    METHODS
    nextSeq - returns the next sequence object
    gen_files - stores a generator that yields the record generator for each file. 
     
    ATTRIBUTES
    .fileGen - the generator constructed by gen_files.
    .fileGen.next() - returns a tuple of (SeqenceRecordObject, fileName, fileNo )
    
    returns tuple(object of all sequence records, fileName , fileNo).
    '''     

    def __init__(self, inFiles = None, fileType = '', dataPath = ''):
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
    
        print 'Processing the following files...'
        for f in inFiles:
            print f
            
        self.numfiles = len(inFiles)
        
        self.fileGen = self.init_files_gen(inFiles)
        self.seqGen = self.init_seq_gen()
        
    def init_files_gen(self, inFiles):
        ''' Constructs the next file generator object '''
        # Generic file handling code
        for fileNo, fileName in enumerate(inFiles):   
            
            self.fileNo = fileNo     
            self.fileName = fileName
                 
            # Check file extensions
            extension = fileName.split('.')[-1] 
            if extension == 'idx':
                # Raise warning if first call to generator
                print 'Warning: Processing .idx files takes ~3 times longer than using .bgzf or .fastq files.'
                print 'Processing .idx files currently unsuported.'
                print 'Accepted formats for cycling through all records = .gz .bgzf and .fastq'
                sys.exit()
               
            elif extension == 'gz' or extension == 'bgzf': 
                import gzip
                try:
                    handle = gzip.open(fileName, "r")
                    yield SeqIO.parse(handle, format = 'fastq')
                except IOError as e:
                    print e
                    print 'Invalid file name, or {0} may not exist.'.format(fileName)
    
            elif extension == 'fastq':
                try:
                    yield SeqIO.parse(fileName, format = 'fastq')
                except IOError as e:
                    print e
                    print 'Invalid file name, or {0} may not exist.'.format(fileName)
            else:
                print 'File extension {0} currently unsupported.'.format(extension) 
                print 'Accepted formats for cycling through all records = .gz .bgzf and .fastq'
    
    def init_seq_gen(self):
        ''' Return next Sequence Record '''
        for rec_file in self.fileGen:
            for record in rec_file:
                yield record
        
    
def findNumRec(filename):
    ''' Find number of records in a fastq file. 
    
    Makes assumption that all records come from same machine and are in fastq format
    
    TODO - use subprocess properly to call grep in parallel for a list of files.
    
    Possible optimisation
    
    '''
    subprocess.call()
    
    filename
    

def getPropN_meanPhred(inFiles = None, fileType = '', dataPath = ''):
    ''' Return the proportion of Ns and the mean Phred scores of the files given'''

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

    print 'Calculating Proportion of Ns and Mean Phred Score per read for the following files...'
    for f in inFiles:
        print f

    # Define vars and outputs
    numFiles = len(inFiles)
    outList = [0] * numFiles
    recordCounts = [0] * numFiles

    # Find the number of records in each file.
    toc = time.time()
    cum_t = 0
    for fileNo, fileName in enumerate(inFiles):
        # Check its a .idx file extention or that one exists with it.
        if fileName.split('.')[-1] == 'idx':
            idxFileName = fileName
        else:
            assert os.path.isfile('.'.join(fileName.split('.')[:-2]) + '.idx'),  'Index file not found for {0}'.format(fileName)
            idxFileName = '.'.join(fileName.split('.')[:-2]) + '.idx'
    
        # Load file 
        S = SeqIO.index_db(idxFileName, format='fastq')
        recordCounts[fileNo] = len(S)
        
        outList[fileNo] = {'propN'  :  np.zeros(recordCounts[fileNo]),
                        'meanPhred' : np.zeros(recordCounts[fileNo])  }
        
        # Cycle through all records 
        for num, seqRec in enumerate(S.itervalues()):
            # Proportion of Ns
            outList[fileNo]['propN'][num] = float(seqRec.seq.count('N')) / len(seqRec)
            # Mean Phred Score
            outList[fileNo]['meanPhred'][num] = np.array(seqRec.letter_annotations['phred_quality']).mean()
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(fileName, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    totalRecords = np.sum(recordCounts)
    totalStats = {'propN':  np.zeros(totalRecords)  ,
                  'meanPhred': np.zeros(totalRecords) }
    start_idx = 0
    end_idx = 0
    for file_idx, statObj in enumerate(outList):
        end_idx += recordCounts[file_idx]
        totalStats['propN'][start_idx:end_idx] = statObj['propN']
        totalStats['meanPhred'][start_idx:end_idx] = statObj['meanPhred']
        start_idx += recordCounts[file_idx]

    # Saving data
    print 'Saving Data......'
    parent_dir  = os.getcwd().split('/')[-1]
    np.save(parent_dir + '_propN', totalStats['propN'])
    np.save(parent_dir + '_meanPhred', totalStats['meanPhred'])
    print 'DONE! Data saved to {0}'.format(os.getcwd())

    total_t = time.time() - toc
    print 'Processed {0} files in {1}'.format(numFiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))

    return totalStats

            
def getGen4SeqRecFiles(inFiles = None, fileType = '', dataPath = ''):
    ''' Return Generator for all Sequence Records stored in a .bgzf file.
    
    Takes care of file input checks. Each next call returns the next file
    containing the next set of Sequence record objects. 
    
    First next() call returns number of files to process as int(). 
    
    Rest returns tuple(object of all sequence records, fileName , fileNo).
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

    print 'Processing the following files...'
    for f in inFiles:
        print f
        
    yield len(inFiles)
    
    # Generic file handling code
    for fileNo, fileName in enumerate(inFiles):
        
        # Check file extensions
        extension = fileName.split('.')[-1] 
        if extension == 'idx':
            # Raise warning if first call to generator
            print 'Warning: Processing .idx files takes ~3 times longer than using .bgzf or .fastq files.'
            print 'Processing .idx files currently unsuported. Accepted formats = .gz .bgzf and .fastq'
            sys.exit()
           
        elif extension == 'gz' or extension == 'bgzf': 
            import gzip
            try:
                handle = gzip.open(fileName, "r")
                yield SeqIO.parse(handle, format = 'fastq'), fileName, fileNo
            except IOError as e:
                print e
                print 'Invalid file name, or {0} may not exist.'.format(fileName)

        elif extension == 'fastq':
            try:
                yield SeqIO.parse(fileName, format = 'fastq'), fileName, fileNo
            except IOError as e:
                print e
                print 'Invalid file name, or {0} may not exist.'.format(fileName)
        else:
            print 'File extension {0} currently unsupported.'.format(extension) 
            print 'Accepted formats = .gz .bgzf and .fastq'
            
            
def getReadLengths2(inFiles = None, fileType = '', dataPath = ''):
    ''' Return histogram of read lengths '''        
    
    #Generator for Sequence Record files
    SeqRecFileGenerator = getGen4SeqRecFiles(inFiles = inFiles, fileType = fileType, dataPath = dataPath)

    print 'Calculating length per read ...'
    
    # Define vars and outputs
    numFiles = SeqRecFileGenerator() # First call returns total number of files
    outList = [0] * numFiles
    recordCounts = [0] * numFiles

    toc = time.time()
    cum_t = 0
    for seqRecs, fileName, fileNo in SeqRecFileGenerator:
        # count number of Records 
        idxFileName = '.'.join(fileName.split('.')[:-2]) + '.idx'
        S = SeqIO.index_db(idxFileName, format='fastq')
        recordCounts[fileNo] = len(S)
        del S
        
        outList[fileNo] = {'length'  :  np.zeros(recordCounts[fileNo]) }
        
        for num, seqR in enumerate(seqRecs):
            outList[fileNo]['length'][num] = len(seqR.seq)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(fileName, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    totalRecords = np.sum(recordCounts)
    totalStats = {'length':  np.zeros(totalRecords)}
    start_idx = 0
    end_idx = 0
    for file_idx, statObj in enumerate(outList):
        end_idx += recordCounts[file_idx]
        totalStats['length'][start_idx:end_idx] = statObj['length']
        start_idx += recordCounts[file_idx]

    # Saving data
    print 'Saving Data......'
    parent_dir  = os.getcwd().split('/')[-1]
    np.save(parent_dir + '_length', totalStats['length'])
    print 'DONE! Data saved to {0}'.format(os.getcwd())

    total_t = time.time() - toc
    print 'Processed {0} files in {1}'.format(numFiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))

    return totalStats

def getReadLengths3(inFiles = None, fileType = '', dataPath = ''):
    ''' Return histogram of read lengths. Uses RecordCyler object'''        
    
    #Generator for Sequence Record files
    SeqRecCycler = RecordCycler(inFiles = inFiles, fileType = fileType, dataPath = dataPath)

    print 'Calculating length per read ...'
    
    # Define vars and outputs
    outList = [0] * SeqRecCycler.numFiles
    recordCounts = [0] * SeqRecCycler.numFiles

    toc = time.time()
    cum_t = 0
    for seqRecs in SeqRecCycler.fileGen:
        fileName = SeqRecCycler.fileGen.fileName
        fileNo = SeqRecCycler.fileGen.fileNo
        # count number of Records, may need to alter this if no idx file
        idxFileName = '.'.join(fileName.split('.')[:-2]) + '.idx'
        S = SeqIO.index_db(idxFileName, format='fastq')
        recordCounts[fileNo] = len(S)
        del S
        
        outList[fileNo] = {'length'  :  np.zeros(recordCounts[fileNo]) }
        
        for num, seqR in enumerate(seqRecs):
            outList[fileNo]['length'][num] = len(seqR.seq)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(fileName, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    totalRecords = np.sum(recordCounts)
    totalStats = {'length':  np.zeros(totalRecords)}
    start_idx = 0
    end_idx = 0
    for file_idx, statObj in enumerate(outList):
        end_idx += recordCounts[file_idx]
        totalStats['length'][start_idx:end_idx] = statObj['length']
        start_idx += recordCounts[file_idx]

    # Saving data
    print 'Saving Data......'
    parent_dir  = os.getcwd().split('/')[-1]
    np.save(parent_dir + '_length', totalStats['length'])
    print 'DONE! Data saved to {0}'.format(os.getcwd())

    total_t = time.time() - toc
    print 'Processed {0} files in {1}'.format(SeqRecCycler.numFiles, 
                                              time.strftime('%H:%M:%S', time.gmtime(total_t)))

    return totalStats
         
def getReadLengths1(inFiles = None, fileType = '', dataPath = ''):
    ''' old read Length Function'''

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

    print 'Calculating length per read for the following files...'
    for f in inFiles:
        print f

    # Define vars and outputs
    numFiles = len(inFiles)
    outList = [0] * numFiles
    recordCounts = [0] * numFiles

    # Find the number of records in each file.
    toc = time.time()
    cum_t = 0
    for fileNo, fileName in enumerate(inFiles):
        # Check its a .idx file extention or that one exists with it.
        if fileName.split('.')[-1] == 'idx':
            idxFileName = fileName
        else:
            assert os.path.isfile('.'.join(fileName.split('.')[:-2]) + '.idx'),  'Index file not found for {0}'.format(fileName)
            idxFileName = '.'.join(fileName.split('.')[:-2]) + '.idx'
    
        # Load file 
        S = SeqIO.index_db(idxFileName, format='fastq')
        recordCounts[fileNo] = len(S)
        
        outList[fileNo] = {'length'  :  np.zeros(recordCounts[fileNo]) }
        
        # Cycle through all records 
        for num, seqRec in enumerate(S.itervalues()):
            # Find Length
            outList[fileNo]['length'][num] = len(seqRec)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(fileName, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    totalRecords = np.sum(recordCounts)
    totalStats = {'length':  np.zeros(totalRecords)}
    start_idx = 0
    end_idx = 0
    for file_idx, statObj in enumerate(outList):
        end_idx += recordCounts[file_idx]
        totalStats['length'][start_idx:end_idx] = statObj['length']
        start_idx += recordCounts[file_idx]

    # Saving data
    print 'Saving Data......'
    parent_dir  = os.getcwd().split('/')[-1]
    np.save(parent_dir + '_length', totalStats['length'])
    print 'DONE! Data saved to {0}'.format(os.getcwd())

    total_t = time.time() - toc
    print 'Processed {0} files in {1}'.format(numFiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))

    return totalStats


def boxPlotPhredPerBase(inFiles = None, fileType = '', dataPath = ''):
    ''' Find the median, upper and lower quartile for the Phred score per base 
    
    Returns the stats and the counter dictionary. 
    
    Counter dictionary may become standard way to store mass Phred/seq bases data.  '''
    
    from collections import Counter
    
    #Generator for Sequence Records
    SeqRecFileGenerator = getGen4SeqRecFiles(inFiles = inFiles, fileType = fileType, dataPath = dataPath)
    numFiles = SeqRecFileGenerator.next()

    toc = time.time()
    cum_t = 0
    
    counterList = [0] * 101
    for i in range(len(counterList)):
        counterList[i] = Counter()
    
    for seqRecs, fileName, fileNo in SeqRecFileGenerator:
        for rec in seqRecs:
            for baseNo, phred in enumerate(rec.letter_annotations['phred_quality']):
                counterList[baseNo][phred] += 1
                
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished file {0} of {1} after {2}'.format(fileNo, numFiles, time.strftime('%H:%M:%S', time.gmtime(loop_t))) 

    # Calculate min, max Q1, Q2, Median and Average
    stats = getStatsFT(counterList)

    total_t = time.time() - toc
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', time.gmtime(total_t)))
            
    return stats, counterList     


def getStatsFT(counterList):
    ''' Calculate min, max, Q1, Q2, Median and Average for a list of frequency table of data 
    
    Returns a structured array the length of the input list '''
    
    # Input check    
    if type(counterList) != list:
        counterList = [counterList] 

    # Define data type    
    dt = np.dtype([ ('min', 'u2'), 
                ('Q1', 'f2'), 
                ('median','f2'), 
                ('ave', 'f2'), 
                ('Q2', 'f2'), 
                ('max', 'u2') ])
    
    numCounters = len(counterList)
    stats = np.zeros(numCounters, dtype = dt)
    
    # Calculate min, max Q1, Q2, Median and Average
    for i, cnt in enumerate(counterList):
        keys = cnt.keys()
        vals = percentileFT(cnt, [0.25, 0.5, 0.75])
        stats[i]['Q1'] = vals[0]
        stats[i]['median'] = vals[1]
        stats[i]['Q2'] = vals[2]
        stats[i]['min'] = min(keys)
        stats[i]['max'] = max(keys)
        ave = float(0)
        for k in keys:
            ave = ave + (k * cnt[k])
        ave = ave / sum(cnt.values()) 
        stats[i]['ave'] = ave

    return stats

def percentileFT(hashCount, percents):
    ''' Calculate the required percentile from the hash count or frequency table given. 
    
    Can take a list of values for percentile, and will return an output for each'''
        
    import math
    
    if type(percents) != list:
        percents = [percents]
    
    N = sum(hashCount.values())
    sorted_keys = hashCount.keys()
    sorted_keys.sort()
    
    outValues = np.zeros(len(percents), dtype = np.uint8)
    
    for i, per in enumerate(percents):
         
        k = (N-1) * per + 1 # f and c must count from 1 not zero
        f = math.floor(k)  
        c = math.ceil(k) 
        
        cum_count = 0
        if f == c: # No remainder so no need to interpolate
            for key in sorted_keys:    
                cum_count += hashCount[key]
                if cum_count >= k:
                    # Found the right one 
                    break
            outValues[i] = key
    
        else: # Must now check to see if we need to interpolate
            floor_key = 0
            ceil_key = 0
            found_floor_val = 0
            for key in sorted_keys:    
                cum_count += hashCount[key]
                if not found_floor_val:
                    if cum_count >= f:  
                        # Found the right one 
                        floor_key = key
                        found_floor_val = 1
                if cum_count >= c: 
                    ceil_key = key
                    break
            
            if floor_key == ceil_key:
                outValues[i] = key
            else:
                # Interpolate
                val1 = floor_key * (c-k)
                val2 = ceil_key * (k-f)
                outValues[i] = val1 + val2     
    
    if len(outValues) == 1:
        return outValues[0]
    else:                        
        return outValues
               
def plotFTstats(stats):
    ''' Plot the min, max, Q1, Q2, Median and Average from a stats record structure '''
    
    f = plt.figure()
    ax = f.add_subplot(1,1,1)
    ax.plot(stats['max'], 'k-', label='max')
    ax.plot(stats['Q2'], 'b.', label = 'Q2')
    ax.plot(stats['median'], 'r.', label = 'median')
    ax.plot(stats['ave'], 'g-', label = 'average')
    ax.plot(stats['Q1'], 'b.', label='Q1')
    ax.plot(stats['min'], 'k-', label = 'min')
    ax.set_ylabel('Phred Score') 
    ax.set_xlabel('Base Position')
    ax.legend(loc=3)
    ax.set_title('Statistics for Phred Scores per Base Position')
    plt.draw()
     

def getMeanPhredPerBase(inFiles = None, fileType = '', dataPath = ''):
    ''' Find the mean Phred score per base '''
    
    #Generator for Sequence Records
    SeqRecGenerator = getGen4SeqRecFiles(inFiles = inFiles, fileType = fileType, dataPath = dataPath)
    numfiles = SeqRecGenerator.next() 

    # Make Data Structure
    running_aves = np.zeros(101)
    
    toc = time.time()
    cum_t = 0
    count = 1
    for fileNo, seqRecIdx in enumerate(SeqRecGenerator):
        for rec in seqRecIdx.itervalues():
            xi = np.array(rec.letter_annotations['phred_quality'])
            running_aves = running_aves + ((xi - running_aves)/ (count)) 
            count += 1
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished file {0} after {1}'.format(fileNo, time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
    
    total_t = time.time() - toc
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', time.gmtime(total_t)))
    
    return running_aves
          

if __name__ == '__main__':

#    dataLoacation = '/Users/chris/Datasets/gazellesZebras/'
#    dataLoacation = '/space/musselle/datasets/gazellesAndZebras/lane6'

#    getPropN_meanPhred(fileType = '*.idx', dataPath = dataLoacation)
#    getReadLengths(fileType = '*.idx', dataPath = dataLoacation)
    
    
#    dataLoacation = '/space/musselle/datasets/gazellesAndZebras/lane6'
#    lane6_PhredAves = getMeanPhredPerBase(fileType = '*.idx', dataPath = dataLoacation)
    
#    dataLoacation = '/space/musselle/datasets/gazellesAndZebras/lane6'
#    lane6_PhredStats, lane6_PhredCounts = boxPlotPhredPerBase(fileType = '*.bgzf', dataPath = dataLoacation)    
#
#    dataLoacation = '/space/musselle/datasets/gazellesAndZebras/lane8'
#    lane8_PhredStats, lane8_PhredCounts = boxPlotPhredPerBase(fileType = '*.bgzf', dataPath = dataLoacation)    
#
#    with open('L6_phredCounts.pkl', 'wb') as f:
#        pkl.dump(lane6_PhredCounts, f )
#    with open('L8_phredCounts.pkl', 'wb') as f:
#        pkl.dump(lane8_PhredCounts, f )
#
#    np.save('L6_phredStats', lane6_PhredStats)
#    np.save('L8_phredStats', lane8_PhredStats)


    dataLoacation = '/space/musselle/datasets/gazellesAndZebras/lane6'
    files = ['lane6_NoIndex_L006_R1_005.fastq.bgzf']
    readL1 = getReadLengths1(inFiles = files, dataPath=dataLoacation)
    readL2 = getReadLengths2(inFiles = files, dataPath=dataLoacation)
    readL3 = getReadLengths3(inFiles = files, dataPath=dataLoacation)

#    lane6_PhredStats = boxPlotPhredPerBase( files , dataPath = dataLoacation)  
    
      


