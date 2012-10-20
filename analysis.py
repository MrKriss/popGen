'''
Created on Oct 9, 2012

@author: chris
'''

from Bio import SeqIO
import numpy as np
import time
import os

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

    # Define vars and outputs
    numFiles = len(inFiles)
    outList = [0] * numFiles
    recordCounts = [0] * numFiles

    # Find the number of records in each file.
    toc = time.time()
    cum_t = 0
    for fileNo, fileName in enumerate(inFiles):
        # Check its a .idx file extention or that one exists with it.
        if fileName.split('.')[-1] == '.idx':
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
    np.save(parent_dir + 'propN', totalStats['propN'])
    np.save(parent_dir + 'lane6_meanPhred', totalStats['meanPhred'])
    print 'DONE! Data saved to {0}'.format(os.getcwd())

    total_t = time.time() - toc
    print 'Processed {0} files in {1}'.format(numFiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))

    return totalStats



if __name__ == '__main__':

#    dataLoacation = '/Users/chris/Datasets/gazellesZebras/'
    dataLoacation = '/space/musselle/datasets/gazellesAndZebras/lane6'

    getPropN_meanPhred(fileType = '*.idx', dataPath = dataLoacation)
    
 