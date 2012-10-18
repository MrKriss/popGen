'''
Created on Oct 9, 2012

@author: chris
'''

from Bio import SeqIO
import numpy as np
import time
import os
import cPickle as pkl
import gzip


def countRecords(fileName):
    ''' Counts number of records in file'''

    idxFileName = '.'.join(fileName.split('.')[:-2]) + '.idx'
    S = SeqIO.index_db(idxFileName, format='fastq')

    return   len(S)

if __name__ == '__main__':

#    dataLoacation = '/Users/chris/Datasets/gazellesZebras/'
    dataLoacation = '/space/musselle/datasets/gazellesAndZebras/'
    os.chdir(dataLoacation)

#    filesLane8_1 = ["lane8/lane8_NoIndex_L008_R1_00%i.fastq.gz" % (i+1) for i in range(9)]
#    filesLane8_2 = ["lane8/lane8_NoIndex_L008_R1_0%i.fastq.gz" % (i+10) for i in range(1)]
    filesLane6_1 = ["lane6/lane6_NoIndex_L006_R1_00%i.fastq.gz" % (i+1) for i in range(9)]
    filesLane6_2 = ["lane6/lane6_NoIndex_L006_R1_0%i.fastq.gz" % (i+10) for i in range(2)]
#    fileNames = filesLane8_1 + filesLane8_2 # + filesLane6_1 + filesLane6_2
    fileNames = filesLane6_1 + filesLane6_2
    numFiles = len(fileNames)


    toc = time.time()
    recordCounts = map(countRecords, fileNames)
    tek = time.time() - toc

    print 'Finished counting records for {0} files in {1}'.format(numFiles,
                time.strftime('%H:%M:%S', time.gmtime(tek)))

    # Construct Data structures
    outList = [0] * numFiles
    for i in range(numFiles):
        outList[i] = {'propN':  np.zeros(recordCounts[i])  ,
                      'meanPhred': np.zeros(recordCounts[i]) }

    # Fill with stats
    tik = time.time()
    cum_t = 0
    for idx, f in enumerate(fileNames):
        # Change Files to .bgzf extenstion
        bgzfFileName = '.'.join(f.split('.')[:-1]) + '.bgzf'
        with gzip.open(bgzfFileName, 'r') as handle:
            S = SeqIO.parse(handle, 'fastq')
            for c, seqR in enumerate(S):
                # Proportion of Ns
                outList[idx]['propN'][c] = float(seqR.seq.count('N')) / len(seqR)
                # Mean Phred Score
                outList[idx]['meanPhred'][c] = np.array(seqR.letter_annotations['phred_quality']).mean()

        loop_t = time.time() - tik - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(bgzfFileName, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    # Construct and fill total Records stat object
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

    # Pickle data dictionary
    print 'Saving Data......'
    np.save('lane6_propN', totalStats['propN'])
    np.save('lane6_meanPhred', totalStats['meanPhred'])
    print 'DONE! Data saved to {0}'.format(os.getcwd())

    total_t = time.time() - toc
    print 'Processed {0} files in {1}'.format(numFiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))
