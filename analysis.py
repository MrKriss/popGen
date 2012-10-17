'''
Created on Oct 9, 2012

@author: chris
'''

from Bio import SeqIO
import numpy as np
import time
import os
import cPickle as pkl

def countRecords(filename):
  ''' Counts number of records in file'''

  S = SeqIO.parse(filename, 'fastq')
  count = 0
  for seq in S:
    count += 1
    
  return count

if __name__ == '__main__':
  
  dataLoacation = '/Users/chris/Datasets/gazellesZebras/'
  os.chdir(dataLoacation)
  
  fileNames = ['lane8_NoIndex_L008_R1_001.fastq', 
           'lane8_NoIndex_L008_R1_002.fastq',
           'lane8_NoIndex_L008_R1_003.fastq',
           'lane8_NoIndex_L008_R1_004.fastq',
           'lane8_NoIndex_L008_R1_005.fastq' ]
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
    S = SeqIO.parse(f, 'fastq')
    for c, seqR in enumerate(S):
      # Proportion of Ns
      outList[idx]['propN'][c] = float(seqR.seq.count('N')) / len(seqR)
      # Mean Phred Score
      outList[idx]['meanPhred'][c] = np.array(seqR.letter_annotations['phred_quality']).mean()
    
    loop_t = time.time() - tik - cum_t
    cum_t += loop_t 
    print 'Finished {0} after {1}'.format(f, time.strftime('%H:%M:%S', time.gmtime(loop_t)))

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
  print 'Saving Data.......'  
  D = {'files' : fileNames, 
       'numRecords' : recordCounts,
       'fileStats': outList,
       'totalStats' : totalStats}
  filename = fileNames[0].split('_')[0] + '-stats.pkl' 
  with open(filename, 'wb') as f:
    pkl.dump(D, f) 
  print 'DONE! Data saved to {0}'.format(os.getcwd() + filename)
  
  total_t = time.time() - toc
  print 'Processed {0} files in {1}'.format(numFiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))
  