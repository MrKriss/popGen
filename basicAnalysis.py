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
from utils.utils import pklSave

from utils.utils import Cycler

def getPropN_meanPhred(inFiles = None, fileType = '', dataPath = ''):
    ''' Calculate the proportion of Ns and the mean Phred scores of the files given 
    
    Returns a Histogram'''

    RecCycler = Cycler(inFiles = inFiles, fileType = fileType, dataPath = dataPath)
    
    print '\nCalculating Proportion of Ns and Mean Phred Score per read.\n'
    
    # Define vars and outputs
    numFiles = RecCycler.numFiles
    outList = [0] * numFiles
    recordCounts = [0] * numFiles

    # Find the number of records in each file.
    toc = time.time()
    cum_t = 0
    for seqRecGen in RecCycler.seqFileGen:
    
        fileName = RecCycler.curFileName
        fileNum = RecCycler.curFileNum
        
        # Check fileName is a .idx file extention or that one exists for it.
        if fileName.endswith('idx'):
            idxFileName = fileName
        else:
            assert os.path.isfile(fileName.split('.')[0] + '.idx'),  'Index file not found for {0}'.format(fileName)
            idxFileName = fileName.split('.')[0] + '.idx'
        
        # count number of Records 
        S = SeqIO.index_db(idxFileName, format='fastq')
        recordCounts[fileNum] = len(S)
        del S
        
        outList[fileNum] = {'propN'  :  np.zeros(recordCounts[fileNum]),
                        'meanPhred' : np.zeros(recordCounts[fileNum])  }

        # Cycle through all records 
        for num, seqRec in enumerate(seqRecGen):
            # Proportion of Ns
            outList[fileNum]['propN'][num] = float(seqRec.seq.count('N')) / len(seqRec)
            # Mean Phred Score
            outList[fileNum]['meanPhred'][num] = np.array(seqRec.letter_annotations['phred_quality']).mean()
        
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
    print '\nProcessed {0} files in {1}'.format(numFiles, time.strftime('%H:%M:%S', time.gmtime(total_t)))

    return totalStats
        
def getReadLengths(inFiles = None, fileType = '', dataPath = ''):
    ''' Return histogram of read lengths. Uses RecordCyler object'''        
    
    #Generator for Sequence Record files
    SeqRecCycler = Cycler(inFiles = inFiles, fileType = fileType, dataPath = dataPath)

    print 'Calculating length per read ...'
    
    # Define vars and outputs
    outList = [0] * SeqRecCycler.numFiles
    recordCounts = [0] * SeqRecCycler.numFiles

    toc = time.time()
    cum_t = 0
    for seqRecs in SeqRecCycler.fileGen:
        fileName = SeqRecCycler.fileGen.fileName
        fileNum = SeqRecCycler.fileGen.fileNum
        # count number of Records, may need to alter this if no idx file
        idxFileName = '.'.join(fileName.split('.')[:-2]) + '.idx'
        S = SeqIO.index_db(idxFileName, format='fastq')
        recordCounts[fileNum] = len(S)
        del S
        
        outList[fileNum] = {'length'  :  np.zeros(recordCounts[fileNum]) }
        
        for num, seqR in enumerate(seqRecs):
            outList[fileNum]['length'][num] = len(seqR.seq)
        
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

def boxPlotPhredPerBase(inFiles = None, fileType = '', dataPath = '', saveprefix = ''):
    ''' Find the median, upper and lower quartile for the Phred score per base 
    
    Returns the stats and the counter dictionary. 
    
    Counter dictionary may become standard way to store mass Phred/seq bases data.  '''
    from collections import Counter

    RecCycler = Cycler(inFiles = inFiles, fileType = fileType, dataPath = dataPath)
    
    print '\nCalculating Box plot stats of phred scores per base position.\n'
    
    # Define vars and outputs
    numFiles = RecCycler.numFiles

    toc = time.time()
    cum_t = 0
    
    counterList = [0] * 101
    for i in range(len(counterList)):
        counterList[i] = Counter()
    
    for seqRecGen in RecCycler.seqFileGen:
        
        fileName = RecCycler.curFileName
        fileNum = RecCycler.curFileNum
        
        for rec in seqRecGen:
            for baseNo, phred in enumerate(rec.letter_annotations['phred_quality']):
                counterList[baseNo][phred] += 1
                
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} \nfile {1} of {2} after {3}'.format(fileName, fileNum, numFiles, time.strftime('%H:%M:%S', time.gmtime(loop_t))) 

    # Calculate min, max Q1, Q2, Median and Average
    stats = getStatsFT(counterList)

    total_t = time.time() - toc
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', time.gmtime(total_t)))
    
    pklfilename = dataPath.split('/')[-1]
            
    pklSave(counterList, '_'.join([pklfilename, saveprefix , 'phredCount']))
    np.save( '_'.join([pklfilename, saveprefix , 'phredStats.npy']) , stats)
    
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
    RecCycler = Cycler(inFiles = inFiles, fileType = fileType, dataPath = dataPath)
    
    print '\nCalculating mean phred score per base position.\n'

    # Make Data Structure
    running_aves = np.zeros(101)
    
    toc = time.time()
    cum_t = 0
    count = 1
    for num, seqRecs in enumerate(RecCycler):
        for rec in seqRecs.itervalues():
            xi = np.array(rec.letter_annotations['phred_quality'])
            running_aves = running_aves + ((xi - running_aves) / (count)) 
            count += 1
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished file {0} after {1}'.format(num, time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
    
    total_t = time.time() - toc
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', time.gmtime(total_t)))
    
    return running_aves


if __name__ == '__main__':

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


    dataLoc = '/space/musselle/datasets/gazellesAndZebras'
    files = ['lane6_NoIndex_L006_R1_005.fastq.bgzf']

#    lane6_PhredStats = boxPlotPhredPerBase( files , dataPath = dataLoacation)  
    
      


