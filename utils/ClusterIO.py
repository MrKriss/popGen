'''
Created on 24 Apr 2013

@author: musselle
'''
import os 
import sys 
from collections import defaultdict, Counter

import numpy as np 


'''
IO support for CDHIT clustering files in python.

Features

    * Iterate through Clusters in a file
    * Return counters for cluster sizes in a file
    * Analise Distribution of sequences within cluster 
    (counts dictionary for how far sequences are away from the representative one)


  
'''


def mismatche2percentage(mismatches, seq_length):
    ''' Convert a number of mismatches to a percentage rounded to the nearest 2 dp '''
    return  round((float(mismatches) / seq_length) * 10000) / 100 

def percentage2mismatch(percentage, seq_length):
    ''' Convert a percentage similarity to a number of mismatches (rounded to integer)'''
    return  round( (percentage / 100.) * seq_length)



def parse(handle, ):
    """ Reads in a CDHIT cluster file and returns an iterator for the clusters. 
    
     - handle   - handle to the file, or the filename as a string
     - mode     - 

    Currently iterates in the order the file is read. 
    Typical usage, opening a file to read in, and looping over the record(s):
    """

    # input checks    
    if type(handle) is str:
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        
    # Setup Data structure?
        
    try:
        with handle as cluster_file:   
            
            for line in cluster_file:              
                line = line.strip()
                
                if line.startswith('>'):
                    # This is start of new cluster
                    if seq_in_cluster and rep_length:
                        # This number is the size of last cluster
                        # Store result
                        ds[str(rep_length)][str(seq_in_cluster)] += 1
                        seq_in_cluster = 0
                        rep_length = 0
                        
                elif line.endswith('*'): 
                    # This is the representative sequence for the cluster
                    rep_length = int(line.split()[1].strip('nt,'))
                    seq_in_cluster += 1
                else:
                    seq_in_cluster += 1
            
            # Got to end of file but still one more cluster to add
            ds[str(rep_length)][str(seq_in_cluster)] += 1
    
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        print "Finished Scanning cluster file."    
    
    
    
    
def summary_counter(self, infile, path='', mode='total', report=True):
    ''' Takes cluster file output by CD-Hit and produces two Counter dictionaries 
    
    output depends on mode specified:
    
    modes:
    counter_per_sequence_length = { 'sequence_length' : Counter(cluster_size) }
    total = Counter(cluster_sizes_for _all_sequences)
    
    '''

    if not infile.endswith('.clstr'):
        infile = infile + '.clstr'
    
    # Data structure to store cluster size info is a DefaultDictionary of Counter dictionaries.
    # ds = { 'seq_len' : Counter(cluster_size)  }
    # empty keys of ds are initialised with a Counter dictionary. 
    
    ds = defaultdict(Counter)

    seq_in_cluster = 0
    rep_length = 0

    print 'Generating cluster summary for  %s ...' % (infile)

    try:
        with open(os.path.join(path,infile), 'rb')  as cluster_file:   
            
            for line in cluster_file:              
                line = line.strip()
                
                if line.startswith('>'):
                    # This is start of new cluster
                    if seq_in_cluster and rep_length:
                        # This number is the size of last cluster
                        # Store result
                        ds[str(rep_length)][str(seq_in_cluster)] += 1
                        seq_in_cluster = 0
                        rep_length = 0
                        
                elif line.endswith('*'): 
                    # This is the representative sequence for the cluster
                    rep_length = int(line.split()[1].strip('nt,'))
                    seq_in_cluster += 1
                else:
                    seq_in_cluster += 1
            
            # Got to end of file but still one more cluster to add
            ds[str(rep_length)][str(seq_in_cluster)] += 1
    
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        print "Finished Scanning cluster file."
    
    # Construct total cluster size counter 
    total_cluster_size_counter = Counter()
    for v in ds.itervalues():
        total_cluster_size_counter.update(v)
        
    # Construct representative sequence length counter 
    seq_len_counter = Counter()
    for k, v in ds.iteritems():
        seq_len_counter[k] += sum(v.values()) 
    
    if report:
        print 'Top 5 Cluster Sizes: ', total_cluster_size_counter.most_common()[:5]
        print 'Top 5 Sequence Lengths: ', seq_len_counter.most_common()[:5]
    
    # Decide what to output    
    if mode == 'total':
        return total_cluster_size_counter
    elif mode == 'by_seqlen':
        return ds
    elif mode =='both':
        return total_cluster_size_counter, ds
        
    
def hist_counters(counters, labels=None, **kwargs):
    ''' Construct a series of histograms from a list of Counter Dictionarys '''
    
    import matplotlib.pyplot as plt
    
    if type(counters) is not list and type(counters) is not tuple:
        counters = [counters]
    
    if labels is not None:
        assert len(labels) == len(counters), "Number of labels must match number of counters."
    
    
    for i in range(len(counters)):
        data = np.array(list(counters[i].elements()), dtype = np.int)

        if labels:
            plt.hist(data, histtype='step', label=labels[i], **kwargs)
        else:
            
            plt.hist(data, histtype='step', label='Counter-'+str(i), **kwargs)
        plt.title("Cluster Size Distribution")
        plt.xlabel("Value")
        plt.ylabel("Frequency")

    plt.legend()
    plt.show()






if __name__ == '__main__':
    pass