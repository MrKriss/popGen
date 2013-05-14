'''
Created on 24 Apr 2013

@author: musselle
'''
import os 
import sys 
import time
from collections import defaultdict, Counter

import numpy as np 
from Bio import SeqIO

'''
IO support for CDHIT clustering files in python.

Features

    * Iterate through Clusters in a file
    * Return counters for cluster sizes in a file
    * Analise Distribution of sequences within cluster 
    (counts dictionary for how far sequences are away from the representative one)

Sample Cluster File:

>Cluster 0
0    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:1148:2096... *
1    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:2328:25159... at +/100.00%
2    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:7411:52830... at +/98.86%
3    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:9529:76269... at +/100.00%
4    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1102:6206:12309... at +/100.00%
5    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1102:8696:31554... at +/100.00%
6    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1103:7706:11162... at +/100.00%
7    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1104:9565:35059... at +/98.86%
8    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1104:15657:55944... at +/100.00%
9    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1105:13547:22466... at +/98.86%
10    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1105:20059:45463... at +/97.73%

'''
"""IDEAS:

Representative sequence record

member descriptions 
index file to lookup sequence records of members

counter of unique clusters 

method to plot frequency of site 

Object tied to file or clusters within file? 

"""

class Cluster(object):
    """ Holds all cluster based information. """
    
    def __init__(self):
        
        # Cluster vars
        self.rep_seq_desc = ""    
        self.rep_seq = ""
        self.rep_phred = None
        self.members_desc = []    
        self.members_seq = []
        self.members_phred = []
        self.size = 0  
        self.id = 0
        self.edit_dists = []
        
        # Start and end Locations on the cluster file in bytes from the beginning  
        self.start_loc = 0
        self.end_loc = 0
            
        self.idx_file_path = ""
        
    def getfromdb(self, items, lookup_db=None):
        """ Lookup to indexed db to retrieve and store the specified items.
        
        items - list of things to fetch for the cluster. Option of either/both 'seq' and 'phred'
                'rep' will fetch just the seq and phred info for the representative sequence.
        
        """
        
        assert 'seq' in items or 'phred' in items or 'rep' in items, "Invalid values for items to lookup"
        
        if not self.idx_file_path:
            raise Exception('No idx_file_path specified.')
        
        if lookup_db is None: 
            print "Loading {} ...".format(os.path.split(self.idx_file_path)[1])
            lookup_db = SeqIO.index_db(self.idx_file_path)  
            print "Loading complete"
        
        start_dir = os.getcwd()
        path = os.path.split(self.idx_file_path)[0]
        
        if start_dir != path:
            os.chdir(path) 

        if 'seq' in items:
            self.rep_seq = lookup_db[self.rep_seq_desc].seq.tostring()
            for elem in self.members_desc:
                self.members_seq.append(lookup_db[elem].seq.tostring())

        if 'phred' in items:
            self.rep_phred = np.array(lookup_db[self.rep_seq_desc].letter_annotations['phred_quality'])
            for elem in self.members_desc:
                self.members_phred.append(np.array(lookup_db[elem].letter_annotations['phred_quality']))

        if 'rep' in items: # Just fetch the representitive sequence info 
            self.rep_seq = lookup_db[self.rep_seq_desc].seq.tostring()
            self.rep_phred = np.array(lookup_db[self.rep_seq_desc].letter_annotations['phred_quality'])
        
        if os.getcwd() != start_dir:
            os.chdir(start_dir)
        
    
def gettop_clusters(clusterfile, x, lower_limit, upper_limit):
    """ Return the top x cluster sizes which have most reads"""
    
    # Get the top cluster sizes that contain the top most x reads
    reads_counter = summary_counter(clusterfile, mode='reads_per_cluster', report=False)
    top_reads_per_cluster = reads_counter.most_common(x)


    for cluster_size in top_reads_per_cluster.iterkeys():
        
        if cluster_size > upper_limit or cluster_size < lower_limit:
            continue
        
        # Now have a target cluster size to work with.
        
        # Generator to return clusters of this size? 
        
    
    
    pass
    """
    Get summary counter for cluster file.
    
    Find most representative clusters/ hump of the clusters
    
     
    """
    
    
    


def mismatche2percentage(mismatches, seq_length):
    ''' Convert a number of mismatches to a percentage rounded to the nearest 2 dp '''
    mismatches = int(mismatches)
    seq_length = int(seq_length)
    return  round((float(mismatches) / seq_length) * 10000) / 100 

def percentage2mismatch(percentage, seq_length):
    ''' Convert a percentage similarity to a number of mismatches (rounded to integer)'''
    percentage = float(percentage)
    seq_length = int(seq_length)
    return  int(round( (percentage / 100.) * seq_length))

def parse(handle, idx_file_path=""):
    """ Reads in a CDHIT cluster file and returns a generator for the clusters. 
    
     - handle   - handle to the file, or the filename as a string
     - idx_file - full path to file containing the index to the original sequence records.
                  If given, the file is stored and passed to each Cluster object generated
                  for possible future lookup. 

    Currently iterates in the order the file is read. 
    Typical usage, opening a file to read in, and looping over the record(s):
    """

    # Setup Data structure
    cluster = Cluster()
    
    # input checks    
    if type(handle) is str:
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        
    if idx_file_path:
        assert type(idx_file_path) is str and os.path.exists(idx_file_path), 'Invalid path.'
        cluster.idx_file_path = idx_file_path
        
    first = True
    
    # Cycle through file     
    try:
        with handle as cluster_file:   
            
            for line in cluster_file:              
                line = line.strip() # Remove white space at start and end
                
                if line.startswith('>'):
                    # This is start of new cluster
                    
                    # Claculate stats for last cluster
                    cluster.size = len(cluster.members_desc) + 1
                    
                    # Store number for next yield 
                    next_cluster_num = line.split()[1]
                    
                    # yeild this cluster
                    if first:
                        first = False
                    else: 
                        yield cluster
                    
                elif line.endswith('*'): 
                    # Record Cluster number as this section only gets executed once per cluster
                    cluster.id = next_cluster_num
                    
                    # This is the representative sequence for the cluster
                    cluster.rep_seq_desc = line.split()[2].strip('>.')
                else:
                    
                    line_parts = line.split()
                    
                    cluster.members_desc.append(line_parts[2].strip('>.'))
                    similarity = line_parts[4].strip('+/%')
                    seq_len = line_parts[1].strip('nt,')
                    cluster.edit_dists.append(percentage2mismatch( 100 - float(similarity), seq_len))
                    
            # Got to end of file but still one more cluster to add
            cluster.size = len(cluster.members_desc) + 1
            yield cluster
    
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        print "Finished Scanning cluster file."    
    
def sortby(handle, reverse=True, mode='cluster_size'):
    """ Reads in a CDHIT cluster file and writes a sorted file for the clusters.
    by their size.  
    
     - handle   - handle to the file, or the filename as a string.
     - reverse  - Python sorts in ascending order by default, if reverse is true
                  the sort will instead be in descending order.  
     - mode     - "cluster_size" = Sort by the cluster sizes
                  "reads_per_clust_size" = Sort by total reads in clusters of each size

    Works by scanning file once through to build an index, sorting the index, then 
    rewriting file accordingly. 
    """

    # input checks    
    if type(handle) is str:
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        
    # Setup Data structure
    # cluster_idxs = [ (size, start_bytes, end_bytes ), ... ]
    cluster_idxs = []
    
    start_time = time.time()
    
    first = True
    # First cycle through file     
    try:
        with handle as cluster_file:   
            
            size = 0
            clust_end = 0
            clust_start = cluster_file.tell()
            
            while True:
                line = cluster_file.readline()  # Reading one line from file including end-of-line char '\n'
                if not line: 
                    break
                
                if first:
                    first = False
                    continue
                
                if line.startswith('>'):
                    # This is start of a new cluster
                    
                    # Record results for last cluster 
                    cluster_idxs.append((size, clust_start, clust_end))
                    # Reset cluster start and size counter
                    clust_start = clust_end
                    size = 0
                else:
                    size += 1
                    clust_end = cluster_file.tell()
                    
            # Got to end of file but still one more cluster to add
            cluster_idxs.append((size, clust_start, clust_end)) 
            
            # Sort index 
            if mode == "cluster_size":
                sorted_cluster_idxs = sorted(cluster_idxs, key=lambda x : x[0], reverse=reverse)    
    
            if mode == "reads_per_clust_size":
                pass
    
            # Rewrite file
            old_filename_parts = cluster_file.name.split('.')
            old_filename_parts[0]  += '-sorted'
            new_filename = '.'.join(old_filename_parts)
            with open(new_filename, 'wb') as sorted_cluster_file :
                
                for tup in sorted_cluster_idxs:
                    
                    cluster_file.seek(tup[1])
                    cluster_text = cluster_file.read(tup[2] - tup[1])
                    sorted_cluster_file.write(cluster_text)
                    
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        t = time.time() - start_time
        print "Finished sorting cluster file after {0}\n".format(time.strftime('%H:%M:%S', time.gmtime(t)))    

    return sorted_cluster_idxs, cluster_idxs

    
def get_seqrec4descriptions(cluster_list, path2idxfile):
    """Return a cluster dictionary containing the sequence record objects for the 
    representative sequence and all cluster members.
    
    Only return minimal info for now incase run into memory issues with large data volume.
    
    rewrote to use cluster classs
    """
    
    seq_record_db = SeqIO.index_db(path2idxfile)  
    
    cluster_seqrecs_list = []
    
    for cluster in cluster_list:
    
        newclust = Cluster()
        newclust.rep_seq = seq_record_db[cluster.rep_seq_desc]
        
        for elem in cluster.members_desc:
            newclust.members_desc.append(seq_record_db[elem])
        
        cluster_seqrecs_list.append(newclust)
        
    return cluster_seqrecs_list

def get_seqrecs(cluster_list, path2idxfile):
    """Return a cluster dictionary containing the sequence record objects for the 
    representative sequence and all cluster members.
    
    Only return minimal info for now incase run into memory issues with large data volume.
    
    rewrote to use cluster classs
    """
    
    start_dir = os.getcwd()
    path = os.path.split(path2idxfile)[0]
    
    if start_dir != path:
        os.chdir(path) 
    
    print "Loading {} ...".format(os.path.split(path2idxfile)[1])
    seq_record_db = SeqIO.index_db(path2idxfile)  
    print "Loading complete"
    
    cluster_seqrecs_list = []
    
    for cluster in cluster_list:
    
        newclust = Cluster()
        newclust.rep_seq = seq_record_db[cluster.rep_seq_desc]
        
        for elem in cluster.members_desc:
            newclust.members.append(seq_record_db[elem].seq.tostring())
        
        cluster_seqrecs_list.append(newclust)
    
    if os.getcwd() != start_dir:
        os.chdir(start_dir)
    
    return cluster_seqrecs_list


def cluster_unique_read_counter(cluster_obj, path2idxfile):
    """ Work out the counts for unique reads within a cluster """
    
    cluster_seqrecs = get_seqrec4descriptions(cluster_obj, path2idxfile)
    
    unique_reads_counter = Counter()
    
    unique_reads_counter[cluster_seqrecs['rep_seq']] += 1

    unique_reads_counter.update(Counter(cluster_seqrecs['members']))
    
    return unique_reads_counter


def summary_counter(cluster_path2file, mode='total', report=True):
    ''' Takes cluster file output by CD-Hit and produces two Counter dictionaries 
    
    output depends on mode specified:
    
    modes:
    counter_per_sequence_length = { 'sequence_length' : Counter(cluster_size) }
    total = Counter(cluster_sizes_for _all_sequences)
    reads_per_cluster = { 'cluster size'  : cluster_size * Num clusters }
    '''

    if not cluster_path2file.endswith('.clstr'):
        cluster_path2file = cluster_path2file + '.clstr'
    
    # Data structure to store cluster size info is a DefaultDictionary of Counter dictionaries.
    # ds = { 'seq_len' : Counter(cluster_size)  }
    # empty keys of ds are initialised with a Counter dictionary. 
    
    ds = defaultdict(Counter)

    seq_in_cluster = 0
    rep_length = 0

    print 'Generating cluster summary for  %s ...' % (os.path.split(cluster_path2file)[1])

    try:
        with open(cluster_path2file, 'rb')  as cluster_file:   
            
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
    elif mode == 'reads_per_cluster':
        reads_per_cluster = {}
        for k,v in total_cluster_size_counter.iteritems():
            reads_per_cluster[k] = int(k) * v
        return reads_per_cluster
    
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

def plot_counters(counters, labels=None, log='xy', xlab="", ylab="", title="", **kwargs):
    ''' Construct a series of scatter plots from a list of Counter Dictionaries '''
    
    import matplotlib.pyplot as plt
    
    # Default plot options
    if 'ms' not in kwargs:
        kwargs['ms'] = 4.0
    if 'marker' not in kwargs:
        kwargs['marker'] = '.'
    if 'mew' not in kwargs:
        kwargs['mew'] = 0.0
    
    if type(counters) is not list and type(counters) is not tuple:
        counters = [counters]
    
    if labels is not None:
        assert len(labels) == len(counters), "Number of labels must match number of counters."
    
    for i in range(len(counters)):
        
        data_xs = [int(k) for k in counters[i].keys()]
        data_ys = counters[i].values()

        if labels:
            plt.plot(data_xs, data_ys, label=labels[i], ls='', **kwargs)
        else:
            plt.plot(data_xs, data_ys, label='Counter-'+str(i), ls='', **kwargs)
        
        ax = plt.gca()
        if 'x' in log:
            ax.set_xscale('log')
        if 'y' in log:
            ax.set_yscale('log')
        
        plt.title(title)
        plt.xlabel(xlab)
        plt.ylabel(ylab)

    plt.legend(numpoints=1, markerscale=8)
    plt.show()


if __name__ == '__main__':
    pass