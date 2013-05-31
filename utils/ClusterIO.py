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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


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

        if 'seq' in items and 'phred' in items:
            self.rep_seq = lookup_db[self.rep_seq_desc].seq.tostring()
            self.rep_phred = np.array(lookup_db[self.rep_seq_desc].letter_annotations['phred_quality'])
            for elem in self.members_desc:
                self.members_seq.append(lookup_db[elem].seq.tostring())
                self.members_phred.append(np.array(lookup_db[elem].letter_annotations['phred_quality']))

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
                    
    def get_unique_seq(self, ignoreup2=6, lookup_db=None):
        """ Work out the counts for unique reads within a cluster """
        
        if not self.members_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(self, ['seq'], lookup_db=lookup_db)
        if not self.rep_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(self, ['rep'], lookup_db=lookup_db)
        
        unique_seq_counter = Counter()
        unique_seq_counter[self.rep_seq[ignoreup2:]] += 1
        
        seqs = [s[ignoreup2:] for s in self.members_seq]        
        unique_seq_counter.update(Counter(seqs))
        self.unique_seq = unique_seq_counter
    
    def get_basefraction(self, lookup_db=None):
        """ Calculate the fraction of nucleotide bases per base location """
        
        # Make sure seq data is available
        if not self.members_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['seq'], lookup_db=lookup_db)
        if not self.rep_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['rep'], lookup_db=lookup_db)
        
        self.basefrac = {} 
        self.basefrac['A'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['T'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['G'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['C'] = np.zeros(len(self.rep_seq)) 
        self.basefrac['N'] = np.zeros(len(self.rep_seq)) 
        
        # Cheack representative sequence
        for i in range(len(self.rep_seq)):
            self.basefrac[self.rep_seq[i]][i] += 1
            
        # Check Cluster members
        for seq in self.members_seq:
            for i in range(len(self.rep_seq)):  
                self.basefrac[seq[i]][i] += 1
        
        for i in range(len(self.rep_seq)):
            self.basefrac['A'][i] /= float(self.size)
            self.basefrac['T'][i] /= float(self.size)
            self.basefrac['G'][i] /= float(self.size)
            self.basefrac['C'][i] /= float(self.size)
            self.basefrac['N'][i] /= float(self.size)
            
    def plot_basefrac(self, ignorefirst=6, xlab="", ylab="", title="", **kwargs):
        """ Visualise the fraction of bases """
        
        import matplotlib.pyplot as plt
    
        # Default plot options
        if 'ms' not in kwargs:
            kwargs['ms'] = 10.0
        if 'marker' not in kwargs:
            kwargs['marker'] = '.'
        if 'mew' not in kwargs:
            kwargs['mew'] = 0.0
    
        plt.figure()
        for k,v in self.basefrac.iteritems():
        
            data_xs = range(1,len(v[ignorefirst:])+1) 
            data_ys = v[ignorefirst:]
            label =  k 
            
            vars()['line' + k], = plt.plot(data_xs, data_ys, label=label, ls='', **kwargs)
    
        # Shink current axis's height by 10% on the bottom
        ax = plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
    
        # Put a legend below current axis
        
        ax.legend([vars()['lineA'], vars()['lineT'], vars()['lineG'], vars()['lineC'], vars()['lineN']],
                  ['A', 'T', 'G', 'C', 'N'], loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=5, numpoints=1, markerscale=3)
        
        plt.ylim(-0.1, 1.1)
        plt.xlim(1, 101)
        plt.title(title)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.show()
        
    def write2fastq(self):
        """ Write cluster to a fastq file format. First seq is representative seq for the cluster """
        pass

def filter_clusters(cluster_filepath, idx_filepath, size_range, output_dirpath):
    """ Writes a subset of cluster sizes to FastQ files 
    
    The representative sequence is the first sequence record written.
    
    """
    
    starting_dir = os.getcwd()
    idx_dir = os.path.split(idx_filepath)[0]
    
    # Check and create directory 
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
        
    cluster_gen = parse(cluster_filepath, idx_filepath)
    seqrec_lookup = SeqIO.index_db(idx_filepath)
    
    size_counter = Counter()
    
    for cluster in cluster_gen:
        # Check if cluster size is within defined range
        if cluster.size >= size_range[0] and cluster.size < size_range[1]:
            
            size_counter[cluster.size] += 1
            
            # Get the sequence records for the cluster 
            seqs = []
            if os.getcwd() != idx_dir:
                os.chdir(idx_dir)
            seqs.append(seqrec_lookup[cluster.rep_seq_desc])
            for member in cluster.members_desc:
                seqs.append(seqrec_lookup[member])
            
            if os.getcwd() != output_dirpath:
                os.chdir(output_dirpath)
            # Write cluster to a file 
            fname = "clustersize{0}-No{1}.fastq".format(str(cluster.size), str(size_counter[cluster.size]))
            output_handle = open(fname, "wb")
            SeqIO.write(seqs, output_handle, "fastq")
        
def filter_clusters2(cluster_filepath, idx_filepath, size_range, output_dirpath):
    """ Writes a subset of cluster sizes to FastQ files 
    
    The representative sequence is the first sequence record written.
    
    make the sequence record instead of passing it.
    
    """
    
    starting_dir = os.getcwd()
    idx_dir = os.path.split(idx_filepath)[0]
    
    # Check and create directory 
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
        
    cluster_gen = parse(cluster_filepath, idx_filepath)
    seqrec_lookup = SeqIO.index_db(idx_filepath)
    
    size_counter = Counter()
    
    for cluster in cluster_gen:
        # Check if cluster size is within defined range
        if cluster.size > size_range[0] and cluster.size < size_range[1]:
            
            size_counter[cluster.size] += 1
            
            # Get the sequence records for the cluster 
            if os.getcwd() != idx_dir:
                os.chdir(idx_dir)
            # Representative sequence first 
            seqrecord = seqrec_lookup[cluster.rep_seq_desc]
            
            if os.getcwd() != output_dirpath:
                os.chdir(output_dirpath)
            # Write cluster to a file 
            fname = "clustersize{0}-No{1}.fastq".format(str(cluster.size), str(size_counter[cluster.size]))
            
            if os.path.isfile(fname):
                output_handle = open(fname, "wb")
                output_handle.close()
            
            output_handle = open(fname, "a")
            SeqIO.write(seqrecord, output_handle, "fastq")
            
            for member in cluster.members_desc :
                
                if os.getcwd() != idx_dir:
                    os.chdir(idx_dir)
                # Representative sequence first 
                seqrecord = seqrec_lookup[member]
            
                if os.getcwd() != output_dirpath:
                    os.chdir(output_dirpath)
                # Write sequence record to file 
                SeqIO.write(seqrecord, output_handle, "fastq")
                
    if os.getcwd() != starting_dir: 
        os.chdir(starting_dir)
        
def filter_subprocess(infile, file_sorted=False, mode=None):
    """" Filtering Steps to obtain clusters in a certain size range and % range of all reads """
    
    # Sort clusters by reads in clusters of that size 
    if not file_sorted:
        infile = sortby(infile, reverse=True, mode='reads_per_cluster')
    
    # Return size only in parsing clusters  

    for c in parse(infile):
        pass


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




"""
Impliment a double pass generator to parse reads:

1. Get the next cluster size, and store start and end places. 

2. If size passes, go back and rescan file to get detailed info for that cluster.

"""




def parse(handle, idx_file_path="", mode="size_only"):
    """ Reads in a CDHIT cluster file and returns a generator for the clusters. 
    
     - handle   - handle to the file, or the filename as a string
     - idx_file - full path to file containing the index to the original sequence records.
                  If given, the file is stored and passed to each Cluster object generated
                  for possible future lookup. 
    - mode:
        'size_only'   =   Only return the cluster size each iteration.

    Currently iterates in the order the file is read. 
    Typical usage, opening a file to read in, and looping over the record(s):
    """
    
    # input checks    
    if type(handle) is str:
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        if handle.closed:
            handle = open(handle.name, handle.mode)
        
    if idx_file_path:
        assert type(idx_file_path) is str and os.path.exists(idx_file_path), 'Invalid path.'
        
    first = True
    # Setup Data structure
    cluster = Cluster()
    cluster.idx_file_path = idx_file_path
    
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
                        # Reset Cluster
                        cluster = Cluster()
                        cluster.idx_file_path = idx_file_path
                    
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
                  "reads_per_cluster" = Sort by total reads in clusters of each size

    Works by scanning file once through to build an index, sorting the index, then 
    rewriting file accordingly. 
    """

    if mode == 'reads_per_cluster':
        reads_per_cluster_counter = summary_counter(handle, mode='reads_per_cluster', report=True)
        
    # Sorting key functions
    key_ops = {
               'cluster_size': lambda x : x[0],
               'reads_per_cluster': lambda x : reads_per_cluster_counter[str(x[0])], 
               }

    # Input checks    
    if type(handle) is str:
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        if handle.closed:
            handle = open(handle.name, handle.mode)

    # Setup Data structure and utility vars
    # cluster_idxs = [ (size, start_bytes, end_bytes ), ... ]
    cluster_idxs = []
    start_time = time.time()
    first = True
    # First cycle through file and store all cluster start and end positions 
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
            
            # Sort index in the desired fashion 
            sorted_cluster_idxs = sorted(cluster_idxs, key=key_ops[mode], reverse=reverse)    
                
            # Rewrite file
            old_filename_parts = cluster_file.name.split('.')
            old_filename_parts[0]  += '-sortedby_' + mode
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

    return handle, sorted_cluster_idxs, cluster_idxs

    
def getfromdb(cluster_list, items, lookup_db=None, path2idxfile=None,):
    """ Adds the specified items to each cluster in list, looked up from the indexed db.
    
    items - list of things to fetch for the cluster. Option of either/both 'seq' and 'phred'
            'rep' will fetch just the seq and phred info for the representative sequence.
            
    """

    if lookup_db is None:
        assert path2idxfile, 'No lookup data base parsed or index file specified'
        lookup_db = SeqIO.index_db(path2idxfile)  
    
    for cluster in cluster_list:
        cluster.getfromdb(items, lookup_db)
        
    return cluster_list


def summary_counter(handle, mode='cluster_size', report=True):
    ''' Takes cluster file output by CD-Hit and produces two Counter dictionaries 
    output depends on mode specified:
    
    - handle   - handle to the file, or the filename as a string.
    - mode:
        'by_seqlen'         = { 'sequence_length' : Counter(cluster_size) }
        'cluster_size'      = Counter(cluster_sizes_for_all_sequence_lengths)
        'reads_per_cluster' = { 'cluster size'  : cluster_size * Num clusters }
    
    '''

    # input checks    
    if type(handle) is str:
        if not handle.endswith('.clstr'):
            handle = handle + '.clstr'
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        if handle.closed:
            handle = open(handle.name, handle.mode)

    # Data structure to store cluster size info is a DefaultDictionary of Counter dictionaries.
    # ds = { 'seq_len' : Counter(cluster_size)  }
    # empty keys of ds are initialised with a Counter dictionary. 
    
    ds = defaultdict(Counter)

    seq_in_cluster = 0
    rep_length = 0

    print 'Generating cluster summary for  %s ...' % (handle.name)

    try:
        with handle  as cluster_file:   
            
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
    if mode == 'cluster_size':
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
    gaz_clustf = '/space/musselle/data/RAD-seq/gazelles-zebras/clusters/gz_allg_allz_95g1_clustered_reads_c95_g1/gz_allg_allz-gazelle-clustered_c95_g1.clstr'
    zeb_clustf = '/space/musselle/data/RAD-seq/gazelles-zebras/clusters/gz_allg_allz_95g1_clustered_reads_c95_g1/gz_allg_allz-zebra-clustered_c95_g1.clstr'
    idx_file = '/space/musselle/data/RAD-seq/gazelles-zebras/processed-data/all_reads_fastq.idx'

    output_dir = '/space/musselle/data/RAD-seq/gazelles-zebras/clusters/gz_allg_allz_95g1_clustered_reads_c95_g1/gazellles_filtered/'

    filter_clusters(gaz_clustf, idx_filepath=idx_file, size_range=(300,300), output_dirpath=output_dir)
#     filter_clusters2(gaz_clustf, idx_filepath=idx_file, size_range=(300,310), output_dirpath=output_dir)
