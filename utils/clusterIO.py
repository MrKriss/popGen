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

from editdist import distance

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

IMPROVEMENT
-----------

seqid is now an integer not the long text description. The seqid can be 
looked up in the main database to retrieve the information contained in 
the description if necessaery. 

No longer need a index file to lookup the file positions in compressed .bgzf files. 
All seq info is contained within a single database.  

'''

class ClusterObj(object):
    """ Holds all cluster based information. """
    
    def __init__(self):
        
        # Cluster vars
        self.rep_seq_id = ""      
        self.rep_seq = ""
        self.rep_phred = None
        self.members_id = []    
        self.members_seq = []
        self.members_phred = []
        self.size = 0  
        self.id = 0
        self.edit_dists = []
        
        # Start and end Locations on the cluster file in bytes from the beginning  
        self.start_loc = 0
        self.end_loc = 0
        
    def getfromdb(self, items, db):
        """ Lookup to main db to retrieve and store the specified items.
        
        items - list of things to fetch for the cluster. Option of either/both 'seq' and 'phred'
                'rep' will fetch just the seq and phred info for the representative sequence.
        
        db    - Reference to a Reads_db database object.
        """
        
        assert 'seq' in items or 'phred' in items or 'rep' in items, "Invalid values for items to lookup"
        
        start_dir = os.getcwd()
        path = os.path.split(db.dbfilepath)[0]
        
        if start_dir != path:
            os.chdir(path) 

        # Setup query
        sql_query = """ SELECT (seq, phred) FROM seqs WHERE seqid = ? """
        if 'seq' in items and 'phred' in items:
            
            # Get rep seq and phred
            record_curs = db.execute(sql_query, (self.rep_seq_id,))
            record = record_curs.fetchone()
            
            self.rep_seq = record['seq']
            phred_ascii = record['phred']
            phred_list = [ord(c) - 33 for c in phred_ascii]
            self.rep_phred = np.array(phred_list)

            # get members seq and phred
            for seqid in self.members_id:
                
                record_curs = db.execute(sql_query, (seqid,))
                record = record_curs.fetchone()
                
                self.members_seq.append(record['seq'])
                phred_ascii = record['phred']
                phred_list = [ord(c) - 33 for c in phred_ascii]
                self.members_phred.append(np.array(phred_list))

        if 'seq' in items:
            # Get rep seq
            record_curs = db.execute(sql_query, (self.rep_seq_id,))
            record = record_curs.fetchone()
            self.rep_seq = record['seq']
            
            # get members seq 
            for seqid in self.members_id:
                record_curs = db.execute(sql_query, (seqid,))
                record = record_curs.fetchone()
                self.members_seq.append(record['seq'])

        if 'phred' in items:
           
            # Get rep phred
            record_curs = db.execute(sql_query, (self.rep_seq_id,))
            record = record_curs.fetchone()
            phred_ascii = record['phred']
            phred_list = [ord(c) - 33 for c in phred_ascii]
            self.rep_phred = np.array(phred_list)

            # get members phred
            for seqid in self.members_id:
                record_curs = db.execute(sql_query, (seqid,))
                record = record_curs.fetchone()
                phred_ascii = record['phred']
                phred_list = [ord(c) - 33 for c in phred_ascii]
                self.members_phred.append(np.array(phred_list))

        if 'rep' in items: # Just fetch the representitive sequence info 
            # Get rep seq and phred
            record_curs = db.execute(sql_query, (self.rep_seq_id,))
            record = record_curs.fetchone()
            
            self.rep_seq = record['seq']
            phred_ascii = record['phred']
            phred_list = [ord(c) - 33 for c in phred_ascii]
            self.rep_phred = np.array(phred_list)
        
    def get_unique_seq(self, ignoreup2=6, db=None):
        """ Work out the counts for unique reads within a cluster """
        
        if not self.members_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['seq'], db=db)
        if not self.rep_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['rep'], db=db)
        
        unique_seq_counter = Counter()
        unique_seq_counter[self.rep_seq[ignoreup2:]] += 1
        
        seqs = [s[ignoreup2:] for s in self.members_seq]        
        unique_seq_counter.update(Counter(seqs))
        self.unique_seq = unique_seq_counter
    
    def get_basefraction(self, db=None):
        """ Calculate the fraction of nucleotide bases per base location """
        
        # Make sure seq data is available
        if not self.members_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['seq'], db=db)
        if not self.rep_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['rep'], db=db)
        
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
        
        
    def write2db(self, db):
        ''' Call the db method to load cluster into database. '''
        db.load_single_clusterobj(self)
        
    def write2clstr(self, handle, db=None):
        """ Write cluster to a .clstr file format like CD-HIT. 
        
        e.g.
        
        >Cluster 0
        0    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:1148:2096... *
        1    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:2328:25159... at +/100.00%
        2    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:7411:52830... at +/98.86%
        3    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1101:9529:76269... at +/100.00%
        4    88nt, >HWI-ST0747:233:C0RH3ACXX:6:1102:6206:12309... at +/100.00%
        
        """
        
        if type(handle) is str:
            if not handle.endswith('.clstr'):
                filename = handle + '.clstr'
            handle = open(handle, 'a')
        else:
            assert type(handle) is file, 'Invalid file handle.'
            if handle.closed:
                handle = open(handle.name, 'a')
        
        # Make sure seq data is available
        if not self.members_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['seq'], db=db)
        if not self.rep_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(['rep'], db=db)
        
        # Derived vars
        seq_len = str(len(self.rep_seq) - 13)
        
        with handle as clust_file:
            
            # Write Header 
            clust_file.write('>Cluster {0}\n'.format(self.id))
            
            # Write representative sequence
            clust_file.write('0\t{length}nt, >{rep_seq_id}... *\n'.format(length=seq_len , rep_seq_desc=self.rep_seq_id))
            
            count = 0
            # Write rest of cluster members
            for idx, desc in enumerate(self.members_id):
                
                count += 1

                # Calculate percentage similarity
                mismatches = distance(self.members_seq[idx][12:], self.rep_seq[12:])
                percentage = 100.00 - mismatch2percentage(mismatches, seq_len)
                
                clust_file.write('{count}\t{length}nt, >{seq_desc}... at +/{percentage:.2f}%\n'.format(
                    count=str(count), length=seq_len, seq_desc=desc, percentage=percentage 
                                                                                                  ))
                
class ClusterFilter(object):
    """ Holds all methods for parsing and filtering a cluster file from CD-HIT
    
    Impliments a double pass generator to parse reads:
    1. Get the next cluster size, and store start and end places. 

    2. If size passes, go back and rescan file to get detailed info for that cluster.
    """
    
    def __init__(self, handle, idx_file_path, db, filter_params,
                 reads_per_cluster_size_counter, output_dirpath=None):

        # IO
        handle = input_check(handle)
        self.handle4parser = handle
        self.handle4get_desc = open(handle.name, 'rb')
        
        if output_dirpath is None:
            self.output_dirpath = os.getcwd()
        else:
            if not os.path.exists(output_dirpath):
                os.makedirs(output_dirpath)
            self.output_dirpath = output_dirpath
            
        # Store internal params        
        self.idx_file_path = idx_file_path
        self.idx_file_dir = os.path.split(idx_file_path)[0]
        self.db = db
        self.filter_params = filter_params
        self.reads_per_cluster_size_counter = reads_per_cluster_size_counter
        # Setup generator
        self.cluster_size_gen = self.__setup_cluster_size_parser()
        
    def __setup_cluster_size_parser(self):
        
        first = True
        # Setup Data structure
        cluster = ClusterObj()        
        # Cycle through file 
        seq_count = 0    
        try:
            with self.handle4parser as cluster_file:  
                
                while True:
                    line = cluster_file.readline()  # Reading one line from file including end-of-line char '\n'
                    if not line: 
                        break
                
                    if line.startswith('>'):
                        # This is start of new cluster
                        # yield last cluster if not first pass
                        if not first:
                            # Record Cluster info 
                            
                            start_of_line_position = cluster_file.tell() - len(line)
                            cluster.size = seq_count
                            cluster.end_loc = start_of_line_position
                            
                            yield cluster
                            
                            # Reset Cluster and counter
                            cluster = ClusterObj()
                            seq_count = 0
                            cluster.start_loc = start_of_line_position
                        else: 
                            first = False
                    else:
                        seq_count += 1
                            
                # Got to end of file but still one more cluster to add
                cluster.size = seq_count
                yield cluster
        
        except IOError:
            print "Error: can\'t find file or read data"
        else:
            print "Finished Scanning cluster file."    

    def get_desc(self, cluster):
        """ Scan back over cluster and get sequence descriptions """
        
        cluster_file = self.handle4get_desc

        # Update info 
        cluster.idx_file_path = self.idx_file_path

        # Load cluster 
        cluster_file.seek(cluster.start_loc)
        
        first = True
        for line in cluster_file:              
            line = line.strip() # Remove white space at start and end
            
            if line.startswith('>'):
                # This is start of new cluster
                
                # yield this cluster
                if first:
                    cluster.id = line.split()[1]
                    first = False
                else: 
                    break
                
            elif line.endswith('*'): 
                # This is the representative sequence for the cluster
                cluster.rep_seq_id = line.split()[2].strip('>.')
            else: 
                line_parts = line.split()
                next_desc = line_parts[2].strip('>.')
                cluster.members_id.append(next_desc)
                
        if os.getcwd() != self.output_dirpath:
            os.chdir(self.output_dirpath)
                
        return cluster    

    def run_filter(self):
        """ Run through file and write desired clusters to fastq file """
 
        size_counter = Counter()
 
        count = 0
 
        for cluster in self.cluster_size_gen:
            
            count += 1 
            
            if cluster.size >= self.filter_params['min_size'] and cluster.size <= self.filter_params['max_size']:
                
                # Get seq descriptions in cluster 
                cluster = self.get_desc(cluster)
                
                # Write to output file
                size_counter[str(cluster.size)] += 1
                
                # Get the sequence records for the cluster 
                seqs = []
                if os.getcwd() != self.idx_file_dir:
                    os.chdir(self.idx_file_dir)
                seqs.append(self.db[cluster.rep_seq_id])
                for member in cluster.members_id:
                    seqs.append(self.db[member])
                
                if os.getcwd() != self.output_dirpath:
                    os.chdir(self.output_dirpath)
                    
                # Write cluster to a file 
                fname = "clustersize{0}-No{1}.fastq".format(str(cluster.size), str(size_counter[str(cluster.size)]))
                output_handle = open(fname, "wb")
                SeqIO.write(seqs, output_handle, "fastq")
                
            elif self.reads_per_cluster_size_counter[str(cluster.size)] < self.filter_params['min_reads']:
                break
        
        print count


#==============================================================================
# Utility functions
#==============================================================================

def input_check(handle):
    """ Check if input is valid and or a string for a file name.
    Always returns a ahndle for an an open file.
    
    handle = input_check(handle)
    
    """
    if type(handle) is str:
        handle = open(handle, 'rb')
    else:
        assert type(handle) is file, 'Invalid file handle.'
        if handle.closed:
            handle = open(handle.name, 'rb')

    return handle 

def mismatch2percentage(mismatches, seq_length):
    ''' Convert a number of mismatches to a percentage rounded to the nearest 2 dp '''
    mismatches = int(mismatches)
    seq_length = int(seq_length)
    return  round((float(mismatches) / seq_length) * 10000) / 100 

def percentage2mismatch(percentage, seq_length):
    ''' Convert a percentage similarity to a number of mismatches (rounded to integer)'''
    percentage = float(percentage)
    seq_length = int(seq_length)
    return  int(round( (percentage / 100.) * seq_length))
                   
def parse(handle, db=None, edit_dist=False, similarity_count=False):
    """ Reads in a CDHIT cluster file and returns a generator for the clusters. 
    
     - handle   - handle to the file, or the filename as a string
     - db       - a Reads_db database object. If present 
     - edit_dist - Only calculated if set to true.
     - similarity_count - Will return a counter for similarity to rep seq 

    Currently iterates in the order the file is read. 
    Typical usage, opening a file to read in, and looping over all record(s):
    """
    
    # input checks    
    handle = input_check(handle)
        
    first = True
    # Setup Data structure
    cluster = ClusterObj()
    
    if similarity_count:
        cluster.similarity_counter = Counter()
    
    # Cycle through file     
    try:
        with handle as cluster_file:   
            
            for line in cluster_file:              
                line = line.strip() # Remove white space at start and end
                
                if line.startswith('>'):
                    # This is start of new cluster
                    
                    # Claculate stats for last cluster
                    cluster.size = len(cluster.members_id) + 1
                    
                    # Store number for next yield 
                    next_cluster_num = line.split()[1]
                    
                    # yeild this cluster
                    if first:
                        cluster.id = line.split()[1]
                        first = False
                    else: 
                        yield cluster
                        # Reset Cluster
                        cluster = ClusterObj()
                        cluster.id = next_cluster_num
                        if similarity_count:
                            cluster.similarity_counter = Counter()
                    
                elif line.endswith('*'): 
                    # This is the representative sequence for the cluster
                    cluster.rep_seq_id = line.split()[2].strip('>.')
                else:
                    
                    line_parts = line.split()
                    
                    cluster.members_id.append(line_parts[2].strip('>.'))
                    if edit_dist:
                        similarity = line_parts[4].strip('+/%')
                        seq_len = line_parts[1].strip('nt,')
                        cluster.edit_dists.append(percentage2mismatch( 100 - float(similarity), seq_len))
                    if similarity_count:
                        similarity = line_parts[4].strip('+/%')
                        cluster.similarity_counter[similarity] += 1
                    
            # Got to end of file but still one more cluster to add
            cluster.size = len(cluster.members_id) + 1
            yield cluster
    
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        print "Finished Scanning cluster file."    

    
def sortby(handle, reverse=True, mode='cluster_size', outfile_postfix=None, 
           cutoff=None, clustersize_min=None, clustersize_max=None):
    """ Reads in a CDHIT cluster file and writes a sorted file for the clusters.
    by their size.  
    
     - handle   - handle to the file, or the filename as a string.
     - reverse  - Python sorts in ascending order by default, if reverse is true
                  the sort will instead be in descending order.  
     - mode     - "cluster_size" = Sort by the cluster sizes
                  "reads_per_cluster" = Sort by total reads in clusters of each size

     - cutoff   - Minimum threshold for sorted attribute. Any clusters with smaller values 
                  are not written to file.

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
    handle = input_check(handle)

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
            if outfile_postfix is None: 
                outfile_postfix = '-sortedby_' + mode
                
            old_filename_parts = cluster_file.name.split('.')
            old_filename_parts[0]  += outfile_postfix
            new_filename = '.'.join(old_filename_parts)
                
            with open(new_filename, 'wb') as sorted_cluster_file :
                
                if cutoff and clustersize_min and clustersize_max:
                    # Filter the clusters that are written 
                    for tup in sorted_cluster_idxs:
                        if key_ops[mode](tup) >= cutoff:
                            if tup[0] >= clustersize_min and tup[0] <= clustersize_max:
                                cluster_file.seek(tup[1])
                                cluster_text = cluster_file.read(tup[2] - tup[1])
                                sorted_cluster_file.write(cluster_text)
                        else:
                            break
                elif cutoff:
                    for tup in sorted_cluster_idxs:
                        # Stop writing after minimum cutoff is reached
                        if key_ops[mode](tup) >= cutoff:
                            cluster_file.seek(tup[1])
                            cluster_text = cluster_file.read(tup[2] - tup[1])
                            sorted_cluster_file.write(cluster_text)
                        else:
                            break
                else:
                    for tup in sorted_cluster_idxs:
                        # Write all clusters to file
                        cluster_file.seek(tup[1])
                        cluster_text = cluster_file.read(tup[2] - tup[1])
                        sorted_cluster_file.write(cluster_text)
                        
    except IOError:
        print "Error: can\'t find file or read data"
    else:
        t = time.time() - start_time
        print "Finished sorting cluster file after {0}\n".format(time.strftime('%H:%M:%S', time.gmtime(t)))    

    return sorted_cluster_file

    
def getfromdb(cluster_list, items, db=None):
    """ Adds the specified items to each cluster in list, looked up from the main db.
    
    items - list of things to fetch for the cluster. Option of either/both 'seq' and 'phred'
            'rep' will fetch just the seq and phred info for the representative sequence.
            
    """
    for cluster in cluster_list:
        cluster.getfromdb(items, db)
        
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
            
    assert mode in ['by_seqlen', 'cluster_size', 'reads_per_cluster'], 'Mode not recognised.'
            
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

