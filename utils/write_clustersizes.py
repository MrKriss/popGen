'''
Created on 6 Jun 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser(description='Procedure to write CDHIT.clstr file to FastQ')



class Cluster(object):
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
            
        self.idx_file_path = ""
        
    def getfromdb(self, items, db=None):
        """ Lookup to indexed db to retrieve and store the specified items.
        
        items - list of things to fetch for the cluster. Option of either/both 'seq' and 'phred'
                'rep' will fetch just the seq and phred info for the representative sequence.
        """
        
        assert 'seq' in items or 'phred' in items or 'rep' in items, "Invalid values for items to lookup"
        
        if not self.idx_file_path:
            raise Exception('No idx_file_path specified.')
        
        if db is None: 
            print "Loading {} ...".format(os.path.split(self.idx_file_path)[1])
            db = SeqIO.index_db(self.idx_file_path)  
            print "Loading complete"
        
        start_dir = os.getcwd()
        path = os.path.split(self.idx_file_path)[0]
        
        if start_dir != path:
            os.chdir(path) 

        if 'seq' in items and 'phred' in items:
            self.rep_seq = db[self.rep_seq_id].seq.tostring()
            self.rep_phred = np.array(db[self.rep_seq_id].letter_annotations['phred_quality'])
            for elem in self.members_id:
                self.members_seq.append(db[elem].seq.tostring())
                self.members_phred.append(np.array(db[elem].letter_annotations['phred_quality']))

        if 'seq' in items:
            self.rep_seq = db[self.rep_seq_id].seq.tostring()
            for elem in self.members_id:
                self.members_seq.append(db[elem].seq.tostring())

        if 'phred' in items:
            self.rep_phred = np.array(db[self.rep_seq_id].letter_annotations['phred_quality'])
            for elem in self.members_id:
                self.members_phred.append(np.array(db[elem].letter_annotations['phred_quality']))

        if 'rep' in items: # Just fetch the representitive sequence info 
            self.rep_seq = db[self.rep_seq_id].seq.tostring()
            self.rep_phred = np.array(db[self.rep_seq_id].letter_annotations['phred_quality'])
        
        if os.getcwd() != start_dir:
            os.chdir(start_dir)
                    
    def get_unique_seq(self, ignoreup2=6, db=None):
        """ Work out the counts for unique reads within a cluster """
        
        if not self.members_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(self, ['seq'], db=db)
        if not self.rep_seq:
            print "Sequence data not present in cluster. Retrieving from data base..."
            self.getfromdb(self, ['rep'], db=db)
        
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
        
    def write2fastq(self):
        """ Write cluster to a fastq file format. First seq is representative seq for the cluster """
        pass



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
        self.lookup_db = db
        self.filter_params = filter_params
        self.reads_per_cluster_size_counter = reads_per_cluster_size_counter
        # Setup generator
        self.cluster_size_gen = self.__setup_cluster_size_parser()
        
    def __setup_cluster_size_parser(self):
        
        first = True
        # Setup Data structure
        cluster = Cluster()        
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
                            cluster = Cluster()
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
                seqs.append(self.lookup_db[cluster.rep_seq_id])
                for member in cluster.members_id:
                    seqs.append(self.lookup_db[member])
                
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


if __name__ == '__main__':
    pass