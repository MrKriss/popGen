'''
Created on 6 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
from subprocess import call, Popen, STDOUT, PIPE
import shlex
import re
import string
from threading  import Thread
from Queue import Queue, Empty

from collections import Counter, defaultdict

class Clustering(object):
    ''' Class to act as a holder of all wrappers for all clustering methods 
    '''
    def __init__(self, config, db, infiles, inpath):

        self.c = config
        self.input_files = infiles
        self.inpath = inpath
        self.db = db
        self.clusterd_postfix = 
             
        # Default Vars for clustering 
        self.default_parameters = { 'c_thresh' : 0.90,
                                    'n_filter' : 8,
                                    'threads' : 1,
                                    'mem' : 0,
                                    'maskN' : False}
        
        
    def run_single_cdhit_clustering(self, **kwargs):
        ''' Runs cd-hit-est over the list of files given for a single set of parameters. ''' 

        inputs_dict = {}
        inputs_dict.update(self.default_parameters)
        inputs_dict.update(kwargs)

        # Use defaults if no others were passed 
        if 'infiles' not in inputs_dict:
            inputs_dict['infile'] = self.input_files
            
        if 'inpath' not in inputs_dict:
            inputs_dict['inpath'] = self.inpath
            
        # (outfiles, outfiles.clstr, outpath, cmd) 
        #    = cluster_cdhit(infiles, inpath=None, outpath=None, 
        #                  outfile_postfix='-clustered', c_thresh=None, 
        #                  n_filter=None, threads=1, mem=0, maskN=True, 
        #                  allvall = False)
        out = cluster_cdhit(**inputs_dict)
        
        # Get cluster size summary counter 
        counter = cluster_summary_counter(infile=out[0], mode='total', report=True)
        
        # Update database
        self.db.add_results .....
        
            db.
        
        
        
        
        return out

    def run_batch_cdhit_clustering(self, batch_parameters, **kwargs):
        ''' Runs cd-hit-est over the list of files given, for each setup listed in
        batch_parameters.
        
        Elements of batch_parameters are dictionaries containing the parameters
        to be changed from the default value. 
        
        Output is a list of length (batch parameters) containing tuples of 
        
        (outfiles1, outfiles2, outpath, cmd) for each set of parameters
        
        '''
        outputs_list = [] 

        for d in batch_parameters:
    
            inputs_dict = {}
            inputs_dict.update(self.default_parameters)
            inputs_dict.update(d)
            inputs_dict.update(kwargs)
    
            # Use defaults if no others were passed 
            if 'infiles' not in inputs_dict:
                inputs_dict['infiles'] = self.input_files
                
            if 'inpath' not in inputs_dict:
                inputs_dict['inpath'] = self.inpath
        
            dirname = self.c.experiment_name + '_clustered_reads'
            outfile_postfix = '-clustered'

            if 'c_thresh' in d:
                dirname = dirname + '_c{}'.format(int(d['c_thresh']*100))
                outfile_postfix = outfile_postfix + '_c{}'.format(int(d['c_thresh']*100))
            if 'n_filter' in d:
                dirname = dirname + '_n{}'.format(d['n_filter'])
                outfile_postfix = outfile_postfix + '_n{}'.format(d['n_filter'])                
            if 'maskN' in d:
                dirname = dirname + '_maskN'
                outfile_postfix = outfile_postfix + '-maskN'
            if 'allvall' in d:
                dirname = dirname + '_g1'
                outfile_postfix = outfile_postfix + '_g1'
                
            inputs_dict['outfile_postfix'] = outfile_postfix
            
            path = os.path.join(self.c.clusters_outpath, dirname)        
            if not os.path.exists(path):
                os.makedirs(path)
            
            inputs_dict['outpath'] = path
            
            out = cluster_cdhit(**inputs_dict)
            outputs_list.append(out)
            
        return outputs_list

def cluster_cdhit(infiles, inpath='', outpath='', outfile_postfix='-clustered', 
                  c_thresh=None, n_filter=None, threads=1, mem=0, maskN=True, 
                  allvall = False):
    ''' Run CD-HIT in parallel on list of fasta files. Each file is clustered seperately.
    
    Other flags used:
    -d 0   --> No limit on description written to cluster file (goes to first space in seq ID). 
    -r 0   --> DO Only +/+ and not -/+ alignment comparisons as reads are done in both directions but on different strands. 
    -s 0.8 --> If shorter sequence is less than 80% of the representative sequence, dont cluster. 
    
    Writes stdout to console and saves to log file in real time. 
    
    '''
    # input check
    if type(infiles) is not list: 
        infiles = [infiles]
          
    # Define re patterns 
    pattern1 = re.compile('^\.')
    pattern2 = re.compile('%$')
    pattern3 = re.compile('^#')
    
    cd_hit_path = os.path.expanduser("~/bin/cd-hit-v4.6.1/")
    
    returned_outfiles_list1 = []
    for f in infiles:
        
        out_filename = f.split('.')[0] + outfile_postfix
        
        returned_outfiles_list1.append(out_filename)
        log_filename = f.split('.')[0] + '-report.log'
        
        infile_path = os.path.join(inpath, f)
        outfile_path = os.path.join(outpath, out_filename)
        logfile_path = os.path.join(outpath, log_filename)
    
        cmd = ('cd-hit-est -i {0} -o {1} -c {2} -n {3} -d 0 -r 0 -s 0.8 -M {4} '
            '-T {5}').format(infile_path, outfile_path, c_thresh, n_filter, mem, threads)   
    
        if maskN:
            cmd = cmd + ' -mask N'
        if allvall:
            cmd = cmd + ' -g 1'
            
        # Wish to write progress to console and also capture summary in a log file
        
        # function to thread
        def enqueue_output(out, queue):
            for line in iter(out.readline, b''):
                queue.put(line)
            out.close()
        
        # Process to run CD-HIT
        proc = Popen(shlex.split(os.path.join(cd_hit_path, cmd)), stdout=PIPE, bufsize=1)
        
        # Setup queue and threading 
        q = Queue()
        t = Thread(target=enqueue_output, args=(proc.stdout, q))
        t.daemon = True # thread dies with the program
        t.start()
    
        with open(logfile_path, 'wb') as logfile:
            while True:
                # read line without blocking
                try:  line = q.get_nowait() # or q.get(timeout=.1)
                except Empty:
                    if proc.poll() != None:
                        break
                else: # got line
                    if (pattern1.match(line) or pattern2.match(line) 
                      or pattern3.match(line)):
                        sys.stdout.write(line)
                        sys.stdout.flush()
                        continue
                    else:
                        line.translate(string.maketrans("\r","\n"))
                        sys.stdout.write(line)
                        sys.stdout.flush()
                        logfile.write(line)
                        logfile.flush()
       
    returned_outfiles_list2 = [0]* len(returned_outfiles_list1)             
    for i, f in enumerate(returned_outfiles_list1):
        returned_outfiles_list2[i] = f + '.clstr'         
    
    return (returned_outfiles_list1, returned_outfiles_list2, outpath, cmd) 


def cluster_summary_counter(infile, mode='by_seqlen', report=True):
    ''' Takes cluster file output by CD-Hit and produces two Counter for the 
    
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

    # Helper vars
    cluster_size_counter = Counter()
    seq_length_counter = Counter()

    seq_in_cluster = 0
    rep_length = 0

    print 'Generating cluster summary for  %s ...' % (infile)

    try:
        with open(infile, 'rb')  as cluster_file:   
            
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
    for k in ds.iterkeys():
        seq_len_counter.update(int(k)) 
    
    if report:
        print 'Top 5 Cluster Sizes: ', total_cluster_size_counter.most_common()[:5]
        print 'Top 5 Sequence Lengths: ', seq_len_counter.most_common()[:5]
    
    # Decide what to output    
    if mode == 'total':
        return total_cluster_size_counter
    elif mode == 'by_seqlen'
        return ds
    
def hist_counter(counter, **kwargs):
    ''' Construct a histogram from a Counter Dictionary '''
    
    data = np.array(list(counter.elements()), dtype = np.int)

    plt.hist(data, histtype='step', **kwargs)
    plt.title("Cluster Size Distribution")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
#    plt.legend()
    plt.show()
    

def plot_counter( bins,cluster_length_bins=None, mincutoff=10, bins=5000, report=True, plot_hist=True,):
    
    
    # Process ds to extract needed data for histergrams
    if cluster_length_bins is not None:
        # Update counters for each bin in cluster_length_bins
        for cluster_size, seq_len_c in ds.iteritems():

            for tup in cluster_length_bins:   
                name = 'Length-' + str(tup)
                 
                if len(tup) == 1:
                    # single number bin
                    target_key = str(tup[0])
                    if seq_len_c.has_key(target_key):
                        vars()[name][cluster_size] += seq_len_c[target_key]
        
                elif len(tup) == 2:
                    # range bin
                    for seq_length in seq_len_c.iterkeys():
                        if int(seq_length) >= tup[0] and int(seq_length) <= tup[1]:
                            vars()[name][cluster_size] += seq_len_c[seq_length]
    
    #===============================================================================
    # Plot results
    #===============================================================================
    
    # TODO Add a level of input checking so if Counter is empty, an empty histogram is plotted
    
    if plot_hist:
    
        plt.figure()
        if cluster_length_bins is not None:   
            for tup in cluster_length_bins:   
                name = 'Length-' + str(tup)
                
                hist_counter(vars()[name], bins=bins, label=name, range=(mincutoff, 10000))
         
        hist_counter(cluster_size_counter, bins=bins, label='All Seq Lengths', range=(mincutoff, 10000))   
        plt.title("Cluster Size Distribution")
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.legend()
        plt.show()
    
    
    
# Deprecated
def cluster_cdhit_para(infile, outfile, c_thresh, n_filter, maskN=True):
    ''' Run CD-HIT in parallel on one large fasta file'''
    
    cd_hit_path = os.path.expanduser("~/bin/cd-hit-v4.6.1/")
    
    if maskN:
        cmd = ('cd-hit-para.pl -i {0} -o {1} -c {2} -n {3}'
              ' --L 1 --S 64 --P "cd-hit-est -mask N"').format(infile, outfile, 
                                                  c_thresh, n_filter)
    else: 
        cmd = ('cd-hit-para.pl -i {0} -o {1} -c {2} -n {3}'
           ' --L 1 --S 64 --P "cd-hit-est"').format(infile, outfile, 
     
                                                  c_thresh, n_filter)

    call(shlex.split(cd_hit_path + cmd))

# Depricated
def summary(infile, data_inpath=None, cluster_sizes=None, seq_lengths=None):
    ''' Display summary of cluster sizes '''
     
    cd_hit_path = os.path.expanduser("~/bin/cd-hit-v4.6.1/")   
     
    if data_inpath is None:
        data_inpath = os.getcwd()

    if not infile.endswith('clstr'):
        infile = infile + '.clstr'
    if cluster_sizes is None:
        cluster_sizes = '1,2-4,5-9,10-19,20-49,50-99,100-299,500-99999'
        
    if cluster_sizes == 'all':
        # All clusters from size 1 to 10,000
        cluster_sizes = ','.join([str(i) for i in xrange(1,10001)])

        
    if seq_lengths is None:
        seq_lengths = '1-88,89,90-150'
       
    cmd =  '{0} {1} {2} {3}'.format(os.path.join(cd_hit_path, 'plot_len1.pl'),
                os.path.join(data_inpath, infile),
                cluster_sizes, seq_lengths)

    cmd = shlex.split(cmd)  
    print 'Running\n' + cmd[0] 
    process = Popen(cmd, stdout=PIPE)
    
    output = process.communicate()[0]
    print output 
    
    return output

if __name__ == '__main__':
    pass