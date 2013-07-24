'''
Created on 6 Dec 2012

@author: musselle
'''
import os 
import sys 
import time
import shlex
import subprocess
import re

import numpy as np 
from collections import Counter, defaultdict


class ClusterClass(object):
    ''' Class to act as a holder of all wrappers for all clustering methods 
    '''
    def __init__(self, infiles=None, inpath=None, config=None, db=None, defaults=None):

        if infiles is not None:
            self.input_files = infiles
        if inpath is not None:
            self.inpath = inpath
        self.clustered_postfix = '-clustered'

        if db:
            self.db = db
        if config:
            self.c = config
             
        # Default Vars for clustering
        self.default_parameters = { 'c_thresh' : 0.90,
                                    'n_filter' : 8,
                                    'threads' : 1,
                                    'mem' : 0,
                                    'maskN' : False}
        if defaults:
            self.default_parameters.update(defaults)
        
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
        out = self.cluster_cdhit(**inputs_dict)
        
        return out

    def run_batch_cdhit_clustering(self, batch_parameters, **kwargs):
        ''' Runs cd-hit-est over the list of files given, for each setup listed in
        batch_parameters.
        
        Elements of batch_parameters are dictionaries containing the parameters
        to be changed from the default value. 
        
        Output is a list of length (batch parameters) containing tuples of 
        
        (returned_outfiles_list, outpath, counters) for each set of parameters
        
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
            
            out = self.cluster_cdhit(**inputs_dict)
            outputs_list.append(out)
            
        return outputs_list

    def cluster_cdhit(self, infiles, inpath='', outpath='', outfile_postfix=None, 
                      c_thresh=None, n_filter=None, threads=1, mem=0, maskN=True, 
                      allvall = False):
        ''' Run CD-HIT in parallel on list of fasta files. Each file is clustered seperately.
        
        Other flags used:
        -d 0   --> No limit on description written to cluster file (goes to first space in seq ID). 
        -r 0   --> DO Only +/+ and not -/+ alignment comparisons as reads are done in both directions but on different strands. 
        -s 0.8 --> If shorter sequence is less than 80% of the representative sequence, dont cluster. 
        
        Writes stdout to console.
        Counter dictionary and summary logfile are generated after each run. 
        
        '''
        if outfile_postfix is None:
            outfile_postfix = self.clustered_postfix
        
        # input check
        if type(infiles) is str: 
            infiles = [infiles]
        
        returned_outfiles_list = []
        counters = []
        for infile in infiles:
            
            start_time = time.time()
            
            out_filename = infile.split('.')[0] + outfile_postfix
            returned_outfiles_list.append(out_filename)
            
            log_filename = infile.split('.')[0] + '-report.log'
            
            infile_path = os.path.join(inpath, infile)
            outfile_path = os.path.join(outpath, out_filename)
            logfile_path = os.path.join(outpath, log_filename)
        
            cmd = ('cd-hit-est -i {0} -o {1} -c {2} -n {3} -d 0 -r 0 -s 0.8 -M {4} '
                '-T {5}').format(infile_path, outfile_path, c_thresh, n_filter, mem, threads)   
        
            if maskN:
                cmd = cmd + ' -mask N'
            if allvall:
                cmd = cmd + ' -g 1'
            
            # Process to run CD-HIT
            subprocess.check_call(shlex.split(os.path.join(self.c.cdhit_path, cmd)))
            
            finish_time = time.time()
            
            # Update database if present 
            if self.db:
                
                # Get cluster size summary counter 
                total_counter, by_seqlen_counter = self.cluster_summary_counter(infile=out_filename, path=outpath,
                                 mode='both', report=True)    
            
                st_idx = cmd.find('-c ')
                CDHIT_parameters = cmd[st_idx:]
            
                # Write summary logfile 
                with open(logfile_path, 'wb') as f:
                    program_name = os.path.join(self.c.cdhit_path, cmd).split(' -i ')[0]
                    f.write('=========================================================\n')
                    f.write('Program     : {0}\n'.format(program_name))
                    f.write('Input File  : {0}\n'.format(infile_path))
                    f.write('Output File : {0}\n'.format(outfile_path))
                    f.write('Commands    : {0}\n'.format(CDHIT_parameters))
                    f.write('\n')
                    f.write('Started     : {0}\n'.format(time.strftime('%a, %d %b %Y, %H:%M:%S', 
                                                            time.gmtime(start_time))))
                    f.write('=========================================================\n')
                    f.write('\n')
                    f.write('                       Report Log\n')
                    f.write('---------------------------------------------------------\n')
                    
                    reads_per_cluster = {key: int(key)*value for key, value in total_counter.iteritems()}
                    total_reads = sum(reads_per_cluster.values())
                    total_clusters = sum(total_counter.values())
                    f.write('Total number of reads     : {0}\n'.format(total_reads))
                    f.write('Total number of clusters  : {0}\n'.format(total_clusters))
                    read_lengths = [int(key) for key in by_seqlen_counter.keys()]
                    f.write('Read length Min and Max    : {0} and {1}\n'.format(min(read_lengths), max(read_lengths)))
                    f.write('Time taken                 : {0}\n'.format(time.strftime('%H:%M:%S', 
                                                            time.gmtime(finish_time - start_time))))
                    f.write('\n')
                    f.write('Top 20 Percentage Reads per cluster \n')
                    f.write('---------------------------------------------------------\n')
                    f.write('Cluster Size    No. Clusters    Total Reads         %    \n')
                    f.write('---------------------------------------------------------\n')
                    top_reads_per_cluster = sorted(reads_per_cluster.iteritems(), 
                                                   key=lambda tup: int(tup[1]), reverse=True)[:20]
                    for tup in top_reads_per_cluster:
                        if total_reads == 0:
                            perc = 0.0
                        else:
                            perc = float(tup[1]) / total_reads
                        
                        f.write("{clust_size: <16}{num_clust: <16}{total_reads: <18d}{percentage:.2%}\n".format(
                              clust_size=tup[0], num_clust=total_counter[tup[0]], total_reads=tup[1], 
                              percentage=perc))

                self.db.add_results_parameters_datafiles(infile, out_filename, total_counter, self.c, CDHIT_parameters)
                
                counters.append(total_counter)
            
        return (returned_outfiles_list, outpath, counters) 


    def cluster_summary_counter(self, infile, path='', mode='total', report=True):
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
        
    
    def hist_counters(self, counters, labels=None, **kwargs):
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