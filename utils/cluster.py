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

class Clustering(object):
    ''' Class to act as a holder of all wrappers for all clustering methods 
    '''
    
    def __init__(self, config, file2cluster):

        self.c = config
        self.input_file = file2cluster
             
        # default Vars for clustering 
        self.default_parameters = { 'c_thresh' : 0.90,
                                    'n_filter' : 8,
                                    'threads' : 1,
                                    'mem' : 0,
                                    'maskN' : False}
        
    def run_single_cdhit_clustering(self, **kwargs):
        ''' Runs a single instance of cd-hit-est ''' 

        # Use defaults if no others were passed 
        if 'infile' not in kwargs:
            kwargs['infile'] = self.input_file
        if 'outfile' not in kwargs:
            kwargs['outfile'] = self.c.experiment_name + '_all_reads'
        if 'c_thresh' not in kwargs:
            kwargs['c_thresh'] = self.default_parameters['c_thresh']
        if 'n_filter' not in kwargs:
            kwargs['n_filter'] = self.default_parameters['n_filter']
        if 'threads' not in kwargs:
            kwargs['threads'] = self.default_parameters['threads']
        if 'mem' not in kwargs:
            kwargs['mem'] = self.default_parameters['mem']
        if 'maskN' not in kwargs:
            kwargs['maskN'] = self.default_parameters['maskN']
            
        cluster_cdhit(**kwargs)

    def run_batch_cdhit_clustering(self, batch_parameters, **kwargs):
        ''' Runs cdhit in serial for each entry in batch_parameters.
        
        Elements of batch_parameters are dictionaries containing the parameters
        to be changed from the default values '''

        for d in batch_parameters:
    
            inputs_dict = {}
            inputs_dict.update(self.default_parameters)
            inputs_dict.update(d)
            inputs_dict.update(kwargs)
    
            dirname = self.c.experiment_name + '_clustered_reads'
            outfile = self.c.experiment_name + '_clustered_reads'

            if 'c_thresh' in d:
                dirname = dirname + '-c{}'.format(int(d['c_thresh']*100))
                outfile = outfile + '-c{}'.format(int(d['c_thresh']*100))
            if 'n_filter' in d:
                dirname = dirname + '-n{}'.format(d['n_filter'])
                outfile = outfile + '-n{}'.format(d['n_filter'])                
            if 'maskN' in d:
                dirname = dirname + '-maskN'
                outfile = outfile + '-maskN'
            
            path = os.path.join(self.c.clusters_outpath, dirname)        
            if not os.path.exists(path):
                os.makedirs(path)
                
            inputs_dict['outfile'] = os.path.join(path, outfile)
            inputs_dict['log_filename'] = os.path.join(path, 'report.log')
        
            if 'infile' not in inputs_dict:
                inputs_dict['infile'] = self.input_file
            
            cluster_cdhit(**inputs_dict)

def cluster_cdhit(infile, outfile, c_thresh, n_filter, threads=1, 
                  mem=0, maskN=True, log_filename='cd-hit-report.log'):
    ''' Run CD-HIT in parallel on one large fasta file
    
    Other flags used:
    -d 0   --> No limit on description written to cluster file (goes to first space in seq ID). 
    -r 0   --> Do only +/+ alignment comparisons
    -s 0.8 --> If shorter sequence is less than 80% of the representative sequence, dont cluster. 
    
    Writes stdout to console and saves to log file in real time. 
    
    '''
    
    pattern1 = re.compile('^\.')
    pattern2 = re.compile('%$')
    
    cd_hit_path = os.path.expanduser("~/bin/cd-hit-v4.6.1/")
    
    cmd = ('cd-hit-est -i {0} -o {1} -c {2} -n {3} -d 0 -r 0 -s 0.8 -M {4} '
            '-T {5}').format(infile, outfile, c_thresh, n_filter, mem, threads)   
    if maskN:
        cmd = cmd + ' -mask N'
  
    proc = Popen(shlex.split(os.path.join(cd_hit_path, cmd)), stdout=PIPE)
    
    # Write progress to console and a summary to a log file
    with open(log_filename, 'wb') as logfile:
        while True:
            out = proc.stdout.readline()
            if pattern1.match(out) or pattern2.match(out):
                sys.stdout.write(out)
                sys.stdout.flush()
                continue
            out.translate(string.maketrans("\r","\n"))
            if not out and proc.poll() != None:
                break
            else:
                sys.stdout.write(out)
                sys.stdout.flush()
                logfile.write(out)
                logfile.flush()
    

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