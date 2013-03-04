'''
Created on 6 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
from subprocess import call, Popen, STDOUT, PIPE
import shlex


class Clustering(object):
    ''' Class to act as a holder of all wrappers for all clustering methods 
    '''
    
    def __init__(self, preprocessor):

        self.c = preprocessor.c
        self.next_inputs = preprocessor.next_input_files
        self.next_path = preprocessor.next_input_path
        
        # default Vars for clustering 
        self.default_parameters = { 'c_thresh' : 0.90,
                                    'n_filter' : 8,
                                    'threads' : 1,
                                    'mem' : 0,
                                    'maskN' : False}
        
        
        

    def run_cdhit_clustering(self, **kwargs): 

        if 'infile' not in kwargs:
            kwargs['infile'] = os.path.join(self.next_input_path, self.next_input_files)
        if 'outfile' not in kwargs:
            kwargs['outfile'] = 'all_reads'
            
        cluster_cdhit(**kwargs)





# default Vars for clustering 
default_vars = { 'c_thresh' : 0.90,
                 'n_filter' : 8,
                 'threads' : 1,
                 'mem' : 0,
                 'maskN' : False}

# Variations to run
clustering_runs = [ { 'c_thresh' : 0.95},
                    { 'c_thresh' : 0.95, 'maskN' : True},
                    { 'c_thresh' : 0.90},
                    { 'c_thresh' : 0.90, 'maskN' : True},
                    { 'c_thresh' : 0.85},
                    { 'c_thresh' : 0.85, 'maskN' : True},
                   ]
                   
for d in clustering_runs:
    
    inputs_dict = {}
    inputs_dict.update(default_vars)
    inputs_dict.update(d)
    
    dirname = experiment_name + '_clustered_reads'
    outfile = experiment_name + '_clustered_reads'
    if 'c_thresh' in d:
        dirname = dirname + '-c{}'.format(int(d['c_thresh']*100))
        outfile = outfile + '-c{}'.format(int(d['c_thresh']*100))
    if 'maskN' in d:
        dirname = dirname + '-maskN'
        outfile = outfile + '-maskN'
    
    path = os.path.join(c.clusters_outpath, dirname)        
    if not os.path.exists(path):
        os.makedirs(path)
        
    path2outfile  = os.path.join(path, outfile)
    inputs_dict['log_filename'] = os.path.join(path, 'report.log')

    Experiment.run_cdhit_clustering(infile=allreads_file, outfile=path2outfile,
              **inputs_dict)

























def cluster_cdhit(infile, outfile, c_thresh, n_filter, threads=1, 
                  mem=0, maskN=True, log_filename='cd-hit-report.log'):
    ''' Run CD-HIT in parallel on one large fasta file
    
    Other flags used:
    -d 0   --> No limit on description written to cluster file (goes to first space in seq ID). 
    -r 0   --> Do only +/+ alignment comparisons
    -s 0.8 --> If shorter sequence is less than 80% of the representative sequence, dont cluster. 
    
    Writes stdout to console and saves to log file in real time. 
    
    '''
    
    cd_hit_path = os.path.expanduser("~/bin/cd-hit-v4.6.1/")
    
    cmd = ('cd-hit-est -i {0} -o {1} -c {2} -n {3} -d 0 -r 0 -s 0.8 -M {4} '
            '-T {5}').format(infile, outfile, c_thresh, n_filter, mem, threads)   
    if maskN:
        cmd = cmd + ' -mask N'
  
    proc = Popen(shlex.split(os.path.join(cd_hit_path, cmd)), stdout=PIPE)
    
    with open(log_filename, 'wb') as logfile:
        while True:
            out = proc.stdout.readline()
            if out == '' and proc.poll() != None:
                break
            if out != '':
                sys.stdout.write(out)
                sys.stdout.flush()
                logfile.write(out)
    

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

    subprocess.call(shlex.split(cd_hit_path + cmd))

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
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    
    output = process.communicate()[0]
    print output 
    
    return output

    
if __name__ == '__main__':
    pass