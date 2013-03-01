'''
Created on 6 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
from subprocess import call, Popen, STDOUT, PIPE
import shlex

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