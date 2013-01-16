'''
Created on 6 Dec 2012

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt
import subprocess as sub
import shlex

def cluster_cdhit(infile, outfile, c_thresh, n_filter, maskN=True):
    ''' Run CD-HIT in parallel on one large fasta file'''
    
    cd_hit_path = '~/bin/cd-hit-v4.6.1/'
    
    if maskN:
        cmd = ('cd-hit-est -i {0} -o {1} -c {2} -n {3}'
              ' -mask N').format(infile, outfile, c_thresh, n_filter)
    else: 
        cmd = ('cd-hit-est -i {0} -o {1} -c {2} -n {3}'
                '').format(infile, outfile, c_thresh, n_filter)

    sub.call(shlex.split(cd_hit_path + cmd))

def cluster_cdhit_para(infile, outfile, c_thresh, n_filter, maskN=True):
    ''' Run CD-HIT in parallel on one large fasta file'''
    
    cd_hit_path = '~/bin/cd-hit-v4.6.1/'
    
    if maskN:
        cmd = ('cd-hit-para.pl -i {0} -o {1} -c {2} -n {3}'
              ' --L 1 --S 64 --P "cd-hit-est -mask N"').format(infile, outfile, 
                                                  c_thresh, n_filter)
    else: 
        cmd = ('cd-hit-para.pl -i {0} -o {1} -c {2} -n {3}'
           ' --L 1 --S 64 --P "cd-hit-est"').format(infile, outfile, 
                                                  c_thresh, n_filter)

    sub.call(shlex.split(cd_hit_path + cmd))

def summary(infile, inpath=None, cluster_sizes=None, seq_lengths=None):
    ''' Display summary of cluster sizes '''
     
    home = os.path.expanduser("~")   
    cd_hit_path = os.path.join(home,'/bin/cd-hit-v4.6.1/')

    if inpath is None:
        inpath = os.getcwd()

    if not infile.endswith('clstr'):
        infile = infile + '.clstr'
    if cluster_sizes is None:
        cluster_sizes = '1,2-4,5-9,10-19,20-49,50-99,100-299,500-99999'
    if seq_lengths is None:
        seq_lengths = '1-88,89,90-150'
       
    cmd =  '{0} {1} {2} {3}'.format(os.path.join(cd_hit_path, 'plot_len1.pl'),
                os.path.join(inpath, infile),
                cluster_sizes, seq_lengths)
  
    print shlex.split(cmd)
    process = sub.Popen(shlex.split(cmd), stdout=sub.PIPE)
    
    output = process.communicate()[0]
    print output 
    
    return output

    


if __name__ == '__main__':
    pass