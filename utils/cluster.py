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
    
    if maskN:
        cmd = ('cd-hit-para.pl -i {0} -o {1} -c {2} -n {3}'
              ' --L 1 --S 64 --P "cd-hit-est -mask N"').format(infile, outfile, 
                                                  c_thresh, n_filter)
    else: 
        cmd = ('cd-hit-para.pl -i {0} -o {1} -c {2} -n {3}'
           ' --L 1 --S 64 --P "cd-hit-est"').format(infile, outfile, 
                                                  c_thresh, n_filter)

    sub.call(shlex.split(cmd))

def summary(infile, cluster_sizes=None, seq_lengths=None):
    ''' Display summary of cluster sizes '''

    if not infile.endswith('clstr'):
        infile = infile + '.clstr'
    if cluster_sizes is None:
        cluster_sizes = '1 ,2-4 ,5-9 ,10-19 ,20-49 ,50-99 ,100-299 ,500-99999'
    if seq_lengths is None:
        seq_lengths = '1-78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90'
    
    cmd = 'plot_len.pl {0} {1} {2}'.format(infile, cluster_sizes, seq_lengths)
    
    process = sub.Popen(shlex.split(cmd), stdout=sub.PIPE)
    
    output = process.communicate()[0]
    print output 
    
    return output

    


if __name__ == '__main__':
    pass