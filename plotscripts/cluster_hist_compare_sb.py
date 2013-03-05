'''
Created on 5 Mar 2013

@author: musselle
'''
import os 
import sys 

from plot_utils import cluster_summary_plot, hist_counter

from preprocess import  Preprocessor, ConfigClass
from cluster import Clustering

import socket 

# Work out where data is stored
if socket.gethostname() == 'yildun':
    prefix = '/space/musselle/datasets'
elif socket.gethostname() == 'luca':
    prefix = '/home/musselle/san/data'
elif socket.gethostname() == 'gg-pc6':
    prefix = '/home/musselle/data'

c = ConfigClass()
c.experiment_name = 'sb'
c.clusters_outpath = os.path.join(prefix,'sticklebacks/clusters')

batch_parameters = [ { 'c_thresh' : 0.95},
                    { 'c_thresh' : 0.95, 'maskN' : True},
                    { 'c_thresh' : 0.90},
                    { 'c_thresh' : 0.90, 'maskN' : True},
                    { 'c_thresh' : 0.85},
                    { 'c_thresh' : 0.85, 'maskN' : True},
                   ]

outfiles_list = [] 

for d in batch_parameters:

    dirname = c.experiment_name + '_clustered_reads'
    outfile = c.experiment_name + '_clustered_reads'

    if 'c_thresh' in d:
        dirname = dirname + '-c{}'.format(int(d['c_thresh']*100))
        outfile = outfile + '-c{}'.format(int(d['c_thresh']*100))
    if 'n_filter' in d:
        dirname = dirname + '-n{}'.format(d['n_filter'])
        outfile = outfile + '-n{}'.format(d['n_filter'])                
    if 'maskN' in d:
        dirname = dirname + '-maskN'
        outfile = outfile + '-maskN'
    
    path = os.path.join(c.clusters_outpath, dirname)            
    outfiles_list.append(os.path.join(path, outfile))
    

for cluster_file in outfiles_list:
    
    name = os.path.split(cluster_file)[1]  
    out = cluster_summary_plot(cluster_file, plot_hist = 0)  
    hist_counter(out[1], bins=5000, range =(1,10000),label=name)
    
    

