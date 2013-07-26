#! /space/musselle/EPD7-3/bin/python 

'''
Created on 6 Jun 2013

@author: musselle

Read a CDHIT cluster file from stdin and filter based on parsed parameters

write those clusters that  

'''

import sys
import os
import socket

# Work out where data is stored
if socket.gethostname() == 'yildun':
    data_prefix = '/space/musselle/data'
    work_prefix = '/home/pgrad/musselle/ubuntu/workspace/'
        
elif socket.gethostname() == 'luca':
    data_prefix = '/home/musselle/san/data'
    work_prefix = '/home/musselle/'
       
elif socket.gethostname() == 'gg-pc6':
    data_prefix = '/home/musselle/data'
    work_prefix = '/home/musselle/'

sys.path.append(os.path.join(work_prefix, 'popGen'))
sys.path.append(os.path.join(work_prefix, 'popGen/utils'))
sys.path.append(os.path.join(work_prefix, 'popGen/plotscripts'))
sys.path.append(os.path.join(work_prefix, 'popGen/runscripts'))
sys.path.append(os.path.join(work_prefix, 'popGen/clf'))

del work_prefix
del data_prefix

import argparse

from utils.ClusterIO import sortby

parser = argparse.ArgumentParser(description='Procedure to Filter a CDHIT.clstr file, writing a copy that only contains passes.')
parser.add_argument('filename', action='store', help='Name of CD-HIT output cluster file.')
parser.add_argument('min', action='store', type=int, help='Minimum value of threshold for cluster size')
parser.add_argument('max', action='store', type=int, help='Maximum value of threshold for cluster size')
parser.add_argument('minreads', action='store', type=int, help='Minimum value of threshold for total reads in clusters of any size')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

args = parser.parse_args()

# Sort file 
sortby(args.filename, reverse=True, mode='reads_per_cluster', outfile_postfix='-filtered', cutoff=args.minreads, 
       clustersize_min=args.min, clustersize_max=args.max)

# Rewrite file
old_filename_parts = args.filename.split('.')
old_filename_parts[0]  += '-sortedby_' + 'reads_per_cluster'
new_filename = '.'.join(old_filename_parts)

