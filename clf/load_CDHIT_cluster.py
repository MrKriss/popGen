#! /usr/bin/env python 

'''
Created on 17 Jul 2013

@author: musselle
'''

import _addpaths
import os, sys, time, gzip, argparse
import shlex, subprocess
  
import numpy as np 

from editdist import distance
from Bio import SeqIO, bgzf 
from collections import Counter, defaultdict

from utils.fileIO import inputfile_check, outputfile_check
from database.reads_db import Reads_db


if __name__ == '__main__':
    
    # Parse arguments

    toc = time.time()

    parser = argparse.ArgumentParser(description='Run Simple clustering on reads in the database.')
    
    # IO parameters
    parser.add_argument('-i',  dest='input', 
                        required=True,
                        help='CDHIT cluster file to load (/path/filename)')
    parser.add_argument('-o',  dest='output', default='clusterfile',
                        help='Database with which to update/load cluster file. (/path/filename).')
    
    parser.add_argument('-p',  dest='tableprefix', default=None,
                        help='Database with which to update/load cluster file. (/path/filename).')
    
    parser.add_argument('-b',  dest='buffer', default=1000000,
                        help='Number of sequences to read (across multiple clusters) before writing to db.')
    
    parser.add_argument('-n',  dest='overwrite', action = 'store_true',
                        help='Overwrite all previously loaded clusters.')
    
    
    parser.add_argument('--min',  dest='fmin', default = 2, type=int,
                        help='Minimum size of clusters to load. Default = 2 to miss out singleton clusters')
    parser.add_argument('--max',  dest='fmax', default = 0, type=int, 
                        help='Maximum size of clusters to load. Default = 2 to miss out singleton clusters')
    
    parser.add_argument('--skipsort',  dest='skipsort', action = 'store_true',
                        help='Maximum size of clusters to load. Default = 2 to miss out singleton clusters')
    
    print sys.argv
    args = parser.parse_args()
    
    if os.path.exists(args.output):
        
        db = Reads_db(args.output)
    else:
        raise Exception('Database file not found.')
    
    # Load cluster file 
    db.load_cluster_file(args.input, args.tableprefix, args.overwrite, args.fmin, 
                         args.fmax, args.skipsort, buffer_max=args.buffer)
    
    total_t = time.time() - toc    
    print >> sys.stderr, 'Loaded cluster file in {0}'.format(
              time.strftime('%H:%M:%S', time.gmtime(total_t)))
    
    
    
    
    
    
    