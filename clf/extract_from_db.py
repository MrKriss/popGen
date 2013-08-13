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

from database.reads_db import Reads_db
  
if __name__ == '__main__':
    
    # Parse arguments

    toc = time.time()

    parser = argparse.ArgumentParser(description='Run Simple clustering on reads in the database.')
    parser.add_argument('-i',  dest='input', required=True,
                        help='Database file where reads are stored (/path/filename)')
    parser.add_argument('-o',  dest='output', 
                        default=sys.stdout,
                        help='Filename for output clusters (/path/filename). Default is to write to stdout.')
    parser.add_argument('-p',  dest='pattern', 
                        default='*',
                        help='Pattern to match to description for fetching records. Default will cycle through all records in database.')
    parser.add_argument('-t',  dest='typeflag', action= 'store_true', 
                        help='Set to match pattern against type column rather than full description column.')
    parser.add_argument('-f',  dest='format', 
                        default='fasta',
                        help='Format of file written to output.')
    
    print sys.argv
    args = parser.parse_args()
    
    # Write records to output
    db = Reads_db(args.input, recbyname=True)
    
    fastafile_handle = db.write_reads(args.pattern, args.output, 
                                      use_type_column=args.typeflag, format='fasta')
    
    