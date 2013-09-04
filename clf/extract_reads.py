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

    parser = argparse.ArgumentParser(description='''Script to extract selected reads from specified database and write them to a fasta or fastq file.''')
    parser.add_argument('-i',  dest='input', required=True,
                        help='Database file where reads are stored (/path/filename)')
    parser.add_argument('-o',  dest='output', 
                        default=sys.stdout,
                        help='Filename for output clusters (/path/filename). Default is to write to stdout.')
    parser.add_argument('-e',  dest='filter_expression', 
                        default=None,
                        help='''SQL expression to filter the query which selects the sequences in the database. Default is to export all sequences in database.
                         Basic query is:
                             SELECT seqid, seq, phred FROM seqs INNER JOIN samples ON seqs.sampleId=samples.sampleId 
                                    WHERE <filter_expression> ''')
    parser.add_argument('-s',  dest='startidx', 
                        default=0,
                        help='Starting base index of DNA sequences that are written to file, used to miss out cutsite if desired.')
    
    parser.add_argument('-f',  dest='format', 
                        default='fasta',
                        help='Format of file written to output.')    
    
    print sys.argv
    args = parser.parse_args()
    
    # Write records to output
    db = Reads_db(args.input, recbyname=True)
    
    fastafile_handle = db.write_reads(args.output, output_format=args.format,
                                      filter_expression=args.filter_expression,
                                      startidx=args.startidx)