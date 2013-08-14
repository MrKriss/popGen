#! /usr/bin/env python 
'''
Created on 3 Jul 2013

@author: musselle
'''
import sys, os
import _addpaths
import argparse
import time

from database.reads_db import Reads_db 

# Load in data to SQLite database 
parser = argparse.ArgumentParser(description='Filter and clean up FastQ files.')
parser.add_argument('-i',  dest='input', default = '-', nargs='+',
                    help='Input file(s) to process (/path/filename). Will accept a glob')
parser.add_argument('-b',  dest='barcodes', required=True, nargs='+',
                    help='Barcodes accociated with input file(s) (/path/filename). Will accept a glob')
parser.add_argument('--buffer',  dest='buffer_max', type=int, default=100000,
                    help='Max number of reads that are written in batch to database. Default=1000000.')

# Output parameters
parser.add_argument('-d',  dest='database_filepath', required=True,
                    help='Path and filename of database to writ reads to.')

parser.add_argument('-f',  dest='force_overwrite', action='store_true', default=False,
                    help='Overwrite previous tables with the same name.')


args = parser.parse_args()

toc = time.time()

if args.input == '-':
    args.input = sys.stdin

# Connect/make database
db = Reads_db(db_file=args.database_filepath, recbyname=True)

if ('seqs' not in db.tables) or (args.force_overwrite == True):
    db.create_seqs_table(overwrite=args.force_overwrite)
    
if ('samples' not in db.tables) or (args.force_overwrite == True):
    db.create_samples_table(overwrite=args.force_overwrite)
    
db.load_seqs(data_files=args.input, barcode_files=args.barcodes, buffer_max=args.buffer_max)

total_t = time.time() - toc    
print >> sys.stderr, 'Loaded processed reads file in {0}'.format(
              time.strftime('%H:%M:%S', time.gmtime(total_t)))

if __name__ == '__main__':
    pass