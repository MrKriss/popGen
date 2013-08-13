#! /usr/bin/env python 
'''
Created on 3 Jul 2013

@author: musselle
'''
import sys, os
import _addpaths
import argparse

from database.reads_db import Reads_db 

# Load in data to SQLite database 
parser = argparse.ArgumentParser(description='Filter and clean up FastQ files.')
parser.add_argument('-i',  dest='input', default = '-', nargs='+',
                    help='Input file(s) to process. (/path/filename) Will accept a glob')
parser.add_argument('-b',  dest='barcodes', required=True, nargs='+',
                    help='Barcodes accociated with input file(s). Will accept a glob')
# Output parameters
parser.add_argument('-d',  dest='database_filepath', required=True,
                    help='Path and filename of database to writ reads to.')

parser.add_argument('-f',  dest='force_overwrite', action='store_true', default=False,
                    help='Overwrite previous tables with the same name.')

args = parser.parse_args()

if args.input == '-':
    args.input = sys.stdin

# Connect/make database
db = Reads_db(db_file=args.database_filepath, recbyname=True)

if ('seqs' not in db.tables) or (args.force_overwrite == True):
    db.create_seqs_table(overwrite=args.force_overwrite)
    
if ('samples' not in db.tables) or (args.force_overwrite == True):
    db.create_samples_table(overwrite=args.force_overwrite)
    
db.load_seqs(data_files=args.input, barcode_files=args.barcodes)

if __name__ == '__main__':
    pass