#! /usr/bin/env python 
'''
Created on 3 Jul 2013

@author: musselle
'''
import os 
import sys 
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

from database.reads_db import Reads_db 

import argparse

# Load in data to SQLite database 
parser = argparse.ArgumentParser(description='Filter and clean up FastQ files.')
parser.add_argument('-i',  dest='input', default = '-', nargs='+',
                    help='Input file(s) to process. (/path/filename) Will accept a glob')
parser.add_argument('-b',  dest='barcodes', required=True, nargs='+',
                    help='Barcodes accociated with input file(s). Will accept a glob')
# Output parameters
parser.add_argument('-d',  dest='database_filepath', required=True,
                    help='Path and filename of database to writ reads to.')


parser.add_argument('-t',  dest='table_name', default='seqs',
                    help='Name of the table to add the sequences too. Default = "seqs"')
parser.add_argument('-f',  dest='force_overwrite', action='store_true', default=False,
                    help='Overwrite previous tables with the same name.')

args = parser.parse_args()

if args.input == '-':
    args.input = sys.stdin

# Connect/make database
db = Reads_db(db_file=args.database_filepath, recbyname=True)

if (args.table_name not in db.tables) or (args.force_overwrite == True):
    db.create_seqs_table(table_name=args.table_name, overwrite=args.force_overwrite)

db.load_seqs(data_files=args.input, barcode_files=args.barcodes, table_name=args.table_name)


if __name__ == '__main__':
    pass