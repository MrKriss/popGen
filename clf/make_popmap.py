#!/usr/bin/env python
# encoding: utf-8

import glob
import os
import subprocess
import sys, argparse, logging
from Bio import SeqIO
import random

import _addpaths

import pandas as pd

# Gather code in a main() function
def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: \n%(message)s", level=loglevel)

    # Numerical values are assigned according to alphabetical order of sup populations
    df1 = pd.read_csv(args.barcodes, sep='\t', header=None, names=['bar', 'file'])
    df2 = df1.copy()
    df2['file'] = [x[17:19] for x in df1['file']]

    unique_sub_pops = df2['file'].unique()
    unique_sub_pops = sorted(unique_sub_pops)

    # mapping dictionary
    D = dict(zip(unique_sub_pops, range(1, len(unique_sub_pops)+1)))

    print D

    # filter inputs to only unique entries
    basenames = [os.path.basename(name).split('.')[0] for name in args.input]
    # get unique baseneamse
    basenames = list(set(basenames))

    f_out = open(args.out_filepath, 'w')

    for name in basenames:
        subpop = name.split('_')
        subpop = subpop[0].strip('1234567890')
        f_out.write( name + '\t' + D[subpop]  + '\n')
        logging.debug('Wrote {} to file:\n{}'.format(name + '\t' + D[subpop]  + '\n', args.out_filepath))

    f_out.close()


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Takes multiple stacks files and creates a population map accoring to their original files names"
                    " defined in "
                    "the barcodes_filenames.txt file")

    parser.add_argument(
        "-i", dest="input",
        required=True,
        nargs='+',
        help="Glob of multiple input processed stacks files inthe format of 'sample_XXXXXX.fq'")

    parser.add_argument(
        "-b", dest="barcodes",
        required=True,
        help="File mapping barcodes to filenames.")

    parser.add_argument(
        "-o",
        dest='out_filepath',
        default='.',
        help="Optional Path to write output file to.")

    parser.add_argument(
        "-v",
        "--verbose",
        help="increase output verbosity",
        action="store_true")

    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)