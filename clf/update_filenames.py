#!/usr/bin/env python
# encoding: utf-8

import glob
import os
import subprocess
import sys, argparse, logging
from Bio import SeqIO
import random

import _addpaths

# Gather code in a main() function
def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: \n%(message)s", level=loglevel)

    # Load in barcode dictionary
    f = open(args.barcodes, 'rb')
    mid2file = {}
    file2mid = {}
    for line in f:
        line = line.strip().split('\t')
        mid2file[line[0]] = line[1] # Midtag \t filename pairs per line. Need to map filename 2 MID
        file2mid[line[1]] = line[0] # Midtag \t filename pairs per line. Need to map filename 2 MID
    f.close()

    for i, name in enumerate(args.input):

        # extract mid tag
        tag, ext = os.path.splitext(os.path.basename(name)).split('_')

        short_filename = mid2file[tag].split('_')[2]

        # Look up what original filename is
        new_name = short_filename + '_' + tag + '.' + ext

        new_filepath = os.path.join(os.path.dirname(name), new_name)

        os.rename(args.input[i], new_filepath)
        logging.debug('Renames {} to\n{}'.format(args.input[i], new_filepath)

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Takes multiple stacks files and renames them accoring to their original files names defined in "
                    "the barcodes_filenames.txt file")

    parser.add_argument(
        "-i", dest="input",
        required=True,
        nargs='+',
        help=" Glob of multiple input stacks files inthe format of 'sample_XXXXXX.fq'")

    parser.add_argument(
        "-b", dest="barcodes",
        required=True,
        help="File mapping barcodes to filenames.")

    parser.add_argument(
        "-o",
        dest='outputpath',
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