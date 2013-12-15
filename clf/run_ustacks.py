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
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)
    path = os.path.split(args.input[0])[0]

    # format cmd string
    # inputs is a glob, filter so each base fileneame occurs only once
    basenames = [os.path.split(name)[1] for name in args.input]

    # get unique baseneamse
    basenames = list(set(basenames))

    basenames = sorted(basenames)

    sqlindex = args.sqlindex_start

    for name in basenames:

        sample_file = os.path.join(path, name)

        pgm_filepath = os.path.expanduser(args.stackspath)
        pgm_filepath = os.path.join(pgm_filepath, 'ustacks')

        # Run ustacks
        cmd = '{pgm_filepath} -t fastq -f {samplefiles} -i {sqlindex} -o {outpath} -m {min_depth} ' \
              '-M {num_mismatch} -p {num_threads} -r -d --max_locus_stacks {max_locus_stacks} ' \
              '--model_type bounded --bound_high 0.1'.format(
            pgm_filepath=pgm_filepath, samplefiles=sample_file, sqlindex=sqlindex, outpath=args.outputpath,
            num_mismatch=args.num_mismatch, min_depth=args.min_depth, num_threads=args.processors,
            max_locus_stacks=args.max_locus_stacks)

        logging.debug("About to run Ustacks with following comandline arguments:\n{}\n".format(str(cmd.split())))
        subprocess.check_call(cmd.split())
        logging.info('Stacks created for {} and written to\n{}\n'.format(str(basenames), args.outputpath))
        sqlindex += 1

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Takes multiple processed fastq files and constructs stacks using ustacks.")

    parser.add_argument(
        "-i", dest="input",
        required=True,
        nargs='+',
        help="Location of input files. Accepts a glob of multiple files.")

    parser.add_argument(
        "-q", dest="sqlindex_start",
        default=1,
        type=int,
        help="Starting index for sqlindex.")

    #parser.add_argument(
    #    "-s", dest="subpops",
    #    default=None,
    #    nargs='+',
    #    help="List of string patterns denoting files that make up separate Subpopulations. ")


    parser.add_argument(
        "-m", dest="min_depth",
        help="Minimum read depth required to form a stack.")

    parser.add_argument(
        "-M", dest="num_mismatch",
        required=True,
        help="Number of mismatches allowed between reads in a stack.")

    parser.add_argument(
        "-L", dest="max_locus_stacks",
        help="Maximum number of stacks that can be combined into a single locus.")

    parser.add_argument(
        "-o",
        dest='outputpath',
        default='.',
        help="Path to write output file to.")

    parser.add_argument(
        "-p", dest="processors", default=1,
        help="Number of processors to run cstacks with.")

    parser.add_argument(
        "--stackspath", dest="stackspath", default='~/bin/',
        help="Location of where stacks binaries are installed.")

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