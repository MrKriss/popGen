#!/usr/bin/env python
# encoding: utf-8

# import modules used here -- sys is a very standard one
import glob
import os
import subprocess
import sys, argparse, logging
from Bio import SeqIO
import random

# Gather code in a main() function
def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    path = os.path.split(args.input[0])[0]

    # format cmd string
    # inputs is a glob, filter so each base fileneame occurs only once
    basenames = [(os.path.split(name)[1]).split('.')[0] for name in args.input]

    # get unique baseneamse
    basenames = list(set(basenames))

    sample_files = ''

    for name in basenames:
        sample_files += ' -s {}'.format(os.path.join(path, name))

    # Run cstacks
    cmd = 'cstacks {samplefiles} -b {batch_id} -o {outpath} -n {num_mismatch} -p {num_threads}'.format(
        samplefiles=sample_files, batch_id=args.batch_id, outpath=args.outputpath, num_mismatch=args.num_mismatch,
        num_threads=args.processors)

    if args.genomic:
        cmd += ' -g'
    if args.appendpath:
        cmd += ' --catalog {}'.format(args.appendpath)

    print cmd
    subprocess.check_call(cmd.split())
    logging.info('Catalogue created for {} and written to\n{}\n'.format(str(basenames), args.outputpath))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Takes multiple stacks files and combines them into a single catalogue.")

    parser.add_argument(
        "-i", dest="input",
        required=True,
        nargs='+',
        help="Location of input stacks files")

    parser.add_argument(
        "-b", dest="batch_id",
        required=True,
        help="ID used for MySQL for this batch.")

    parser.add_argument(
        "-s", dest="subpops",
        default=None,
        nargs='+',
        help="List of string patterns denoting files that make up separate Subpopulations. ")

    parser.add_argument(
        "-n", dest="num_mismatch",
        required=True,
        help="Number of mismatches allowed between sample tags when generating the catalog.")

    parser.add_argument(
        "-g", dest="genomic",
        help="base catalog matching on genomic location, not sequence identity.")

    parser.add_argument(
        "-a", dest="appendpath",
        help="Provide the path to an existing catalog. cstacks will add data to this existing catalog.")

    parser.add_argument(
        "-o",
        dest='outputpath',
        default='.',
        help="Path to write output file to.")

    parser.add_argument(
        "-p", dest="processors", default=1,
        help="Number of processors to run cstacks with.")

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