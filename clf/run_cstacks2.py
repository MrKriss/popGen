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

    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)
    args_str = str([x for x in dir(args) if not x.startswith('_')])
    logging.debug('Argumnets passed:\n{}'.format(args_str))


    path = os.path.split(args.input[0])[0]

    # format cmd string
    # inputs is a glob, filter so each base fileneame occurs only once
    basenames = ['.'.join((os.path.split(name)[1]).split('.')[:-2]) for name in args.input]

    # get unique baseneamse
    basenames = list(set(basenames))


    if not args.catalog_path:
        # Initialise the catalogue with the first sample
        name = basenames.pop()
        sample_file = ' -s {}'.format(os.path.join(path, name))

        pgm_filepath = os.path.expanduser(args.stackspath)
        pgm_filepath = os.path.join(pgm_filepath, 'cstacks')


        # Run cstacks
        cmd = '{pgm_filepath} {samplefile} -b {batch_id} -o {outpath} -n {num_mismatch} -p {num_threads}'.format(
            pgm_filepath=pgm_filepath, samplefiles=sample_file, batch_id=args.batch_id, outpath=args.outputpath,
            num_mismatch=args.num_mismatch, num_threads=args.processors)

        logging.debug("\nAbout to run cstacks to construct initial catalogue with following commandline "
                      "arguments:\n{}\n".format(str(cmd.split())))
        subprocess.check_call(cmd.split())

        args.catalog_path = os.path.join(args.outputpath, 'batch_' + args.batch_id + 'catalog')
        logging.info("\nFinished constructing catalogue and written to {}".format(args.catalog_path))

    for name in basenames:

        sample_file += ' -s {}'.format(os.path.join(path, name))

        pgm_filepath = os.path.expanduser(args.stackspath)
        pgm_filepath = os.path.join(pgm_filepath, 'cstacks')

        # Run cstacks
        cmd = '{pgm_filepath} {samplefiles} -b {batch_id} -n {num_mismatch} -p {num_threads}' \
              '--catalog {catalog_path}'.format(pgm_filepath=pgm_filepath, samplefiles=sample_file,
                    batch_id=args.batch_id, num_mismatch=args.num_mismatch, num_threads=args.processors,
                    catalog_path=args.catalog_path)

        if args.genomic:
            cmd += ' -g'

        logging.debug("\nAbout to run cstacks to construct initial catalogue with following commandline "
                      "arguments:\n{}".format(str(cmd.split())))
        subprocess.check_call(cmd.split())
        logging.info("\nFinished appending {} to catalogue and written to\n{}".format(name, args.catalog_path))

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Takes multiple stacks files and INCREMENTALLY combines them into a single catalogue.")

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
        "-c", dest="catalog_path",
        default=None,
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