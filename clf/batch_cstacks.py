#!/usr/bin/env python
# encoding: utf-8
"""

SCript to run batch processes on


"""

# import modules used here -- sys is a very standard one
import os, sys, argparse, logging
import subprocess
import glob

from collections import defaultdict

# Gather code in a main() function
def main(args, loglevel):
    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    # Load in barcode dictionary
    f = open(args.barcodes, 'rb')
    mid2file = {}
    file2mid = {}
    for line in f:
        line = line.strip().split('\t')
        mid2file[line[0]] = line[1] # Midtag \t filename pairs per line. Need to map filename 2 MID
        file2mid[line[1]] = line[0] # Midtag \t filename pairs per line. Need to map filename 2 MID
    f.close()

    filenames = file2mid.keys()

    stack_files_dict = defaultdict(list)

    for subpop in args.subpops:

        # Filenames and barcodes corresponding to subpopulation
        files = [fname for fname in filenames if subpop in fname]
        barcodes = [file2mid[fname] for fname in filenames]

        # Find all stacks correesponding to the subpopulation
        for b in barcodes:
            stack_files = glob.glob(os.path.join(args.stackspath, 'sample_' + b + '*'))
            stack_files_dict[subpop] = [os.path.split(fp)[1] for fp in stack_files]

        # Run cstacks to form a 'superparant' catalogue
        # cstacks -b batch_id -s sample_file [-s sample_file_2 ...] [-o path] [-n num] [-g] [-p num_threads] [--catalog path] [-h]




    for i in range(100):


    for


    cmd = 'cstacks -b {} '





    print "Hello there."
    logging.info("You passed an argument.")
    logging.debug("Your Argument: %s" % args.argument)


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Wrapper for executing the ustacks program of the STACKS pipeline over multiple files. '
                    'Parameters are the same')

    parser.add_argument(
        "-i",
        help="Input file glob to parse to the ",
        metavar="ARG",
        nargs='+')

    parser.add_argument(
        "-p",
        help="Path to directory containing all stacks.",
        metavar="stackspath",
        nargs='+')


    parser.add_argument(
        "-subpops",
        default=None,
        nargs='+',
        help="List of all barcodes for all file names.")

    parser.add_argument(
        "barcodes",
        help="List of all barcodes for all file names.")

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



