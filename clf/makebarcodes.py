#!/usr/bin/env python
# encoding: utf-8

# import modules used here -- sys is a very standard one
import os
import sys, argparse, logging
from Bio import SeqIO
import random

# Gather code in a main() function
def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    num_files = len(args.inputs)
    bases = 'ATGC'
    barcode_dict = {}

    # For each file, generate a barcode to asign to it.
    for i in range(num_files):

        filename = args.inputs[i]

        while True:
            # Generate barcode
            tag = ''.join(random.choice(bases) for x in range(6))
            # check its different from all others
            if not tag in barcode_dict:
                barcode_dict[tag] = filename
                break

    # write barcodes only (for input to stacks)
    barcodes_only = open(os.path.join(args.outputpath, 'barcodes_only.txt'), 'wb')
    barcodes_only.writelines("\n".join(barcode_dict.keys()))
    barcodes_only.flush()
    barcodes_only.close()

    logging.info('Written barcodes to processing {}'.format(barcodes_only.name))

    # Write barcode filename pairs for reference later
    barcodes_filenams = open(os.path.join(args.outputpath, 'barcodes_filenames.txt'), 'wb')
    for k, v in barcode_dict.iteritems():
        barcodes_filenams.write(k + '\t' + str(v) + '\n')
    barcodes_filenams.flush()
    barcodes_filenams.close()

    logging.info('Written barcode-filename pairs to processing {}'.format(barcodes_filenams.name))

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Randomly generates a sequence for MID barcodes for each file parsed, and writes the pairs to a "
                    "tab delimited text file, 1 pair per line.")

    parser.add_argument(
        "inputs",
        help="List of inputs files to generate barcodes for.",
        nargs='+')

    parser.add_argument(
        "-o",
        dest='outputpath',
        default='',
        help="Path to write output files to.")

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