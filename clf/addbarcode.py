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

    # Read in barcodes
    f = open(args.barcodes, 'rb')
    barcode_dict ={}
    for line in f:
        line = line.strip().split('\t')
        barcode_dict[line[1]] = line[0] # Midtag \t filename pairs per line. Need tomap filename 2 MID
    f.close()

    assert args.input in barcode_dict

    # Append to all sequences in the description
    seqgen = SeqIO.parse(open(args.input), 'fastq')

    outfile = open(os.path.join(args.outputpath, args.input + '.temp'), 'wb')

    writebuffer = []
    buffer_count = 0
    buffer_limit = 10000

    for seq in seqgen:
        seq.description += barcode_dict[args.input]
        writebuffer.append(seq)
        buffer_count += 1
        if buffer_count > buffer_limit:
            SeqIO.write(writebuffer, outfile, 'fastq')
            writebuffer = []
            buffer_count = 0

    if writebuffer:
        SeqIO.write(writebuffer, outfile, 'fastq')

    outfile.flush()
    outfile.close()

    logging.info('Finished processing {}. Written to {}'.format(args.input, outfile.name))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Takes a single preprocessed file from Floragenics and inserts a dummy barcode on the front so they"
                    "can be rerun through a preprocessing step using STACKS.",)

    parser.add_argument(
        "input",
        help="Single input fastq file to process")

    parser.add_argument(
        "barcodes",
        help="List of all barcodes for all file names.")

    parser.add_argument(
        "-o",
        dest='outputpath',
        default='',
        help="Path to write output file to.")

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