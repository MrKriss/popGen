#!/usr/bin/env python
# encoding: utf-8

# import modules used here -- sys is a very standard one
import os
import sys, argparse, logging
from Bio import SeqIO, Seq
import random

# Gather code in a main() function
def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    # Read in barcodes
    f = open(args.barcodes, 'rb')
    barcode_dict = {}
    for line in f:
        line = line.strip().split('\t')
        barcode_dict[line[1]] = line[0] # Midtag \t filename pairs per line. Need tomap filename 2 MID
    f.close()

    assert args.input in barcode_dict, "{} not found in barcode_dict".format(args.input)

    # Append to all sequences in the description
    seqgen = SeqIO.parse(open(args.input), 'fastq')

    outfile = open(os.path.join(args.outputpath, args.input + '.temp'), 'wb')

    writebuffer = []
    buffer_count = 0
    buffer_limit = 10000

    for rec in seqgen:

        seq_str = rec.seq.tostring()

        # Letter annotations must be removed before editing rec.seq
        temp_var = rec.letter_annotations
        rec.letter_annotations['phred_quality'] = [40]*len(barcode_dict[args.input]) \
                                                  + rec.letter_annotations['phred_quality']

        rec.seq = Seq.Seq(barcode_dict[args.input] + seq_str)

        writebuffer.append(rec)
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
        help="List of all barcode-filename pairs for all file names, one per line.")

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