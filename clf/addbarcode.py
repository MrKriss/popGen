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
    barcode_dict ={}

    # For each file, generate a barcode, append it to all sequence descriptions in the file
    for i in range(num_files):

        filename = args.inputs[i]


        while True:
            # Generate barcode
            tag = ''.join(random.choice(bases) for x in range(6))
            # check its different from all others
            if not tag in barcode_dict:
                barcode_dict[tag] = filename
                break

                # Append to all sequences in the description
        seqgen = SeqIO.parse(open(filename), 'fastq')

        outfile = open(os.path.join(args.outputpath, filename + '.temp'), 'wb')

        buffer = []
        buffer_count = 0
        buffer_limit = 10000

        for seq in seqgen:

            seq.description += tag
            buffer.append(seq)
            buffer_count += 1
            if buffer_count > buffer_limit:
                SeqIO.write(buffer, outfile, 'fastq')
                buffer = []
                buffer_count = 0

        if buffer:
            SeqIO.write(buffer, outfile, 'fastq')

        outfile.flush()
        outfile.close()

        logging.info('Finished processing {}. Written to {}'.format(filename, outfile.name))

    barcodes = open(os.path.join(args.outputpath, 'barcodes.txt'), 'wb')
    barcodes.writelines(barcode_dict.keys())
    barcodes.flush()
    barcodes.close()



# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Takes a list of preprocessed files from Floragenics and inserts a dummy barcode on the front so they can be "
                    "rerun through a preprocessing step using STACKS.",)

    parser.add_argument(
        "inputs",
        help="List of inputs files to process",
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