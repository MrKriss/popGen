#!/usr/bin/env python
# encoding: utf-8

# import modules used here -- sys is a very standard one
import os
import sys, argparse, logging
from Bio import SeqIO, Seq

# Gather code in a main() function
def main(args, loglevel):
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    # Append to all sequences in the description
    seqgen = SeqIO.parse(open(args.input), 'fastq')

    outfile = open(os.path.join(args.outputpath, os.path.split(args.input)[1] + '.bar_add'), 'wb')

    writebuffer = []
    buffer_count = 0
    buffer_limit = 10000

    file_basename = os.path.splitext(os.path.split(args.input)[1])[0]
    bar = file_basename.split('_')[1]

    for rec in seqgen:

        # Letter annotations must be removed before editing rec.seq
        phred_scores = rec.letter_annotations['phred_quality']
        rec.letter_annotations = {}

        seq_str = rec.seq.tostring()
        rec.seq = Seq.Seq(bar + seq_str)

        temp_var =  [40]*len(bar) + phred_scores
        dict.__setitem__(rec._per_letter_annotations, "phred_quality", temp_var)

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

    logging.info('Finished Adding barcodes to {}. Written to {}'.format(args.input, outfile.name))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Takes a single preprocessed file from stacks process_RADtags and inserts a dummy barcode on the front so they"
                    "can be inserted into the database. Barcodes are take from the filenames",)

    parser.add_argument(
        "input",
        help="Single input fastq file to process")

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