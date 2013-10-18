#!/usr/bin/env python
# encoding: utf-8
"""
clf.make_suparparent.py

@author:     cmusselle

@license:    license

@contact:    user_email

"""

# import modules used here -- sys is a very standard one
import os, argparse, logging
import subprocess
import glob
from Bio import SeqIO

from collections import Counter

# Gather code in a main() function
def main(args, loglevel):
    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    args_str = str([x for x in dir(args) if not x.startswith('_')])
    logging.debug('Argumnets passed:\n{}'.format(args_str))

    # Goal #
    # Find unique sequences in the fastq file that have between x and y coverage

    # Method:
        # Two pass: first to count occurances of reads, second to write to file

    # INPUT: processed sample_AAATG.fq file

    # Populate counter
    read_counter = Counter()
    seqGen = SeqIO.parse(args.inputfile_path, 'fastq')
    for seqRec in seqGen:
        s = seqRec.seq.tostring()
        read_counter[s] += 1

    # Trim counter
    to_delete = []
    for k, v in read_counter.iteritems():
        if v < args.min or v > args.max:
            to_delete.append(k)

    logging.info('Found {} unique reads.'.format(len(read_counter)))

    for k in to_delete:
        del read_counter[k]

    logging.info('{} remain after filtering.'.format(len(read_counter)))

    # Write to file, using a buffer for speed
    if os.path.exists(args.outputfile_path):
        os.remove(args.outputfile_path)
    outfile = open(args.outputfile_path, 'a')

    seqGen = SeqIO.parse(args.inputfile_path, 'fastq')
    seqRec_buffer = []
    buf_count = 0
    write_count = 0
    read_count = 0
    for seqRec in seqGen:
        read_count += 1
        s = seqRec.seq.tostring()
        if s in read_counter:
            seqRec_buffer.append(seqRec)
            buf_count += 1
            if buf_count >= 1000:
                # write batch to file
                c = SeqIO.write(seqRec_buffer, outfile, 'fasta')
                write_count += c
                # Reset buffer
                seqRec_buffer = []
                buf_count = 0
        else:
            continue

    # Flush remainder of buffer
    if seqRec_buffer:
        # write batch to file
        c = SeqIO.write(seqRec_buffer, outfile, 'fasta')
        write_count += c

    logging.info('Wrote {} reads out of {} to unitag reference.\n{} skipped due to thresholds.'.format(
                                    write_count, read_count, read_count-write_count))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Constructs a unitag reference from all unique sequences in a specified input file.')

    parser.add_argument(
        "-i", "--inputfile_path",
        required=True,
        help="Location and file name of fastq file to use to construct unitag reference.")

    parser.add_argument(
        "-o", "--outputfile_path",
        required=True,
        help="Location and file name of output file to write unitag reference to.")

    #parser.add_argument(
    #    "-b", "--barcodes",
    #    required=True,
    #    help="Barcode file to use for mapping mid to filenames.")

    #parser.add_argument(
    #    "-s", "--subpops",
    #    default=None,
    #    nargs='+',
    #    help="List of string patterns denoting files that make up separate Subpopulations. ")
    #
    #parser.add_argument(
    #    "-u", "--sup_parent_path",
    #    required=True,
    #    help="Location to write superparent output")
    #
    #parser.add_argument(
    #    "-t", "--stack_path",
    #    required=True,
    #    help="Location to write ustacks output")
    #
    #parser.add_argument(
    #    "-p", "--processors",
    #    help="Number of processors to run ustacks with.",
    #    default=1)

    parser.add_argument(
        "-m", "--min",
        type=int,
        help="Minimum depth of coverage allowed for a unique sequence to be used in unitag reference.",
        default=None)

    parser.add_argument(
        "-M", "--max",
        type=int,
        help="Maximum depth of coverage allowed for a unique sequence to be used in unitag reference",
        default=None)

    parser.add_argument(
        "-v", "--verbose",
        help="increase output verbosity",
        action="store_true")

    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)