#!/usr/bin/env python
# encoding: utf-8
"""
clf.make_suparparent.py

@author:     cmusselle

@license:    license

@contact:    user_email
"""

import os, argparse, logging
import subprocess
import glob
from Bio import SeqIO

from collections import Counter
import pandas as pd
import cPickle as pkl


# PROCEDURE #
#-----------#
# Find top x individuals with read coverage.
# For each:
    # Construct a counter for read depth.
    # Store read depth frequencies in descending order
# Plot all frequencies

# Gather code in a main() function
def main(args, loglevel):

    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)
    args_str = str([x for x in dir(args) if not x.startswith('_')])
    logging.debug('Argumnets passed:\n{}'.format(args_str))

    # Load in barcode dictionary
    f = open(args.barcodes, 'rb')
    mid2file = {}
    file2mid = {}
    for line in f:
        line = line.strip().split('\t')
        mid2file[line[0]] = line[1] # Midtag \t filename pairs per line. Need to map filename 2 MID
        file2mid[line[1]] = line[0] # Midtag \t filename pairs per line. Need to map filename 2 MID
    f.close()

    # Find top samples for retained read count
    df = pd.read_csv('processed_radtags.tsv', sep='\t')
    df = df.sort(column='Retained Reads', ascending=0)
    top_samples = df.ix[0:5, ['File', 'Retained Reads']]

    # Populate counters for retained reads
    counter_list = []
    for row in top_samples.iterrows():

        filename = row[1]['File']
        infilepath = os.path.join(os.path.abspath(args.inputpath), 'sample_' + file2mid[filename] + '.fq')

        # Populate counter
        read_counter = Counter()
        seqGen = SeqIO.parse(infilepath, 'fastq')
        for seqRec in seqGen:
            s = seqRec.seq.tostring()
            read_counter[s] += 1

        counter_list.append(read_counter)

    # pickle data to disk
    logging.info('Begining to write to pkl file')
    pkl.dump(counter_list, open(args.outfile_path + '.pkl', 'w'))
    logging.info('Finished writing to pkl file')

    #
    ## Populate counter
    #read_counter = Counter()
    #seqGen = SeqIO.parse(args.inputfile_path, 'fastq')
    #for seqRec in seqGen:
    #    s = seqRec.seq.tostring()
    #    read_counter[s] += 1
    #
    ## Trim counter
    #to_delete = []
    #for k, v in read_counter.iteritems():
    #    if v < args.min or v > args.max:
    #        to_delete.append(k)
    #
    #logging.info('Found {} unique reads.'.format(len(read_counter)))
    #
    #for k in to_delete:
    #    del read_counter[k]
    #
    #logging.info('{} remain after filtering.'.format(len(read_counter)))
    #
    ## Write to file, using a buffer for speed
    #if os.path.exists(args.outputfile_path):
    #    os.remove(args.outputfile_path)
    #outfile = open(args.outputfile_path, 'a')
    #
    #seqGen = SeqIO.parse(args.inputfile_path, 'fastq')
    #seqRec_buffer = []
    #buf_count = 0
    #write_count = 0
    #read_count = 0
    #for seqRec in seqGen:
    #    read_count += 1
    #    s = seqRec.seq.tostring()
    #    if s in read_counter:
    #        seqRec_buffer.append(seqRec)
    #        buf_count += 1
    #        if buf_count >= 1000:
    #            # write batch to file
    #            c = SeqIO.write(seqRec_buffer, outfile, 'fasta')
    #            write_count += c
    #            # Reset buffer
    #            seqRec_buffer = []
    #            buf_count = 0
    #    else:
    #        continue
    #
    ## Flush remainder of buffer
    #if seqRec_buffer:
    #    # write batch to file
    #    c = SeqIO.write(seqRec_buffer, outfile, 'fasta')
    #    write_count += c
    #
    #logging.info('Wrote {} reads out of {} to unitag reference.\n{} skipped due to thresholds.'.format(
    #                                write_count, read_count, read_count-write_count))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Calculate distributions of read depth for top x samples with highest total read count.')

    parser.add_argument(
        "-p", desp="inputpath",
        required=True,
        help="Input path for files preprocessed with stacks.")

    #parser.add_argument(
    #    "-i", desp="inputfiles",
    #    required=True,
    #    nargs='+',
    #    help="Filenames of fastq file to asses.")

    parser.add_argument(
        "-o", dest="outfile_path",
        default='read_depth_data',
        help="Location and file name of output file to write read depth data to.")

    parser.add_argument(
        "-x", dest="topx",
        default=5, type=int,
        help="Number of top samples to plot.")

    parser.add_argument(
        "-b", dest="barcodes",
        required=True,
        help="File storing filename barcode pairs.")

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