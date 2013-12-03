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

    # Goal #
    # Find unique sequences in the fastq file that have between x and y coverage

    # Method:
        # Two pass: first to count occurances of reads, second to write to file and report statistics

    # INPUT: processed sample_AAATG.fq file

    # OUTPUT: Unitag sample.fq file + logfile.log

# Gather code in a main() function
def main(args, loglevel):

    # Defined Path Vars
    out_path = os.path.split(args.outfile_path)[0]
    in_path = os.path.split(args.infile_path)[0]

    # Setup Logging
    #--------------

    # Log file to record everything
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename= os.path.join(out_path, 'logfile.log'),
                    filemode='a')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    # Using root now goes to both, using console just goes to console

    # Log parameters passed
    args_str = str([x for x in dir(args) if not x.startswith('_')])
    logging.debug('Argumnets passed:\n{}'.format(args_str))

    # Populate counter
    read_counter = Counter()
    seqGen = SeqIO.parse(args.infile_path, 'fastq')
    for seqRec in seqGen:
        s = seqRec.seq.tostring()
        read_counter[s] += 1

    # Trim counter
    too_many = []
    too_few = []
    for k, v in read_counter.iteritems():
        if v < args.min:
            too_few.append(k)
        elif v > args.max:
            too_many.append(k)

    total_reads = sum(read_counter.values())
    retained_reads = total_reads - len(too_few) - len(too_many)

    # Log Stats
    logging.info('\n-----------------------------------------'
                 '\nUnitag Constructed with:'
                 '\nMIN = {min}'
                 '\nMAX = {max}'
                 '\n-----------------------------------------'
                 '\nTotal Reads:\t\t{0}'
                 '\nTotal Unique Reads:\t{1}'
                 '\nReads Fewer Than MIN:\t{2}\t({3:.2%})'
                 '\nReads Greater than MAX\t{4}\t({5:.2%})'
                 '\nRetained Reads:\t\t{6}\t{7:.2%}'.format(
                        sum(read_counter.values()),
                        len(read_counter),
                        len(too_few), float(len(too_few))/ total_reads,
                        len(too_many), float(len(too_many))/ total_reads,
                        retained_reads, float(retained_reads)/ total_reads),
                        min=args.min, max=args.max)

    for k in too_few:
        del read_counter[k]
    for k in too_many:
        del read_counter[k]

    # Write to file, using a buffer for speedup
    #---------------------------------------
    # Each unique sequence is only written once.

    # Check if it already exists
    if os.path.exists(args.outfile_path):
        os.remove(args.outfile_path)
    outfile = open(args.outfile_path, 'a')

    seqGen = SeqIO.parse(args.infile_path, 'fastq')
    seqRec_buffer = []
    buf_count = 0
    write_count = 0
    read_count = 0

    # Full list of seqs to write, once each.
    list_of_seqs = read_counter.keys()

    for seqRec in seqGen:
        read_count += 1
        s = seqRec.seq.tostring()
        if s in list_of_seqs:

            # Remove seq from list
            del list_of_seqs[list_of_seqs.index(s)]

            # Add to buffer
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

    logging.info('\nWrote {} reads out of {} to unitag reference.\n{} skipped due to thresholds.'.format(
                                    write_count, read_count, read_count-write_count))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Constructs a unitag reference from all unique sequences in a specified input file.')

    parser.add_argument(
        "-i", dest="infile_path",
        required=True,
        help="Location and file name of fastq file to use to construct unitag reference.")

    parser.add_argument(
        "-o", dest="outfile_path",
        required=True,
        help="Location and file name of output file to write unitag reference to.")

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
        help="Increase output verbosity",
        action="store_true")

    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)