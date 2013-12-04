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
                    filename= os.path.join(out_path, 'unitageref-logfile.log'),
                    filemode='a')
    # define a Handler which writes messages to the sys.stderr based on verbosity level
    console = logging.StreamHandler()
    console.setLevel(loglevel)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(asctime)s %(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    # Using root now goes to both, using console just goes to console

    # Log parameters passed
    args_str = str([x for x in dir(args) if not x.startswith('_')])
    logging.debug('Arguments passed:\n{}'.format(args_str))

    # Populate counter
    read_counter = Counter()
    seqGen = SeqIO.parse(args.infile_path, 'fastq')
    for seqRec in seqGen:
        s = seqRec.seq.tostring()
        read_counter[s] += 1

    logging.debug('Finished populating counter, about to trim counter.')

    # Trim counter
    too_many = set()
    too_few = set()
    for k, v in read_counter.iteritems():
        if v < args.min:
            too_few.add(k)
        elif v > args.max:
            too_many.add(k)

    total_reads = sum(read_counter.values())

    for k in too_few:
        del read_counter[k]
    for k in too_many:
        del read_counter[k]

    logging.debug('Finished trimming counter.')

    # count reads and calculate percentages
    all_unique_read_set = set(read_counter.keys())
    total_unique_reads = len(all_unique_read_set)
    retained_unique_reads = total_unique_reads - len(too_few) - len(too_many)
    retained_read_set = all_unique_read_set.difference(too_few, too_many)

    print len(retained_read_set)
    print retained_unique_reads

    assert len(retained_read_set) == retained_unique_reads



    # calculate sum of reads in retained portion
    total_reads_in_retained_unique_reads = 0
    for key, value in read_counter.iteritems():
        if key not in too_few and key not in too_many:
            total_reads_in_retained_unique_reads += value

    # Log Stats
    logging.info("""\n-----------------------------------------
                 Unitag Constructed with:
                 MIN = {min}
                 MAX = {max}
                 Written to {filepath}
                 --------------------------------------------
                 Total Reads:\t\t{total}
                 Total Unique Reads:\t{total_u}
                 Unique Reads Fewer Than MIN:\t{fewer}\t({fewer_p:.2%})
                 Unique Reads Greater than MAX\t{more}\t({more_p:.2%})
                 Unique Reads Retained:\t\t{retained}\t{retained_p:.2%}
                 Total Reads in Retained Unique Reads:\t{retained_t}\t{retained_t_p}""".format(
                        min=args.min, max=args.max,
                        filepath=args.outfile_path,
                        total=total_reads,
                        total_u=len(read_counter),
                        fewer=len(too_few), fewer_p=float(len(too_few))/ total_unique_reads,
                        more=len(too_many), more_p=float(len(too_many))/ total_unique_reads,
                        retained=retained_unique_reads,
                        retained_p=float(retained_unique_reads)/ total_unique_reads,
                        retained_t=total_reads_in_retained_unique_reads,
                        retained_t_p=float(total_reads_in_retained_unique_reads)/ total_reads
                                                                            )
                )

    logging.debug('About to start writing to file...')

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

    logging.info("""\nWrote {} unique reads out of {} total reads to unitag reference.'
                 '{} unique reads skipped due to thresholds.
                 ----------------------------------------------------------------------""".format(
                                    write_count, read_count, total_unique_reads-write_count))


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
        "-m", dest="min",
        type=int,
        help="Minimum depth of coverage allowed for a unique sequence to be used in unitag reference.",
        default=None)

    parser.add_argument(
        "-M", dest="max",
        type=int,
        help="Maximum depth of coverage allowed for a unique sequence to be used in unitag reference",
        default=None)

    parser.add_argument(
        "-v", dest="verbose",
        help="Increase output verbosity",
        action="store_true")

    args = parser.parse_args()

    # Setup logging
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    main(args, loglevel)