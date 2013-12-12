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

import _addpaths

# Gather code in a main() function
def main(args, loglevel):

    # Setup Logging
    #--------------
    # Log file to record everything
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename= os.path.join(args.out_path, 'run_bowtie-logfile.log'),
                    filemode='a')
    # define a Handler which writes messages to the sys.stderr based on verbosity level
    console = logging.StreamHandler()
    console.setLevel(loglevel)
    # set a format which is simpler for console use
    formatter = logging.Formatter('\n%(asctime)s %(name)s: %(levelname)-8s \n%(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    # Using root now goes to both, using console just goes to console

    # Log parameters passed
    args_str = str([x for x in dir(args) if not x.startswith('_')])
    logging.debug('Arguments passed:\n{}'.format(args_str))


    # MAIN BODY
    #-----------
    # Get generators for filepaths
    logging.debug("Subpops passed {}".format(str(args.in_files)))

    num_infiles = len(args.in_files)

    sqlindex = args.sqlindex_start

    for i, filepath in enumerate(args.in_files):

        # n = max number of mismatchs between
        # l = number of bases over which to test for mismatches
        # --best = report in best to worst order
        logging.info("Processing file {} of {}".format(i, num_infiles))
        print filepath

        sg = SeqIO.parse(filepath, 'fastq')
        len_seq = len(sg.next().seq)

        bowtie_ouput_filename = os.path.splitext(os.path.basename(filepath))[0] + '.bowtie'
        bowtie_ouput_filepath = os.path.join(args.alignment_path, bowtie_ouput_filename)

        bowtie_cmd = 'bowtie -n 3 -l {seqlen} -k {numhits} --best -q -p {threads} {moreargs} {index} {input} {output}'.format(
            seqlen=len_seq, numhits=args.numhits, threads=args.processors, moreargs=args.bowtie_args,
            index=args.index, input=filepath, output=bowtie_ouput_filepath)

        logging.debug("About to run Bowtie with following comandline arguments:\n{}\n".format(str(bowtie_cmd.split())))
        subprocess.check_call(bowtie_cmd.split())
        logging.info("Finished Aligning {} with Bowtie".format(filepath))

        pstacks_pgm = os.path.expanduser(args.stackspath)
        pstacks_pgm = os.path.join(pstacks_pgm ,'pstacks')

        pstacks_cmd = '{pstacks} -p {threads} -t bowtie -f {input} -o {output} -i {sqlindex} -m {min_depth}'.format(
            pstacks=pstacks_pgm, threads=args.processors, input=bowtie_ouput_filepath, output=args.out_path, sqlindex=sqlindex,
            min_depth=args.min_depth)

        if args.model_type:
            pstacks_cmd += ' --model_type {}'.format(args.model_type)
        if args.bound_low:
            pstacks_cmd += ' --bound_low {}'.format(args.bound_low)
        if args.bound_high:
            pstacks_cmd += ' --bound_high {}'.format(args.bound_high)


        logging.debug("About to run pstacks with following comandline arguments:\n{}\n".format(str(pstacks_cmd.split())))
        subprocess.check_call(pstacks_cmd.split())
        sqlindex += 1
        logging.info("Finished Running pstacks for {}".format(filepath))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given a set of files and a bowtie index for a unitag ref, this script runs bowtie to align each '
                    'sample_file to that reference sequence.')

    parser.add_argument(
        "-i", dest="in_files", nargs='+',
        help="Input file set (via a glob) for all processed files to be aligned.")

    parser.add_argument(
        "-a", dest="alignment_path",
        required=True,
        help="Location to write bowtie alignments output")

    parser.add_argument(
        "-o", dest="out_path",
        required=True,
        help="Location to write pstacks output")

    parser.add_argument(
        "-x", dest="index",
        required=True,
        help="File path to bowtie index file")

    parser.add_argument(
        "-q", dest="sqlindex_start",
        default=1,
        type=int,
        help="Starting index for sqlindex.")

    parser.add_argument(
        "-m", dest="min_depth",
        default=1,
        type=int,
        help="Minimum depth of coverage to be considered a stack in pstacks.")

    parser.add_argument(
        "-p", dest="processors", default=1,
        help="Number of processors to run bowtie with.")

    parser.add_argument(
        "-k", dest="numhits", default=1,
        help="Number of alignments to report for each read aligned with bowtie.")

    parser.add_argument(
        "-A", dest="bowtie_args", default='',
        help="Further Arguments to pass to bowtie.")

    parser.add_argument(
        "--stackspath", dest="stackspath", default='~/bin/',
        help="Location of where stacks binaries are installed.")

    parser.add_argument(
        "--model_type", dest="model_type", default='',
        help="Model type to use for snp calls. 'snp', 'bounded' or 'fixed'. ")

    parser.add_argument(
        "--bound_high", dest="bound_high", default='',
        help="upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).")

    parser.add_argument(
        "--bound_low", dest="bound_low", default='',
        help="lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).")



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