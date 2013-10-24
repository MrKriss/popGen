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


# Gather code in a main() function
def main(args, loglevel):
    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    logging.debug('Argumnets passed:\n{}'.format(str(dir(args))))


    # Get generators for filepaths
    if 'all' in args.subpops:

        # Generator to get all processed files.
        file_gen = (x for x in glob.glob(os.path.join(args.infilepath, "sample_*")))
        file_gen_list = [file_gen]
    else:
        # Load in barcode dictionary
        f = open(args.barcodes, 'rb')
        mid2file = {}
        file2mid = {}
        for line in f:
            line = line.strip().split('\t')
            mid2file[line[0]] = line[1] # Midtag \t filename pairs per line. Need to map filename 2 MID
            file2mid[line[1]] = line[0] # Midtag \t filename pairs per line. Need to map filename 2 MID
        f.close()

        # List of raw filenames
        filenames = file2mid.keys()

        file_gen_list = []
        for i, subpop in enumerate(args.subpops):

            # Filenames and barcodes corresponding to subpopulation
            files = [fname for fname in filenames if subpop in fname]
            barcodes = [file2mid[fname] for fname in files]

            file_gen = (glob.glob(os.path.join(args.infilepath, 'sample_' + b + '*')) for b in barcodes)
            file_gen_list.append(file_gen)

    # Align file with bowtie to index, then run on pstacks
    sqlindex = 1
    for gen in file_gen_list:
        for filepath in gen:

            # n = max number of mismatchs between
            # l = number of bases over which to test for mismatches
            # --best = report in best to worst order

            sg = SeqIO.parse(filepath, 'fastq')
            len_seq = len(sg.next().seq)

            bowtie_ouput_filename = os.path.splitext(os.path.basename(filepath))[0] + '.bowtie'
            bowtie_cmd = 'bowtie -n 3 -l {seqlen} -k {numhits} --best -q -p {threads} {moreargs} {index} {input} {output}'.format(
                seqlen=len_seq, numhits=args.numhits, threads=args.processors, moreargs=args.bowtie_args,
                index=args.index, input=filepath, output=bowtie_ouput_filename)

            logging.debug("About to run Bowtie with following comandline arguments:\n{}\n".format(str(bowtie_cmd.split())))
            subprocess.check_call(bowtie_cmd.split())
            logging.info("Finished Aligning {} with Bowtie".format(filepath))


            pstacks_cmd = 'pstacks -p {threads} -t bowtie -f {input} -o {output} -i {sqlindex}'.format(
                threads=args.processors, input=bowtie_ouput_filename, output=args.out_path, sqlindex=sqlindex)
            logging.debug("About to run pstacks with following comandline arguments:\n{}\n".format(str(pstacks_cmd.split())))
            subprocess.check_call(pstacks_cmd.split())
            sqlindex += 1
            logging.info("Finished Running pstacks for {}".format(filepath))


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given a subpopulation and a bowtie index for a unitag ref, this script runs bowtie to align each '
                    'sample_file for the subpopulation to that reference sequence.')

    parser.add_argument(
        "-s", dest="subpops",
        required=True,
        nargs='+',
        help="List of string patterns denoting files that make up separate Subpopulations. 'all' will align all samples"
             "from all subpopulations. ")

    parser.add_argument(
        "-i", dest="infilepath",
        help="Input file path for all processed files if not specifying a subpop.")

    parser.add_argument(
        "-x", dest="index",
        required=True,
        help="Location of superparent input files")

    parser.add_argument(
        "-b", dest="barcodes",
        help="Barcode file to use for mapping mid to filenames.")

    parser.add_argument(
        "-o", dest="out_path",
        required=True,
        help="Location to write ustacks output")

    parser.add_argument(
        "-p", dest="processors", default=1,
        help="Number of processors to run bowtie with.")

    parser.add_argument(
        "-k", dest="numhits", default=1,
        help="Number of alignments to report for each read aligned with bowtie.")

    parser.add_argument(
        "-a", dest="bowtie_args", default='',
        help="Further Arguments to pass to bowtie.",
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