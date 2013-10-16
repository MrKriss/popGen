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

# Gather code in a main() function
def main(args, loglevel):
    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    logging.debug('Argumnets passed:\n{}'.format(str(dir(args))))

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

    for i, subpop in enumerate(args.subpops):

        # Filenames and barcodes corresponding to subpopulation
        files = [fname for fname in filenames if subpop in fname]
        barcodes = [file2mid[fname] for fname in files]

        assert len(barcodes) == len(files)

        # Find all processed files correesponding to the subpopulation
        processed_files = []
        for b in barcodes:
            processed_files.extend(glob.glob(os.path.join(args.processed_files_path, 'sample_' + b + '*')))

        # Concatonate files
        merged_filepath = os.path.join(args.sup_parent_path, 'superparent_' + subpop)
        cat_cmd = 'cat {}'.format(' '.join(processed_files))
        print cat_cmd
        with open(merged_filepath, 'wb') as merged_file:
            subprocess.check_call(cat_cmd.split(), stdout=merged_file)
        logging.info('Superparent created for {} and written to\n{}\n'.format(subpop, args.sup_parent_path))

        # Run ustacks on the Superparent
        # ustacks -t file_type -f file_path [-d] [-r] [-o path] [-i id] [-m min_cov] [-M max_dist] [-p num_threads]
        # [-R] [-H] [-h]
        ustacks_cmd = 'ustacks -t fastq -f {} -o {} -i {} -p {} -r -d'.format(merged_filepath,
                                                                              args.stack_path,
                                                                              i,
                                                                              args.processors)
        # Add extra options if present
        if args.mindepth:
            ustacks_cmd += '-m {}'.format(args.mindepth)
        if args.maxdist_between_stacks:
            ustacks_cmd += '-M {}'.format(args.maxdist_between_stacks)
        if args.maxdist_second2primary:
            ustacks_cmd += '-N {}'.format(args.maxdist_second2primary)

        logging.debug("About to run ustacks with following comandline arguments:\n{}\n".format(str(ustacks_cmd.split())))
        subprocess.check_call(ustacks_cmd.split())
        logging.info("Finished processing {} of {} subpopulations.".format(i, len(args.subpops)))

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Constructs a "superparent" from all sequences in a specified subpopulation by concatonating them '
                    'and running on ustacks. This "superparent" can then be used to create a single catalogue using '
                    'cstacks.')

    parser.add_argument(
        "-i", "--processed_files_path",
        required=True,
        help="Location of preprocessed input files")

    parser.add_argument(
        "-b", "--barcodes",
        required=True,
        help="Barcode file to use for mapping mid to filenames.")

    parser.add_argument(
        "-s", "--subpops",
        default=None,
        nargs='+',
        help="List of string patterns denoting files that make up separate Subpopulations. ")

    parser.add_argument(
        "-u", "--sup_parent_path",
        required=True,
        help="Location to write superparent output")

    parser.add_argument(
        "-t", "--stack_path",
        required=True,
        help="Location to write ustacks output")

    parser.add_argument(
        "-p", "--processors",
        help="Number of processors to run ustacks with.",
        default=1)

    parser.add_argument(
        "-m", "--mindepth",
        help="Minimum depth of coverage required to create a stack in ustacks (default 2).",
        default=None)

    parser.add_argument(
        "-M", "--maxdist_between_stacks",
        help="Maximum distance (in nucleotides) allowed between stacks in ustacks (default 2).",
        default=None)

    parser.add_argument(
        "-N", "--maxdist_second2primary",
        help="Maximum distance allowed to align secondary reads to primary stacks in ustacks (default: M + 2).",
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