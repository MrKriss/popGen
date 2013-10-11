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

# Gather code in a main() function
def main(args, loglevel):
    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)

    for i, subpop in enumerate(args.subpops):

        # Find all raw files that you want to merge and concatonate them
        output_filepath = os.path.join(args.sup_parent_path, 'superparent_' + subpop)
        cat_cmd = 'cat {}/*{}* > {}'.format(args.input_path, subpop, output_filepath)
        subprocess.check_call(cat_cmd.split())
        logging.info('Superparent created for {} and written to\n{}\n'.format(subpop, args.sup_parent_path))

        # Run ustacks on the Superparent
        # ustacks -t file_type -f file_path [-d] [-r] [-o path] [-i id] [-m min_cov] [-M max_dist] [-p num_threads] [-R] [-H] [-h]

        ustacks_cmd = 'ustacks -t fastq -f {} -o {} -i {} -p {} -r -d'.format(args.sup_parent_path,
                                                                              args.stack_path,
                                                                              subpop,
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
        logging.info("Finished processing {} of {} subpopulations.".format(i ,len(args.subpops)))

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Constructs a "superparent" from all sequences in a specified subpopulation by concatonating them '
                    'and running on ustacks. This "superparent" can then be used to create a single catalogue using '
                    'cstacks.')

    parser.add_argument(
        "input_path",
        help="Location of input files")

    parser.add_argument(
        "sup_parent_path",
        help="Location to write superparent output")

    parser.add_argument(
        "stack_path",
        help="Location to write superparent output")

    parser.add_argument(
        "-subpops",
        default=None,
        nargs='+',
        help="List of string patterns denoting files that make up separate Subpopulations. ")

    parser.add_argument(
        "-p",
        help="Number of processors to run ustacks with.",
        metavar="processors",
        default=1)

    parser.add_argument(
        "-m",
        help="Minimum depth of coverage required to create a stack in ustacks (default 2).",
        metavar="mindepth",
        default=None)

    parser.add_argument(
        "-M",
        help="Maximum distance (in nucleotides) allowed between stacks in ustacks (default 2).",
        metavar="maxdist_between_stacks",
        default=None)

    parser.add_argument(
        "-N",
        help="Maximum distance allowed to align secondary reads to primary stacks in ustacks (default: M + 2).",
        metavar="maxdist_second2primary",
        default=None)

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