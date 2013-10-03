#!/usr/bin/env python
# encoding: utf-8
"""
clf.batch_ustacks -- Run STACKS ustack program over a batch of files 

@author:     cmusselle
              
@license:    license

@contact:    user_email

"""


# import modules used here -- sys is a very standard one
import os, sys, argparse, logging
import subprocess


# Gather code in a main() function
def main(args, loglevel):
    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)





    # TODO Replace this with your actual code.





    print "Hello there."
    logging.info("You passed an argument.")
    logging.debug("Your Argument: %s" % args.argument)


# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Wrapper for executing the ustacks program of the STACKS pipeline over multiple files. '
                    'Parameters are the same')

    # TODO Specify your real parameters here.
    parser.add_argument(
        "-i",
        help="Input file glob to parse to the ",
        metavar="ARG",
        nargs='+')

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



