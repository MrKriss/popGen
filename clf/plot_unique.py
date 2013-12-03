#!/usr/bin/env python
# encoding: utf-8
"""
clf.plot_unique.py

@author:     cmusselle

@license:    license

@contact:    user_email
"""

import os, argparse, logging
import subprocess
import glob
from Bio import SeqIO
import numpy as np

from collections import Counter
import pandas as pd
import cPickle as pkl
import matplotlib.pyplot as plt

# PROCEDURE #
#-----------#
# Load in picked results
# Plot a log-log histogram for all frequencies

def plot_counters_hist(counters, bin_width=100, labels=None, log='xy', xlab="", ylab="",
                       title="", calc_rpc=False, **kwargs):
    ''' Construct a series of histogram plots from a list of Counter Dictionaries.

    counter dictionaries should be number of total reads per cluster size'''

    import matplotlib.pyplot as plt

    # Default plot options
#     if 'ms' not in kwargs:
#         kwargs['ms'] = 4.0
#     if 'marker' not in kwargs:
#         kwargs['marker'] = '.'
#     if 'mew' not in kwargs:
#         kwargs['mew'] = 0.0

    if type(counters) is not list and type(counters) is not tuple:
        counters = [counters]

    if labels is not None:
        assert len(labels) == len(counters), "Number of labels must match number of counters."

    for i, counter in enumerate(counters):

        # Compact to hist counter
        hist_counter = Counter()
        for seq, freq in counter.iteritems():
            bin = np.floor(float(freq)/bin_width)
            hist_counter[bin] += freq

        # Extract data
        bags = sorted(hist_counter.items())
        x_data, y_data = zip(*bags)
        x_data = np.array(x_data) * bin_width
        y_data = np.array(y_data)

        plt.step(x_data, y_data)

    # Plot formating
    ax = plt.gca()
    if 'x' in log:
        ax.set_xscale('log')
    if 'y' in log:
        ax.set_yscale('log')

    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    plt.legend(numpoints=1, markerscale=8)
    plt.show()

    return x_data, y_data

# Gather code in a main() function
def main(args, loglevel):

    # Setup Logging
    logging.basicConfig(format="%(levelname)s: %(message)s", level=loglevel)
    args_str = str([x for x in dir(args) if not x.startswith('_')])
    logging.debug('Argumnets passed:\n{}'.format(args_str))

    # Load data
    D = pkl.load(open(args.infile_path, 'r'))
    df = pd.read_csv(os.path.join(os.path.split(args.infile_path)[0], 'top10_read_count'))

    plot_hists(D, df)

def plot_hists(D, df, max_n=None, bins=10000):

    import matplotlib.ticker as ticker

    if not max_n:
        max_n = len(df)

    # Get individual with max reads and set histogram bins and max range
    x_max = max(D[df.loc[0,0]].values())

    print 'Bin width = ', np.floor(x_max / float(bins))

    for i, row in df.iloc[0:max_n, :].iterrows():

        key = row['Sample']
        x = key.split('_')
        y = x[-1].split('-')
        lab = x[2] + '-' + y[-1]
        # Plot
        plt.ion()
        plt.hist(D[key].values(), histtype='step', log=1, bins=bins, range=(0,x_max), label=lab, align='mid')

    plt.title('Unique Read Depths Distribution')
    plt.xlabel('Read Depth')
    plt.ylabel('Frequency')
    plt.legend()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_xlim(0, x_max)

    x = [0, 1, 10, 100, 1000, 10000, 100000]
    x = np.asarray(x) + 1
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x-1)))
    ax.xaxis.set_major_locator(ticker.FixedLocator(x))

def plot_hists2(D, df, max_n=None, bins=10000):

    import matplotlib.ticker as ticker

    if not max_n:
        max_n = len(df)

    # Get individual with max reads and set histogram bins and max range
    x_max = max(D[df.loc[0,0]].values())

    print 'Bin width = ', np.floor(x_max / float(bins))

    for i, row in df.iloc[0:max_n, :].iterrows():

        key = row['Sample']
        x = key.split('_')
        y = x[-1].split('-')
        lab = x[2] + '-' + y[-1]
        # Plot
        plt.ion()
        plt.hist(D[key].values(), histtype='step', log=1, bins=bins, range=(0,x_max), label=lab, align='mid')

    plt.title('Unique Read Depths Distribution')
    plt.xlabel('Read Depth')
    plt.ylabel('Frequency')
    plt.legend()
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_xlim(0, x_max)

    x = [0, 1, 10, 100, 1000, 10000, 100000]
    x = np.asarray(x) + 1
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x-1)))
    ax.xaxis.set_major_locator(ticker.FixedLocator(x))

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot distributions of read depth for top x samples with highest total read count.')

    parser.add_argument(
        "-i", dest="infile_path",
        required=True,
        help="Input path for file containing pickled results.")

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