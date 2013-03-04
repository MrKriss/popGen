'''
Created on 28 Feb 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 

from plot_utils import cluster_summary_plot, hist_counter

infileL6 = '/space/musselle/datasets/gazelles-zebras/clusters/L6clustered_reads'
infileL8 = '/space/musselle/datasets/gazelles-zebras/clusters/L8clustered_reads'

outL6 = cluster_summary_plot(infileL6, plot_hist = 0)
outL8 = cluster_summary_plot(infileL8, plot_hist = 0)

hist_counter(outL6[1], bins=5000, range =(1,10000),label='Lane 6')
hist_counter(outL8[1], bins=5000, range =(1,10000),label='Lane 8')

