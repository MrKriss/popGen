'''
Created on 9 Jan 2013

@author: musselle
'''

# Relative path additions 
""" In a file called _mypath.py"""

import os, sys

thisdir = os.path.dirname(__file__)

utildir = os.path.join(thisdir, '../utils')
runscriptdir = os.path.join(thisdir, '../runscripts')

if utildir not in sys.path:
    sys.path.insert(0, utildir)

if runscriptdir not in sys.path:
    sys.path.insert(0, runscriptdir)

""" then put import _mypath at the top of each script """