'''
Created on 3 Dec 2012

@author: musselle
'''
import os 
import sys 

# Absolute Python path additions
sys.path.insert(0, '/path/to/directory/one')
sys.path.insert(0, '/path/to/directory/two')

# Relative path additions 
""" In a file called _mypath.py"""
import os, sys
thisdir = os.path.dirname(__file__)
libdir = os.path.join(thisdir, '../relative/path/to/lib/from/bin')

if libdir not in sys.path:
    sys.path.insert(0, libdir)
""" then put import _mypath at the top of each script """















if __name__ == '__main__':
    pass