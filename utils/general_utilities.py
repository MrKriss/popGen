'''
Created on Feb 28, 2013

General utility recipies across all projects

@author: chris
'''

import sys

def set_trace():
    from IPython.core.debugger import Pdb 
    Pdb(color_scheme='Linux').set_trace(sys._getframe().f_back)

def debugf(f, *args, **kwargs):
    from IPython.core.debugger import Pdb
    pdb = Pdb(color_scheme='Linux')
    return pdb.runcall(f, *args, **kwargs)
