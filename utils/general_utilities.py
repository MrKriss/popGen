'''
Created on Feb 28, 2013

General utility recipies across all projects

@author: chris
'''

import sys
import os

def set_trace():
    from IPython.core.debugger import Pdb 
    Pdb(color_scheme='Linux').set_trace(sys._getframe().f_back)

def debugf(f, *args, **kwargs):
    from IPython.core.debugger import Pdb
    pdb = Pdb(color_scheme='Linux')
    return pdb.runcall(f, *args, **kwargs)

def print_attr(obj):
    """ Print out all non private methods and attributes of an abject """
    for attr in dir(obj):
        if attr.startswith('__'):
            continue
        else:
            print attr
            print getattr(obj, attr)
            print '\n'
            
def get_data_prefix():
    '''Return directory prefix where data is stored'''    
    import socket

    if socket.gethostname() == 'yildun':
        prefix = '/space/musselle/data/RAD-seq'
    elif socket.gethostname() == 'luca':
        prefix = '/home/musselle/san/data/RAD-seq'
    elif socket.gethostname() == 'gg-pc6':
        prefix = '/home/musselle/data/RAD-seq'

    return prefix


import string
def search_file(filename, search_path, pathsep=os.pathsep):
    """ Given a search path, find file with requested name """
    for path in string.split(search_path, pathsep):
        candidate = os.path.join(path, filename)
        if os.path.exists(candidate): return os.path.abspath(candidate)
    return None