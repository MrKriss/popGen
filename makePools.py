'''
Created on 24 Oct 2012

@author: musselle
'''

import os
from Bio import SeqRecord

def make_MIDdict(infiles = None, filepattern = '', datapath = ''):
    ''' Function to load in MIDs from a list of files into a dictionary '''
    
    if datapath:
        os.chdir(datapath)
  
    # Handle multiple types of input
    if not infiles:
        # Fetch files by file types
        assert filepattern, 'No files listed and No file type specified.'
        import glob
        infiles = glob.glob(filepattern)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]

    # Run through files and store barcodes in a Dictionary object.
    # keys are starting tags (MID (6 BP) + cutsite (6BP))
    tags = {}
    
    for filename in infiles:
        with open(filename, 'rb') as f:
            for line in f:
                elem = line.split()
                tags[elem[0]] = elem[1] 
    return tags


                 
def distance():
    ''' Distance function to compare a sequence read to the list of barcode records'''
    pass


if __name__ == '__main__':
    
    
    datapath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    tags = loadIDs(filepattern = '*.txt', datapath = datapath)
    
    