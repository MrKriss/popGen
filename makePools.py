'''
Created on 24 Oct 2012

@author: musselle
'''

import os
from Bio import SeqRecord

def loadIDs(inFiles = None, fileType = '', dataPath = ''):
    ''' Function to load in MIDs from a list of files into a dictionary '''
    
    if dataPath:
        os.chdir(dataPath)
  
    # Handle multiple types of input
    if not inFiles:
        # Fetch files by file types
        assert fileType, 'No files listed and No file type specified.'
        import glob
        inFiles = glob.glob(fileType)
    elif type(inFiles) == str:
        # Convert to list
        inFiles = [inFiles]

    # Run through files and store barcodes in a Dictionary object.
    # keys are starting tags (MID (6 BP) + cutsite (6BP))
    tags = {}
    
    for fileName in inFiles:
        with open(fileName, 'rb') as f:
            for line in f:
                elem = line.split()
                tags[elem[0]] = elem[1] 
                
    return tags

                 
def distance():
    ''' Distance function to compare a sequence read to the list of barcode records'''
    pass


if __name__ == '__main__':
    
    
    dataPath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    tags = loadIDs(fileType = '*.txt', dataPath = dataPath)
    
    