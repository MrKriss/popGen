'''
Created on 4 Feb 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

import re

'''
The good news is that as of CASAVA version 1.8, Illumina sequencers will
produce FASTQ files using the standard Sanger encoding.
'''


# Function to check what type of fatsq format is used 

def check_fastq_encoding(filename):
    ''' Checks file to see which fastq format is used for the encoding
    
    Inputs: fastq file name (.fastq extention not required)
    
    Outputs: str
    - "fastq" means Sanger style FASTQ files using PHRED scores and an ASCII
    offset of 33 (e.g. from the NCBI Short Read Archive and Illumina 1.8+).
    These can potentially hold PHRED scores from 0 to 93.
    
     - "fastq-solexa" means old Solexa (and also very early Illumina) style FASTQ
    files, using Solexa scores with an ASCII offset 64. These can hold Solexa
    scores from -5 to 62.
    
     - "fastq-illumina" means newer Illumina 1.3 to 1.7 style FASTQ files, using
    PHRED scores but with an ASCII offset 64, allowing PHRED scores from 0
    to 62.
        
    Patterns: 
       
    Sanger or Illumina1.8 are only formats to contain characters [!-:]+  
    Solexa is only format to contain characters [;-?]+ other than Sanger
    
    Illumina1.3 is only format to contain characters [@-A]+ other than Sanger and Solexa   
    Illumina1.5 = [B-h]+ not [!-A]
       
       
     TODO:
     Currently has an error in the logic. Must fix. Needs to sample many sequences
     not just first character.
       
       
    '''
    # compile regex expression filters
    sanger_only_scores = re.compile(r'[!-:]+')
    solexa_scores = re.compile(r'[;-?]+')

    scan_limit = 1000
    count = 0

    # Open file 
    with open(filename, 'rb') as f:
        for i, line in enumerate(f):
            if not (i+1) % 4 : # stop at every 4th line
               
                if sanger_only_scores.match(line):
                    # Its in Sanger Format
                    print 'Quality Encoding detected: Standard Sanger format (Illumina 1.8+)'
                    ans = 'fastq'
                    break                  
                elif solexa_scores.match(line):
                    # Its in Solexa Format
                    print 'Quality Encoding detected: Solexa format'
                    ans = 'fastq-solexa'
                    break
                else:
                    count += 1
                    if count > scan_limit:
                        print 'Most likely Encoding after Sampling %s reads:' % scan_limit
                        print 'Illumina 1.3 - 1.7 format'
                        ans = 'fastq-illumina'
                        break
    return ans

def check_seqid_format(filename):
    ''' Finds and returns the sequence id format used 
    
    Possibilities:
      
    ## Illumina ##
    
    @HWUSI-EAS100R:6:73:941:1973#0/1
    
    HWUSI-EAS100R    the unique instrument name
    6       flowcell lane
    73      tile number within the flowcell lane
    941     'x'-coordinate of the cluster within the tile
    1973    'y'-coordinate of the cluster within the tile
    #0      index number for a multiplexed sample (0 for no indexing)
    /1      the member of a pair, /1 or /2 (paired-end or mate-pair reads only)  
    
    Versions of the Illumina pipeline since 1.4 appear to use #NNNNNN instead 
    of #0 for the multiplex ID, where NNNNNN is the sequence of the multiplex tag.
    
    ## Casava 1.8 ##
    
    @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
    
    EAS139   the unique instrument name
    136      the run id
    FC706VJ  the flowcell id
    2        flowcell lane
    2104     tile number within the flowcell lane
    15343    'x'-coordinate of the cluster within the tile
    197393   'y'-coordinate of the cluster within the tile
    1        the member of a pair, 1 or 2 (paired-end or mate-pair reads only)
    Y        Y if the read fails filter (read is bad), N otherwise
    18       0 when none of the control bits are on, otherwise it is an even number
    ATCACG   index sequence

    ## NCBI SRA ##
    
    @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
    GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
    +SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
    IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
    
    SRR001666.1                         NCBI-assigned identifier
    071112_SLXA-EAS1_s_7:5:1:817:345    as per the illumina format above
    length=36                           length of sequence
    
    '''

    # Compile regexs
    illumina_pattern = re.compile(r'@[\w-]+:\d+:\d+:\d+:\d+')
    casava_pattern = re.compile(r'@[\w-]+:\d+:\w+:\d+:\d+:\d+:\d+ [12]:[YN]:\d+:[ATCG]*')
    ncbi_sra_pattern = re.compile(r'@\w+\.\d ([\w-]+:\d+:\d+:\d+:\d+|\d+_\d+_\d+_\d+) length=\d+')

    # Open file 
    with open(filename, 'rb') as f:
        seqid = f.readline()
      
        if illumina_pattern.match(seqid):
            # Its in Illumina Format
            print 'Seq ID: Illumina format'
            ans = 'illumina'                  
        elif casava_pattern.match(seqid):
            # Its in Casava 1.8 Format
            print 'Seq ID: Casava 1.8 format'
            ans = 'casava'
        elif ncbi_sra_pattern.match(seqid):
            # Its in NCBI SRA Format
            print 'Seq ID: NCBI SRA format'
            ans = 'ncbi_sra'  
        else:
            print 'Seq ID not recognised'
            ans = None
    return ans

if __name__ == '__main__':
    
    # Files to Test 
    files = ['/space/musselle/datasets/sticklebacks/SRR034310.fastq',
             '/space/musselle/datasets/gazellesAndZebras/lane6/lane6_NoIndex_L006_R1_004.fastq']

    for f in files:
        check_fastq_encoding(f)
        check_seqid_format(f)




