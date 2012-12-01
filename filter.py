'''
Created on 19 Nov 2012

@author: musselle
'''
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

from basicAnalysis import Cycler
from Bio import SeqIO, bgzf
import gzip
import time


def filter_reads(infiles=None, filetype='', datapath='', filterfunc=None, outdir='filtered_reads'):
    ''' Filter reads based on criteria 
    
    Default is to use Machine Specific read filter, specific to Casava 
    1.8 Illumina output format at current
    
    filterfunc must take in a sequence record object, and return a boolean
    
    '''   
    starting_dir = os.getcwd()

    # Define filter function if not given 
    if filterfunc is None:
        def fitlerfunc(rec):
            ''' Use machine specific filter         
            N = was not pickup up by machine filter i.e. passed
            Y = was flagged by machine filter i.e. fail 
            '''
            return rec.description.split()[1].split(':')[1] == 'N'
         
    # Have to create two separate generators to return passes and fails 
    # as copying a generator object is not possible.
    RecCycler1 = Cycler(infiles=infiles, filetype=filetype, datapath=datapath)
    RecCycler2 = Cycler(infiles=infiles, filetype=filetype, datapath=datapath)

    if not os.path.isdir('./' + outdir):
        os.mkdir(outdir)
    
    # timings
    toc = time.time()
    cum_t = 0
    for recordgen1 in RecCycler1.seqfilegen:
        # Both gens are identical at this point so iterating over one is the 
        # same length as the other
        
        # Load second record generator
        recordgen2 = RecCycler2.seqfilegen.next()

        print 'Processing {0}'.format(RecCycler1.curfilename) 
        
        # Construct file names
        name = RecCycler1.curfilename.split('.')  
        pass_filename = [name[0] + '-pass'] + name[1:] 
        pass_filename = '.'.join(pass_filename)
        fail_filename = [name[0] + '-fail']  + name[1:]
        fail_filename = '.'.join(fail_filename)
        name = '.'.join(name)
        
        # Setup Generators              
        passgen = (rec for rec in recordgen1 if filterfunc(rec))
        failgen = (rec for rec in recordgen2 if not filterfunc(rec))
        
        os.chdir(outdir)
        
        if name.endswith('.bgzf'):
            pass_filehdl = bgzf.BgzfWriter(pass_filename)
            fail_filehdl = bgzf.BgzfWriter(fail_filename)
        elif name.endswith('.fastq'):
            pass_filehdl = open(pass_filename, 'wb')
            fail_filehdl = open(fail_filename, 'wb')
        elif name.endswith('.gz'):
            pass_filehdl = gzip.open(pass_filename, 'wb')
            fail_filehdl = gzip.open(fail_filename, 'wb')
        else:
            print 'Input file format not supported'
            sys.exit()
        
        print 'Writing passes to file {0} ....'.format(pass_filename)
        numwritten = SeqIO.write( passgen , pass_filehdl , 'fastq')
        pass_filehdl.close()
        print '{0} records written'.format(numwritten)
        
        print 'Writing fails to file  {0} ....'.format(fail_filename)
        numwritten = SeqIO.write( failgen , fail_filehdl , 'fastq')
        fail_filehdl.close()
        print '{0} records written'.format(numwritten)
        
        os.chdir('..')
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished file {0} after {1}'.format(RecCycler1.curfilenum, 
                                time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
        
    total_t = time.time() - toc    
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                        time.gmtime(total_t)))
    os.chdir(starting_dir)

def setup_filter(target_dict):
    ''' Function to return a filter function defined using the given dictionary
    
    Keys are the variables to filter, and the values are the minimum thresholds 
    
    '''

    if len(target_dict) == 1:
        if 'phred' in target_dict:
            # Define filterfunc
            def f(rec):
                ''' filter function for phred '''
                return np.array(rec.letter_annotations['phred_quality']).mean() > target_dict['phred']
            return f
    
        elif 'propN' in target_dict:
            # Define filterfunc
            def f(rec):
                ''' filter function fro propN'''
                return float(rec.seq.count('N')) / len(rec.seq) < target_dict['propN']
            return f
        
    elif len(target_dict) == 2:
        if 'phred' in target_dict and 'propN' in target_dict:
            # Define filterfunc
            def f(rec):
                ''' filter function fro phred and propN'''
                A = np.array(rec.letter_annotations['phred_quality']).mean() > target_dict['phred']
                B = float(rec.seq.count('N')) / len(rec.seq) > target_dict['propN']
                return A and B
            return f
        else:
            raise Exception('Target variables not set as ''phred'' and ''propN''')
    else:
        raise Exception('Number of target values > 2')

if __name__ == '__main__':
    
    dataloc = '/space/musselle/datasets/gazellesAndZebras'
    files = 'testdata_1percent.bgzf'
    outdir = 'machinefiltertest'
    filter_reads(infiles=files, datapath=dataloc, outdir=outdir)
    
    outdir = 'propNfiltertest'
    f = setup_filter({'propN' : 0.1})
    filter_reads(infiles=files, datapath=dataloc, outdir=outdir, filterfunc=f)
    
    outdir = 'phredfiltertest'
    f = setup_filter({'phred' : 15})
    filter_reads(infiles=files, datapath=dataloc, outdir=outdir, filterfunc=f)
    
    outdir = 'phredpropN_filtertest'
    f = setup_filter({'phred' : 15, 'propN' : 0.05})
    filter_reads(infiles=files, datapath=dataloc, outdir=outdir, filterfunc=f)
    
    
    