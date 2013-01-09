'''
Created on 19 Nov 2012

@author: musselle
'''
import os
import sys
import time
import gzip
from subprocess import call

import numpy as np
from Bio import SeqIO, bgzf

from utils import smartopen, Cycler, make_MIDdict

def filter_reads(infiles=None, filepattern='', inpath='', filterfunc=None, 
                 outdir='filtered_reads'):
    ''' Filter reads based on criteria 
    
    Default is to use Machine Specific read filter, specific to Casava 
    1.8 Illumina output format at current
    
    filterfunc must take in a sequence record object, and return a boolean
    
    '''   
    # Path variables
    starting_dir = os.getcwd()    
    outpath = os.path.join(starting_dir, outdir) 
    if not os.path.isdir(outpath):
        os.mkdir(outdir)

    # Define filter function if not given 
    if filterfunc is None:
        def filterfunc(rec):
            ''' Use machine specific filter         
            N = was not pickup up by machine filter i.e. passed
            Y = was flagged by machine filter i.e. fail 
            '''
            return rec.description.split()[1].split(':')[1] == 'N'
         
    # Have to create two separate generators to return passes and fails 
    # as copying a generator object is not possible.
    RecCycler1 = Cycler(infiles=infiles, filepattern=filepattern, inpath=inpath)
    RecCycler2 = Cycler(infiles=infiles, filepattern=filepattern, inpath=inpath)
    
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
        pass_filename = os.path.join(outpath, '.'.join(pass_filename))
        fail_filename = [name[0] + '-fail']  + name[1:]
        fail_filename = os.path.join(outpath, '.'.join(fail_filename))
        name = '.'.join(name)
        
        # Setup Generators              
        passgen = (rec for rec in recordgen1 if filterfunc(rec))
        failgen = (rec for rec in recordgen2 if not filterfunc(rec))
        
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
        
        print 'Writing passes to \n{0} ....'.format(pass_filename)
        numwritten = SeqIO.write( passgen , pass_filehdl , 'fastq')
        pass_filehdl.close()
        print '{0} records written'.format(numwritten)
        
        print 'Writing fails to \n{0} ....'.format(fail_filename)
        numwritten = SeqIO.write( failgen , fail_filehdl , 'fastq')
        fail_filehdl.close()
        print '{0} records written'.format(numwritten)
        
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
                B = float(rec.seq.count('N')) / len(rec.seq) < target_dict['propN']
                return A and B
            return f
        else:
            raise Exception('Target variables not set as ''phred'' and ''propN''')
    else:
        raise Exception('Number of target values > 2')

def process_MIDtag(infiles=None, barcodes=None, filepattern=False, 
                   barcode_pattern=False, inpath='', barcode_path='',
                   outfile_postfix='-clean', outdir='cleaned_data'):
    ''' Goes through fastq files and corrects any errors in MIDtag 
    
    '''
    import editdist as ed

    # Construct Tag dictionary
    MIDdict = make_MIDdict(infiles=barcodes, filepattern=barcode_pattern,
                           inpath=barcode_path)
    # Setup Record Cycler
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, inpath=inpath)
    
    keys = [key[:6] for key in MIDdict.iterkeys()]
    keys.sort()
    
    # Make ouput directory if required
    outpath = os.path.join(inpath, outdir)
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
 
    class Read_Corrector_Class():
    
        def __init__(self):
            ''' Class for a Generator to yield reads that pass a quality check or are corrected 
            
            Keeps track of number of skipped and corrected reads as class attributes.
            '''
            self.skipped_count = 0
            self.corrected_count = 0
    
        def ok_reads_gen(self, recgen, keys):
            
            for rec in recgen:
                recMID = str(rec.seq[:6])
                if recMID not in keys:
                    # Sequencing error in the tag. Work out nearest candidate.
                    distvec = np.array([ed.distance(recMID, key) for key in keys]) 
                    min_dist_candidates = [keys[idx] for idx in np.where(distvec == distvec.min())[0]]
                    if len(min_dist_candidates) > 1:
                        # Muliple candidates. True MID is Ambiguous 
    #                    print ('Multiple minimum distances. ' 
    #                    'MID could not be resolved between\n{0}' 
    #                    '  and \n{1}').format(recMID, min_dist_candidates)
    #                    print 'Skipping read.'
                        self.skipped_count += 1
                        continue
                    elif len(min_dist_candidates) == 1:
                        # Correct the erroneous tag with the candidate.
                        # Letter annotations must be removed before editing seq.
                        temp_var = rec.letter_annotations
                        rec.letter_annotations = {}
                        # Change seq to mutableseq
                        rec.seq = rec.seq.tomutable()
                        rec.seq[:6] = min_dist_candidates[0]
                        rec.seq = rec.seq.toseq()
                        rec.letter_annotations.update(temp_var)
                        self.corrected_count += 1
                        
                yield rec
            
    # main loop
    total_numwritten = 0
    total_numskipped = 0
    total_numcorrected = 0
    toc = time.time()
    cum_t = 0
    for seqfile in RecCycler.seqfilegen:
            
        # Make reads class for generator
        ReadCorrector = Read_Corrector_Class()
        
        # Set output filename and handle
        filename = RecCycler.curfilename
        filename = filename.split('.')
        filename[0] = filename[0] + outfile_postfix
        filename = '.'.join(filename)
        
        outfile_path = os.path.join(outpath, filename)
        
        output_filehdl = bgzf.BgzfWriter(outfile_path, mode='wb')
        
        numwritten = SeqIO.write(ReadCorrector.ok_reads_gen(seqfile, keys), 
                                 output_filehdl, 'fastq')
        output_filehdl.flush()        
        output_filehdl.close()        

        print ('{0} records written, of which ' 
        '{1} were corrected').format(numwritten, ReadCorrector.corrected_count)
        total_numwritten += numwritten
        total_numcorrected += ReadCorrector.corrected_count
        print '{0} records skipped'.format(ReadCorrector.skipped_count)
        total_numskipped += ReadCorrector.skipped_count
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(filename, 
                        time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    print 'Total records written: {0}'.format(total_numwritten)
    print 'Total records skipped: {0}'.format(total_numskipped)
    print 'Total of {0} tags corrected.'.format(total_numcorrected)
            
    total_t = time.time() - toc    
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                        time.gmtime(total_t)))


def gz2bgzf(infiles=None, filepattern=False, inpath='', SQLindex=True):
    ''' Convert the list of files from .gz to .bgzf,
    And also produce an SQL index if needed. 
    
    infiles accepts str of file name of list of str for filenames. 
    If not specified will look at file type and glob the result to infiles. 
  
    '''
    if inpath:
        os.chdir(inpath)
  
    # Handle multiple types of input for infiles
    assert infiles is not None, 'No files listed or file pattern specified.'         
    if filepattern:
        # Fetch files by file types using glob
        import glob 
        infiles = glob.glob(infiles)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]
  
    start_time = time.time() 
    
    for filename in infiles: 
        toc = time.time()
        
        # Checks for type of input
        if filename.split('.')[-1] == '.gz':
            # Drop .gz and append .bgzf
            f2read = gzip.open(filename)
            bgzfFileName = '.'.join(filename.split('.')[:-1]) + '.bgzf'  
        elif filename.split('.')[-1] == '.fastq':
            f2read = open(filename)
            # Append .bgzf
            bgzfFileName = filename + '.bgzf'

        print "Producing BGZF output from {0}...".format(filename)
        w = bgzf.BgzfWriter(bgzfFileName, 'wb')
        while True:
            data = f2read.read(65536)
            w.write(data)
            if not data:
                break
        w.close()
        print '{0} written successfully'.format(bgzfFileName)

        conv_t = time.time() - toc 
        print 'Finished converting {0}\n after {1}\n'.format(filename, time.strftime('%H:%M:%S', time.gmtime(conv_t)))
        
        if SQLindex == True:
            makeSQLindex(bgzfFileName)
      
    total_t = time.time() - start_time
    print 'Finished all processing {0} files in {1}'.format(len(infiles), time.strftime('%H:%M:%S', time.gmtime(total_t)))
  
def makeSQLindex(infiles=None, filepattern=False, inpath=''):
    ''' Creates an SQL index out of either an uncompressed file or a compressed .bgzf file 
    
    if infiles is list, goes through all file names in list
    
    '''
    if inpath:
        os.chdir(inpath)
  
    # Handle multiple types of input for infiles
    assert infiles is not None, 'No files listed or file pattern specified.'         
    if filepattern:
        # Fetch files by file types using glob
        import glob 
        infiles = glob.glob(infiles)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]

    for filename in infiles: 
        tak = time.time()
        print 'Writing SQL index file for {0} ...'.format(filename)
        idx_filename = filename.split('.')[0] + '.idx'
        SeqIO.index_db(idx_filename, filename , 'fastq')
        print '{0} written successfully'.format(idx_filename)
        idx_t = time.time() - tak
        print 'Finished Indexing to {0}\n after {1}\n'.format(idx_filename, time.strftime('%H:%M:%S', time.gmtime(idx_t)))

def file2fasta(filename):
    ''' Convert fastq file to fasta file '''
    
    handle = smartopen(filename)
    out_filename = filename.split('.')[0] + '.fasta'
    
    count = SeqIO.convert(handle, 'fastq', out_filename, 'fasta')
    
    print 'Converted {0} records to file\n{1}'.format(count, out_filename)

def reads2fasta(infiles=None, filepattern=False, inpath='', 
                outfile='outfile.fasta', outpath = ''):
    '''Writes the reads (without the MID tag) to one large fasta file 
    for clustering'''
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, inpath=inpath)
    count = 0
    outfile_part_list = []

    print ('Removing MID tags and converting {0} files to'
           ' fasta format').format(RecCycler.numfiles)

    # Generator to trip off tag.
    for seqfilegen in RecCycler.seqfilegen:
        
        read_gen = (rec[12:] for rec in seqfilegen)
        # File name 
        outfile_part = 'output_part' + str(count) + '.fasta'
        outfile_part_list.append(outfile_part)
        count += 1
        with open(os.path.join(outpath, outfile_part), 'wb') as f:
            write_count = SeqIO.write(read_gen, f, 'fasta')
            print 'Wrote {0} records to file\n{1}'.format(count, outfile_part)
    
    # Combine output parts into one big file
    cmd = ['cat'] + outfile_part_list 
    with open(os.path.join(outpath, outfile), 'wb') as f:
        print 'Running "{0}" and saving to\n{1}'.format(cmd, os.path.join(outpath, outfile))
        call(cmd, stdout=f) 
        print 'Done'         

if __name__ == '__main__':
    
#===============================================================================
#  Test for data filters
#===============================================================================
    dataloc = '/space/musselle/datasets/gazellesAndZebras'
    files = 'testdata_1percent.bgzf'
    outdir = 'machinefiltertest'
    filter_reads(infiles=files, inpath=dataloc, outdir=outdir)
    
    outdir = 'propNfiltertest'
    f = setup_filter({'propN' : 0.1})
    filter_reads(infiles=files, inpath=dataloc, outdir=outdir, filterfunc=f)
    
    outdir = 'phredfiltertest'
    f = setup_filter({'phred' : 15})
    filter_reads(infiles=files, inpath=dataloc, outdir=outdir, filterfunc=f)
    
    outdir = 'phredpropN_filtertest'
    f = setup_filter({'phred' : 15, 'propN' : 0.05})
    filter_reads(infiles=files, inpath=dataloc, outdir=outdir, filterfunc=f)
    
#===============================================================================
# Test Process MID Tags 1 
#===============================================================================
    
    import glob
    from cluster import cluster_cdhit, summary
    
    LANE = '6'
    starting_dir = os.getcwd()

    # Set paths and file patterns 
    inpath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    os.chdir(inpath)
#    raw_files = glob.glob('*[0-9].fastq.bgzf')
#    raw_files.sort()
    
    raw_files = ['lane6_NoIndex_L006_R1_002.fastq.bgzf']
    
    outdir = 'L6_phredprop_filtered'

    # Update names and path
    filtered_files = []
    for name in raw_files:
        temp = name.split('.')
        temp[0] = temp[0] + '-pass'
        temp = '.'.join(temp) 
        filtered_files.append(temp)
    filtered_inpath = inpath + '/' + outdir
    
    cleaned_file_postfix = '-clean' 
    cleaned_outdir = 'cleaned_data'
    barcode_pattern = '*[' + LANE + '].txt'
    
    process_MIDtag(infiles=filtered_files, barcodes=barcode_pattern,
               barcode_pattern=True, inpath=filtered_inpath, 
               barcode_path=barpath, outfile_postfix=cleaned_file_postfix, 
               outdir=cleaned_outdir)
    
    # Update names and path
    cleaned_files = []
    for name in filtered_files:
        temp = name.split('.')
        temp[0] = temp[0] + cleaned_file_postfix
        temp = '.'.join(temp) 
        cleaned_files.append(temp) 
    cleaned_inpath = filtered_inpath + '/' + cleaned_outdir
    
    #===============================================================================
    # Cluster Data 
    #===============================================================================
    allreads_file = 'lane' + LANE + 'allreads-clean.fasta'
    reads2fasta(infiles=cleaned_files, inpath=cleaned_inpath, outfile=allreads_file)
    
    # Variables 
    c_thresh = 0.9
    n_filter = 8
    
    clustered_file = 'lane' + LANE + 'clustered_reads'
    cluster_cdhit(infile=allreads_file, outfile=clustered_file, c_thresh=c_thresh, n_filter=n_filter)
    
    # Display Summary
    summary(clustered_file)
    