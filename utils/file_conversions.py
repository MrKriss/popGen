'''
Created on Oct 15, 2012

@author: chris
'''
import os 
from Bio.bgzf import BgzfWriter
import gzip
from Bio import SeqIO, bgzf
import time
from utils import smartopen, Cycler, make_MIDdict

import numpy as np


def gz2bgzf(infiles=None, filepattern=False, datapath='', SQLindex=True):
    ''' Convert the list of files from .gz to .bgzf,
    And also produce an SQL index if needed. 
    
    infiles accepts str of file name of list of str for filenames. 
    If not specified will look at file type and glob the result to infiles. 
  
    '''
    if datapath:
        os.chdir(datapath)
  
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
        w = BgzfWriter(bgzfFileName, 'wb')
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
  
def makeSQLindex(infiles=None, filepattern=False, datapath=''):
    ''' Creates an SQL index out of either an uncompressed file or a compressed .bgzf file 
    
    if infiles is list, goes through all file names in list
    
    '''
    if datapath:
        os.chdir(datapath)
  
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
    
    handle = smartopen(filename)
    out_filename = filename.split('.')[0] + '.fasta'
    
    count = SeqIO.convert(handle, 'fastq', out_filename, 'fasta')
    
    print 'Converted {0} records to file\n{1}'.format(count, out_filename)

def reads2fasta(infiles=None, filepattern=False, datapath=''):
    '''Writes the reads (without the MID tag) to a fasta file for clustering'''
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, datapath=datapath)



    
    handle = smartopen(filename)
    out_filename = filename.split('.')[0] + '.fasta'
    
    count = SeqIO.convert(handle, 'fastq', out_filename, 'fasta')
    
    print 'Converted {0} records to file\n{1}'.format(count, out_filename)


def process_MIDtag(infiles=None, barcodes=None, filepattern=False, 
                   barcode_pattern=False, datapath='', barcode_path='',
                   outfile_postfix='-clean', outdir='cleaned_data'):
    ''' Goes through fastq files and corrects any errors in MIDtag 
    
    '''
    import editdist as ed

    # Construct Tag dictionary
    MIDdict = make_MIDdict(infiles=barcodes, filepattern=barcode_pattern,
                           datapath=barcode_path)
    # Setup Record Cycler
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, datapath=datapath)
    
    keys = sorted(MIDdict.keys())
    
    # Make ouput directory if required
    outpath = datapath + '/' + outdir
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
 
    def ok_reads_gen(recgen, keys):
        ''' Generator to yield reads that pass a quality check or are corrected 
        sucessfuly.
        '''
        ok_reads_gen.skipped_count = 0
        ok_reads_gen.corrected_count = 0
        
        for rec in recgen:
            recMID = str(rec.seq[:6])
            if recMID not in MIDdict:
                # Sequencing error in the tag. Work out nearest candidate.
                distvec = np.array([ed.distance(recMID, key[:6]) for key in keys]) 
                min_dist_candidates = [ keys[idx][:6] for idx in np.where(distvec == distvec.min())[0]]
                if len(min_dist_candidates) > 1:
                    # Muliple candidates. True MID is Ambiguous 
#                    print ('Multiple minimum distances. ' 
#                    'MID could not be resolved between\n{0}' 
#                    '  and \n{1}').format(recMID, min_dist_candidates)
#                    print 'Skipping read.'
                    ok_reads_gen.skipped_count += 1
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
                    ok_reads_gen.corrected_count += 1
                    
            yield rec
        
    # main loop
    total_numwritten = 0
    total_numskipped = 0
    total_numcorrected = 0
    toc = time.time()
    cum_t = 0
    for seqfile in RecCycler.seqfilegen:
            
        # Make reads generator
        ok_reads = ok_reads_gen(seqfile, keys)
        
        # Set output filename and handle
        filename = RecCycler.curfilename
        filename = filename.split('.')
        filename[0] = filename[0] + outfile_postfix
        filename = ''.join(filename)
        output_filename = filename
        output_filehdl = bgzf.BgzfWriter(output_filename, mode='w')
                   
        if os.getcwd() != outpath:
            os.chdir(outpath)
        
        numwritten = SeqIO.write(ok_reads, output_filehdl, 'fastq')
        print '{0} records written, of which \
        {1} were corrected'.format(numwritten, ok_reads.corrected_count)
        total_numwritten += numwritten
        total_numcorrected += ok_reads.corrected_count
        print '{0} records skipped'.format(ok_reads.skipped_count)
        total_numskipped += ok_reads.skipped_count
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished {0} after {1}'.format(filename, 
                        time.strftime('%H:%M:%S', time.gmtime(loop_t)))

    print 'Total records written: {0}'.format(numwritten)
    total_numwritten += numwritten
    print 'Total records skipped: {0}'.format(ok_reads.skipped_count)
    total_numskipped += ok_reads.skipped_count
    print 'Total of {0} tags corrected.'.format(total_numcorrected)
            
    total_t = time.time() - toc    
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                        time.gmtime(total_t)))
            
                
if __name__ == '__main__':
    
    datapath = '/space/musselle/datasets/gazellesAndZebras/testdata'
    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    
    files_test = 'testdata_1percent.bgzf'
    
#    files_L6 = 'lane6_NoIndex_L006_R1_001.fastq.bgzf'
#    files_L8 = 'lane8_NoIndex_L006_R1_001.fastq.bgzf'
#    
    process_MIDtag(infiles=files_test, barcodes='*.txt', 
                   barcode_pattern=True, datapath=datapath, 
                   barcode_path=barpath)
    
    
#    gz2bgzf(None, '*.gz', datapath = datapath + 'lane6/')
#    gz2bgzf(None, '*.gz', datapath = datapath + 'lane8/')




