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

from subprocess import call
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

def reads2fasta(infiles=None, filepattern=False, inpath='', 
                outfile='outfile.fasta', outpath = ''):
    '''Writes the reads (without the MID tag) to one large fasta file 
    for clustering'''
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, datapath=inpath)
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
    
    keys = [key[:6] for key in MIDdict.iterkeys()]
    keys.sort()
    
    # Make ouput directory if required
    outpath = os.path.join(datapath, outdir)
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
            

def process_MIDtag2(infiles=None, barcodes=None, filepattern=False, 
                   barcode_pattern=False, datapath='', barcode_path='',
                   outfile_postfix='-clean', outdir='cleaned_data'):
    ''' Goes through fastq files and corrects any errors in MIDtag 
    
    This version writes to fasta files, not compressed with bgzf. 
    
    '''
    import editdist as ed

    # Construct Tag dictionary
    MIDdict = make_MIDdict(infiles=barcodes, filepattern=barcode_pattern,
                           datapath=barcode_path)
            
    # Setup Record Cycler
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, datapath=datapath)
    
    keys = [key[:6] for key in MIDdict.iterkeys()]
    keys.sort()
    
    # Make ouput directory if required
    outpath = os.path.join(datapath, outdir)
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
        
        # Replace postfix
        if filename[-1] == 'bgzf':
            if filename[-2] == 'fastq':
                output_filename = '.'.join(filename[:-1])
            else:    
                output_filename = '.'.join(filename[:-1] + ['fastq'])
                
        outfile_path = os.path.join(outpath, output_filename)
        
#        output_filehdl = bgzf.BgzfWriter(output_filename, mode='wb')
        output_filehdl = open(outfile_path, mode='wb')
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
            
                
if __name__ == '__main__':
    
    
    
    
    
#===========================================================================
# Test for functions
#===========================================================================
#    datapath = '/space/musselle/datasets/gazellesAndZebras/testdata'
#    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
#    
#    files_test = 'testdata_1percent.bgzf'
    
#    files_L6 = 'lane6_NoIndex_L006_R1_001.fastq.bgzf'
#    files_L8 = 'lane8_NoIndex_L006_R1_001.fastq.bgzf'
#    
#    gz2bgzf(None, '*.gz', datapath = datapath + 'lane6/')
#    gz2bgzf(None, '*.gz', datapath = datapath + 'lane8/')
    
#    process_MIDtag2(infiles=files_test, barcodes='*.txt', 
#                   barcode_pattern=True, datapath=datapath, 
#                   barcode_path=barpath)
    
#===============================================================================
# Test Process MID Tags 2
#===============================================================================
    
    import glob
    from cluster import cluster_cdhit, summary
    
    LANE = '6'
    starting_dir = os.getcwd()

    # Set paths and file patterns 
    datapath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    os.chdir(datapath)
#    raw_files = glob.glob('*[0-9].fastq.bgzf')
#    raw_files.sort()
    
    raw_files = ['lane6_NoIndex_L006_R1_003.fastq.bgzf']
    
    outdir = 'L6_phredprop_filtered'

    # Update names and path
    filtered_files = []
    for name in raw_files:
        temp = name.split('.')
        temp[0] = temp[0] + '-pass'
        temp = '.'.join(temp) 
        filtered_files.append(temp)
    filtered_datapath = datapath + '/' + outdir
    
    cleaned_file_postfix = '-clean' 
    cleaned_outdir = 'cleaned_data'
    barcode_pattern = '*[' + LANE + '].txt'
    
    process_MIDtag2(infiles=filtered_files, barcodes=barcode_pattern,
               barcode_pattern=True, datapath=filtered_datapath, 
               barcode_path=barpath, outfile_postfix=cleaned_file_postfix, 
               outdir=cleaned_outdir)
    
    # Update names and path
    cleaned_files = []
    for name in filtered_files:
        
        temp = name.split('.')
                
        # Replace postfix
        if name.endswith('fastq.bgzf'):
            temp[0] = temp[0] + cleaned_file_postfix
            temp = '.'.join(temp[:-1])            
        else:
            temp = '.'.join(temp) 
        cleaned_files.append(temp) 
        
    cleaned_datapath = filtered_datapath + '/' + cleaned_outdir
    
    #===============================================================================
    # Cluster Data 
    #===============================================================================
    allreads_file = 'lane' + LANE + 'allreads-clean.fasta'
    reads2fasta(infiles=cleaned_files, datapath=cleaned_datapath, outfile=allreads_file)
    
    # Variables 
    c_thresh = 0.9
    n_filter = 8
    
    clustered_file = 'lane' + LANE + 'clustered_reads'
    cluster_cdhit(infile=allreads_file, outfile=clustered_file,
                  c_thresh=c_thresh, n_filter=n_filter)
    
    # Display Summary
    summary(clustered_file)
    
    
#===============================================================================
# Test Process MID Tags 1 
#===============================================================================
    
    import glob
    from cluster import cluster_cdhit, summary
    
    LANE = '6'
    starting_dir = os.getcwd()

    # Set paths and file patterns 
    datapath = '/space/musselle/datasets/gazellesAndZebras/lane' + LANE
    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    os.chdir(datapath)
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
    filtered_datapath = datapath + '/' + outdir
    
    cleaned_file_postfix = '-clean' 
    cleaned_outdir = 'cleaned_data'
    barcode_pattern = '*[' + LANE + '].txt'
    
    process_MIDtag(infiles=filtered_files, barcodes=barcode_pattern,
               barcode_pattern=True, datapath=filtered_datapath, 
               barcode_path=barpath, outfile_postfix=cleaned_file_postfix, 
               outdir=cleaned_outdir)
    
    # Update names and path
    cleaned_files = []
    for name in filtered_files:
        temp = name.split('.')
        temp[0] = temp[0] + cleaned_file_postfix
        temp = '.'.join(temp) 
        cleaned_files.append(temp) 
    cleaned_datapath = filtered_datapath + '/' + cleaned_outdir
    
    #===============================================================================
    # Cluster Data 
    #===============================================================================
    allreads_file = 'lane' + LANE + 'allreads-clean.fasta'
    reads2fasta(infiles=cleaned_files, datapath=cleaned_datapath, outfile=allreads_file)
    
    # Variables 
    c_thresh = 0.9
    n_filter = 8
    
    clustered_file = 'lane' + LANE + 'clustered_reads'
    cluster_cdhit(infile=allreads_file, outfile=clustered_file,
                  c_thresh=c_thresh, n_filter=n_filter)
    
    # Display Summary
    summary(clustered_file)