'''
Created on 19 Mar 2013

@author: musselle
'''

import os
import sys
import time
import gzip
from subprocess import call
from collections import Counter
import glob

import numpy as np
from Bio import SeqIO, bgzf
import editdist as ed

from utils import Cycler

#===============================================================================
# Cluster IO file 

#===============================================================================
# Possibly usefull
#===============================================================================

def get_desc(self, cluster, items):
        """ Scan back over cluster and get all details """
        
        # Input checks  
        handle = input_check(self.handle)  

        # Store flag info
#         do_rep = 'rep' in items
#         do_phred = 'phred' in items
#         do_seq = 'seq' in items
#         do_dist = 'dist' in items

        # Update info 
        cluster.idx_file_path = self.idx_file_path

        with handle as cluster_file:
             
            # Load cluster 
            cluster_file.seek(cluster.start_loc)
            
            first = True
            for line in cluster_file:              
                line = line.strip() # Remove white space at start and end
                
                if line.startswith('>'):
                    # This is start of new cluster
                    
                    # yield this cluster
                    if first:
                        cluster.id = line.split()[1]
                        first = False
                    else: 
                        break
                    
                elif line.endswith('*'): 
                    # This is the representative sequence for the cluster
                    cluster.rep_seq_id = line.split()[2].strip('>.')
                    
#                     # Fill in info
#                     if do_seq or do_rep:
#                         if os.getcwd() != self.idx_file_dir:
#                             os.chdir(self.idx_file_dir)
#                         cluster.rep_seq = self.lookup_db[cluster.rep_seq_id].seq.tostring()
#                     if do_phred or do_rep:
#                         if os.getcwd() != self.idx_file_dir:
#                             os.chdir(self.idx_file_dir)
#                         cluster.rep_phred = np.array(self.lookup_db[cluster.rep_seq_id].letter_annotations['phred_quality'])
                else: 
                    
                    line_parts = line.split()
                    
                    next_desc = line_parts[2].strip('>.')
                    cluster.members_id.append(next_desc)
#                     if do_seq:
#                         if os.getcwd() != self.idx_file_dir:
#                             os.chdir(self.idx_file_dir)
#                         cluster.members_seq.append(self.lookup_db[next_desc].seq.tostring())
#                     if do_phred:
#                         if os.getcwd() != self.idx_file_dir:
#                             os.chdir(self.idx_file_dir)
#                         cluster.members_phred.append(np.array(self.lookup_db[next_desc].letter_annotations['phred_quality']))
#                     if do_dist: 
#                         similarity = line_parts[4].strip('+/%')
#                         seq_len = line_parts[1].strip('nt,')
#                         cluster.edit_dists.append(percentage2mismatch(100 - float(similarity), seq_len))
                    
            if os.getcwd() != self.output_dirpath:
                os.chdir(self.output_dirpath)
                    
            return cluster    

def filter_clusters(cluster_filepath, idx_filepath, size_range, output_dirpath):
    """ Writes a subset of cluster sizes to FastQ files 
    
    The representative sequence is the first sequence record written.
    
    """
    starting_dir = os.getcwd()
    idx_dir = os.path.split(idx_filepath)[0]
    
    # Check and create directory 
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
        
    cluster_gen = parse(cluster_filepath, idx_filepath)
    seqrec_lookup = SeqIO.index_db(idx_filepath)
    
    size_counter = Counter()
    
    for cluster in cluster_gen:
        # Check if cluster size is within defined range
        if cluster.size >= size_range[0] and cluster.size < size_range[1]:
            
            size_counter[cluster.size] += 1
            
            # Get the sequence records for the cluster 
            seqs = []
            if os.getcwd() != idx_dir:
                os.chdir(idx_dir)
            seqs.append(seqrec_lookup[cluster.rep_seq_id])
            for member in cluster.members_id:
                seqs.append(seqrec_lookup[member])
            
            if os.getcwd() != output_dirpath:
                os.chdir(output_dirpath)
            # Write cluster to a file 
            fname = "clustersize{0}-No{1}.fastq".format(str(cluster.size), str(size_counter[cluster.size]))
            output_handle = open(fname, "wb")
            SeqIO.write(seqs, output_handle, "fastq")
        
def filter_clusters2(cluster_filepath, idx_filepath, size_range, output_dirpath):
    """ Writes a subset of cluster sizes to FastQ files 
    
    The representative sequence is the first sequence record written.
    
    make the sequence record instead of passing it.
    
    """
    
    starting_dir = os.getcwd()
    idx_dir = os.path.split(idx_filepath)[0]
    
    # Check and create directory 
    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath)
        
    cluster_gen = parse(cluster_filepath, idx_filepath)
    seqrec_lookup = SeqIO.index_db(idx_filepath)
    
    size_counter = Counter()
    
    for cluster in cluster_gen:
        # Check if cluster size is within defined range
        if cluster.size > size_range[0] and cluster.size < size_range[1]:
            
            size_counter[cluster.size] += 1
            
            # Get the sequence records for the cluster 
            if os.getcwd() != idx_dir:
                os.chdir(idx_dir)
            # Representative sequence first 
            seqrecord = seqrec_lookup[cluster.rep_seq_id]
            
            if os.getcwd() != output_dirpath:
                os.chdir(output_dirpath)
            # Write cluster to a file 
            fname = "clustersize{0}-No{1}.fastq".format(str(cluster.size), str(size_counter[cluster.size]))
            
            if os.path.isfile(fname):
                output_handle = open(fname, "wb")
                output_handle.close()
            
            output_handle = open(fname, "a")
            SeqIO.write(seqrecord, output_handle, "fastq")
            
            for member in cluster.members_id :
                
                if os.getcwd() != idx_dir:
                    os.chdir(idx_dir)
                # Representative sequence first 
                seqrecord = seqrec_lookup[member]
            
                if os.getcwd() != output_dirpath:
                    os.chdir(output_dirpath)
                # Write sequence record to file 
                SeqIO.write(seqrecord, output_handle, "fastq")
                
    if os.getcwd() != starting_dir: 
        os.chdir(starting_dir)


#===============================================================================


def filter_reads(infiles=None, filepattern='', data_inpath='', filterfunc=None, 
                 outdir='filtered_reads', keepfails=False ):
    ''' Filter reads based on criteria 
    
    Default is to use Machine Specific read filter, specific to Casava 
    1.8 Illumina output format at current
    
    filterfunc must take in a sequence record object, and return a boolean
    
    keepfails - if true saves all records that failed the filter too.
    
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
    RecCycler1 = Cycler(infiles=infiles, filepattern=filepattern, data_inpath=data_inpath)
    RecCycler2 = Cycler(infiles=infiles, filepattern=filepattern, data_inpath=data_inpath)
    
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

def filter_reads_pipeline(infiles=None, filepattern='', data_inpath='', filterfuncs=None, 
                          outdir='filtered_reads', log_fails=False):
    
    ''' Filter reads based on criteria specified by a sequence of functions given
     in filterfuncs.
    
    Default if filterfuncs=None is to use Machine Specific read filter, specific to Casava 
    1.8 Illumina output format at current.
    
    indirvidual functions in filterfunc must take in a sequence record object, 
    and return a boolean.
    
    differs from filter_reads in that this does not save records that fail a filter,
    only those that pass, but if log_fails = True, will log all reads that fail 
    any of the filter functions given, along with the filter number it failed at:
    logfile example:
    
    (rec.id) <TAB> (enumeration of filter failed from list given)  
    e.g.
    HWI-ST0747:233:C0RH3ACXX:6:2109:3962:10857    2
    HWI-ST0747:233:C0RH3ACXX:6:1159:3231:10224    0
    ...
    
    
    Final output file is written in same file format as the input file. 
    
    Supported formats: .bgzf, .fastq and .gz
    
    '''
    #===============================================================================
    # Setup variables and utility functions 
    #===============================================================================
       
    # Path variables
    starting_dir = os.getcwd()
    if data_inpath:
        outpath = os.path.join(data_inpath, outdir) 
    else:        
        outpath = os.path.join(starting_dir, outdir) 
    
    if not os.path.isdir(outpath):
        if data_inpath:
            os.chdir(data_inpath)
        os.mkdir(outdir)
        os.chdir(starting_dir)
    
    # Setup Counters
    c = Counter()
    read_count = [0]
        
    # Define filter function loop
    if log_fails: 
        # Set up log file
        print '\nSaving Failed reads to log file: fails.log'
        if os.path.isfile(os.path.join(starting_dir, outdir, 'fails.log')):
            print 'Overwriting previous Log file'
            logfile = open(os.path.join(outpath, 'fails.log'), 'wb')
        else:
            logfile = open(os.path.join(outpath, 'fails.log'), 'wb')
               
    def make_pass_gen(recordgen, filterfuncs):
        ''' An easier to read implimentation of a generator to only yield 
        records that pass all filter functions specified'''
    
        for rec in recordgen:
            for n, func in enumerate(filterfuncs):
                # Run filter functions on record       
                if func(rec) == False:
                    c.update([n])
                    read_count[0] += 1 
                    if log_fails:
                        logfile.write('%s\t%s\n' % (rec.id, n))
                    break # move on to next record and don't do else:
            else: # if runs through all filterfuncs, run this too. 
                read_count[0] += 1 
                yield rec # else yield record when all filters are passed        
    
    # Define single illumina machine filter function if none given 
    if filterfuncs is None:
        def illumina_filter(rec):
            ''' Use machine specific filter         
            N = was not pickup up by machine filter i.e. passed
            Y = was flagged by machine filter i.e. fail 
            '''
            return rec.description.split()[1].split(':')[1] == 'N'
        filterfuncs = [illumina_filter]
        
    # Single Record Cycler      
    RecCycler = Cycler(infiles=infiles, filepattern=filepattern, data_inpath=data_inpath)    

    #===========================================================================
    # MainLoop
    #===========================================================================
      
    # timings
    toc = time.time()
    cum_t = 0 
    for recordgen in RecCycler.seqfilegen:
           
        print '\nProcessing {0}'.format(RecCycler.curfilename)
                               
        # Generator initiated 
        # Record only returned if it passes all filter functions in the list 
        passgen = make_pass_gen(recordgen, filterfuncs)
                          
        # Construct file names
        name = RecCycler.curfilename.split('.')  
        pass_filename = [name[0] + '-pass'] + name[1:] 
        pass_filename = os.path.join(outpath, '.'.join(pass_filename))
        name = '.'.join(name)
    
        if name.endswith('.bgzf'):
            pass_filehdl = bgzf.BgzfWriter(pass_filename)
        elif name.endswith('.fastq'):
            pass_filehdl = open(pass_filename, 'wb')
        elif name.endswith('.gz'):
            pass_filehdl = gzip.open(pass_filename, 'wb')
        else:
            print 'Input file format not supported: %s' % name
            sys.exit()
                   
        print 'Writing passes to \n{0} ....'.format(pass_filename)
        numwritten = SeqIO.write(passgen , pass_filehdl , 'fastq')
        pass_filehdl.close()
        print '{0} records written'.format(numwritten)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print 'Finished file {0} after {1}'.format(RecCycler.curfilenum, 
                                time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
    print '\nFilter stats'
    print '\nFilter No.\tHits'    
    for n in range(len(filterfuncs)):
        print '%s\t\t%s' % (n,c[n])

    print '\nTotal No. Reads Processed:  %s' % read_count[0]
    print '\nTotal No. filtered:  %s (%.2f %%)' % (sum(c.values()), 
                                                    100 * (sum(c.values())/ 
                                                           float(read_count[0])))
    # Clean up
    if log_fails:     
        logfile.flush()
        logfile.close()
    
    total_t = time.time() - toc    
    print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                        time.gmtime(total_t)))
    os.chdir(starting_dir)

def setup_illumina_filter():
    ''' Returns function based on illumina machine filter         
    N = was not pickup up by machine filter i.e. passed
    Y = was flagged by machine filter i.e. fail 
    '''
    # Define filter function
    def f(rec):
        ''' filter function for phred '''
        return rec.description.split()[1].split(':')[1] == 'N'
    return f  
    
def setup_phred_filter(value):
    ''' Returns a filter function based on the mean phred of each read.
    
    If a read has a mean phred of less than the given value, the filter returns
     false.   
    '''
    # Define filter function
    def f(rec):
        ''' filter function for phred '''
        return np.array(rec.letter_annotations['phred_quality']).mean() > value
    return f  
    
def setup_propN_filter(value):
    ''' Returns a filter function based on the proportion of Ns in the read.

    If a read has a higher proportion of Ns than the given value, the filter 
    returns false.
    '''
    # Define filter function
    def f(rec):
        ''' Filter function for propN'''
        return float(rec.seq.count('N')) / len(rec.seq) < value
    return f
            
def setup_cutsite_filter(target_cutsite, mindist=1):
    ''' Returns a filter function based on the match of the read cutsite to the 
    target_cutsite given.
    
    Reads that differ in edit distance by more than mindist, cause filter to 
    return false 
    
    '''
    cutsite_length = len(target_cutsite)
    # Define filterfunc
    def f(rec):
        ''' Filter function for cutsite'''
        
        cutsite = rec.seq[6: 6 + cutsite_length].tostring()
        cutsite_dist = ed.distance(target_cutsite, cutsite)
        
        return cutsite_dist <= mindist
    return f
    
def setup_overhang_filter(target_cutsite, overhang, mindist=0):
    ''' Returns a filter function based on the overhang part of the cutsite. 
    
    The cut site should end with the specified overhang. Those that dont are likely 
    to be genetic contaminants which have been inedvertantly sequenced, and 
    therefore should be discarded. 
       
    Reads that mismatch in the overhang region by more than mindist, cause the 
    filter to return false. 
    '''   
    cutsite_length = len(target_cutsite)
    overhang_length = len(overhang)
    # Define filterfunc
    def f(rec):
        ''' Filter function for cutsite'''
        
        cutsite = rec.seq[6: 6 + cutsite_length].tostring() 
        if cutsite.endswith(overhang):
            return True
        else:
            overhang_dist = ed.distance(cutsite[-overhang_length:], overhang)
            return overhang_dist <= mindist
    return f



        def f(rec):
            ''' Filter function for cutsite'''
            
            # Must calculate MIDlength, but this may vary between files
            if self.c.barcode_files_setup == 'individual':
                
                fname = self.current_file.split('.')[0].split('-')[0]
                
                if f.target_file is None or f.target_file != fname:
                    f.target_file = fname
                    f.MIDlength =  len(self.c.MIDtags[f.target_file].keys()[0])
                
            if self.c.barcode_files_setup == 'global':
                if f.MIDlength is None: 
                    f.MIDlength = len(self.c.MIDtags.keys()[0])
                
            cutsite = rec.seq[f.MIDlength: f.MIDlength + cutsite_length].tostring()
            cutsite_dist = ed.distance(target_cutsite, cutsite)
            
            return cutsite_dist <= max_edit_dist
        
        f.target_file = None
        f.MIDlength = None

def make_MIDdict(infiles=None, filepattern=False, data_inpath=''):
    ''' Function to load in MIDs from a list of files into a dictionary '''

    # Handle multiple types of input for infiles
    assert infiles is not None, 'No files listed or file pattern specified.'         
    if filepattern:
        # Fetch files by file types using glob
        import glob 
        st_dir = os.getcwd()
        os.chdir(data_inpath)
        infiles = glob.glob(infiles)
        os.chdir(st_dir)
    elif type(infiles) == str:
        # Convert to list
        infiles = [infiles]

    # Run through files and store barcodes in a Dictionary object.
    # keys are starting tags (MID (6 BP) + cutsite (6BP))
    tags = {}
    
    for filename in infiles:
        with open(os.path.join(data_inpath,filename), 'rb') as f:
            for line in f:
                elem = line.split()
                tags[elem[0]] = elem[1] 
    return tags



def process_MIDtag(infiles=None, barcodes=None, filepattern=False, 
                   barcode_pattern=False, data_inpath='', barcode_path='',
                   outfile_postfix='-clean', outdir='', 
                   MIDtag_len = 6, max_edit_dist = 1, cutsite_len = 6):
    ''' Goes through files and corrects any errors in MIDtag and cutsite
    
    TODO:
    Currently MIDtag_len and cutsite_len are unused. make these constants in 
    the larger scope of the pipeline to account for different MID tag/ cut site 
    possibilities.
    
    '''

    # Construct Tag dictionary
    MIDdict = make_MIDdict(infiles=barcodes, filepattern=barcode_pattern,
                           data_inpath=barcode_path)
    # Setup Record Cycler
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, data_inpath=data_inpath)
    
    keys = [key[:6] for key in MIDdict.iterkeys()]
    keys.sort()
    
    cutsite = 'TGCAGG' 
    
    # Make ouput directory if required
    outpath = os.path.join(data_inpath, outdir)
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
 
    class Read_Corrector_Class():
    
        def __init__(self):
            ''' Class for a Generator to yield reads that pass a quality check or are corrected 
            
            Keeps track of number of skipped and corrected reads as class attributes.
            '''
            self.skipped_count = 0
            self.corrected_count = 0
    
        def clean_reads_gen(self, recgen, keys):
            
            for rec in recgen:
                
                #===============================================================
                # Analise MID tag section
                #===============================================================
                recMID = str(rec.seq[:6])
                if recMID not in keys:
                    # Sequencing error in the tag. Work out nearest candidate.
                    distvec = np.array([ed.distance(recMID, key) for key in keys]) 
                    
                    min_dist_candidates = []
                    distvec_min = distvec.min()
                    for elem in distvec:
                        if elem == distvec_min and elem <= max_edit_dist:
                            min_dist_candidates.append(keys[elem])

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
                        # Letter annotations must be removed before editing rec.seq
                        temp_var = rec.letter_annotations
                        rec.letter_annotations = {}
                        # Change seq to mutableseq
                        rec.seq = rec.seq.tomutable()
                        rec.seq[:6] = min_dist_candidates[0]
                        rec.seq = rec.seq.toseq()
                        rec.letter_annotations.update(temp_var)
                        self.corrected_count += 1
                        
                #===============================================================
                # Analise Cut site section
                #===============================================================
                rec_cutsite = str(rec.seq[6:12])
                if rec_cutsite != cutsite:
                    # Sequencing error in the cutsite. Correct if less than max_edit_dist
                    
                    cutsite_dist = ed.distance(rec_cutsite, cutsite)

                    if cutsite_dist <= max_edit_dist:
                        # Correct the erroneous cutsite with the actual cutsite.
                        # Letter annotations must be removed before editing rec.seq
                        temp_var = rec.letter_annotations
                        rec.letter_annotations = {}
                        # Change seq to mutableseq
                        rec.seq = rec.seq.tomutable()
                        rec.seq[6:12] = cutsite
                        rec.seq = rec.seq.toseq()
                        rec.letter_annotations.update(temp_var)
                        self.corrected_count += 1
                    else:
                        # Amount of error in cute site is too large to correct 
                        # May also be from contaminants. So read is skipped altogether 
                        self.skipped_count += 1
                        continue
                # Note: rec only yielded if the MID tag and the cutsite are 
                # sucessfully cleaned         
                yield rec
            
    #===========================================================================
    # MAIN LOOP
    #===========================================================================
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
        
        # write records that are cleaned up
        numwritten = SeqIO.write(ReadCorrector.clean_reads_gen(seqfile, keys), 
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
 

def trim_reads(infiles=None, filepattern=False, data_inpath='', 
                out_filename='out_filename.fasta', outpath='', m=12, n=1):
    ''' Trims off the MID tag of each read, as well as the last 'n' bases.
    Writes the trimed reads to one large fasta file for clustering
    
    m = length of MID tag plus cutsite used
    '''
    
    start_dir = os.getcwd() 
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, data_inpath=data_inpath)
    count = 0
    outfile_part_list = []

    print ('Removing MID tags, triming reads and converting {0} files to'
           ' fasta format').format(RecCycler.numfiles)

    # Generator to trim off MID tag and end of read.
    for seqfilegen in RecCycler.seqfilegen:
        
        read_gen = (rec[m:-n] for rec in seqfilegen)
        # File name 
        outfile_part = 'output_part' + str(count) + '.fasta'
        outfile_part_list.append(outfile_part)
        count += 1
        with open(os.path.join(outpath, outfile_part), 'wb') as f:
            write_count = SeqIO.write(read_gen, f, 'fasta')
            print 'Wrote {0} records to file\n{1}'.format(write_count, outfile_part)
    
    # Combine output parts into one big file
    if outpath:
        os.chdir(outpath)
    cmd = ['cat'] + outfile_part_list 
    with open(os.path.join(outpath, out_filename), 'wb') as f:
        print 'Running "{0}" and saving to\n{1}'.format(cmd, os.path.join(outpath, out_filename))
        call(cmd, stdout=f) 
    
    print 'Done, cleaning up temp files ....'
    for f in outfile_part_list:
        os.remove(f)   
    os.chdir(start_dir)

    print '\nPreprocessing Complete.'
    