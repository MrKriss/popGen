'''
Created on 19 Nov 2012

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

import pdb

from utils import smartopen, Cycler, make_MIDdict

class ConfigClass(object):
    pass

class Workflow(object):
    
    def __init__(self, c=None):
        
        # Setup config file
        if c:
            self.c = c
        else:
            self.c = ConfigClass()
            # Defaults for other parameters
            self.filter_functions = None 
            self.c.barcode_files_setup = None   
        
        # Set raw files for next input 
        self.next_input_files = c.raw_input_files
        self.next_input_path = c.inpath
        
        # Input checks
        if c.barcode_files_setup == 'individual':
            # One barcode file per input file with matching names
            barnames = [b.split('.')[0] for b in c.barcode_files]
            filenames = [f.split('.')[0] for f in c.raw_input_files]
            for fname in filenames:
                if fname not in barnames:
                    raise Exception('Set to individual barcode files, yet at least one input'
                    'file name does not match the given barcode file names')
                               
        # Define MID tag Dictionary
        if c.barcode_files_setup == 'individual':
            # Process each file with its own list of barcodes
            # tags_per_filename = {'filename' :  {'MIDtag' : 'individual' }}
            tags_per_filename = {}
            
            for input_filename in c.raw_input_files:
                name = input_filename.split('.')[0]
                tags_per_filename[name] = {}
                with open(os.path.join(c.barpath,name) + '.txt', 'rb') as f:
                    for line in f:
                        elem = line.split()
                        tags_per_filename[name][elem[0]] = elem[1]
            c.MIDtags = tags_per_filename 
        elif c.barcode_files_setup == 'global':
            # Process each file with reference to a global list of barcodes
            # global_tags = {'MIDtag' : 'individual' }
            global_tags = {}
            
            for bar_filename in c.barcode_files:
                name = bar_filename.split('.')[0]
                with open(os.path.join(c.barpath,name) + '.txt', 'rb') as f:
                    for line in f:
                        elem = line.split()
                        global_tags[elem[0]] = elem[1]
            c.MIDtags = global_tags 
        else:
            raise Exception('Barcode file usage not specified.'
            'Set c.barcode_files_setup to "individual" or "global"')
           
                
    def set_input_files(self, infiles, file_pattern, inpath):
        
        starting_dir = os.getcwd()
        
        if file_pattern:
            if inpath:
                os.chdir(inpath)
            else:
                os.chdir(self.c.paths.inpath)
            raw_files = glob.glob('infiles')
            raw_files.sort()
            self.c.next_input_files = raw_files
        else:
            if not type(infiles) == list:
                infiles = [infiles]
            self.c.next_input_files = infiles
        
        if os.getcwd() != starting_dir:
            os.chdir(starting_dir)
    
    def filter_reads_pipeline(self):
    
        ''' Filter reads based on criteria specified by a sequence of functions given
         in self.filter_functions.
        
        Indirvidual functions in self.filter_functions must take in a sequence record object, 
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
        c = self.c
         
        # Input checks       
        if self.filter_functions is None :
            raise Exception('No filter functions defined')
        
        #===============================================================================
        # Setup variables and utility functions 
        #===============================================================================
           
        # Path variables
        starting_dir = os.getcwd()
        
        if not os.path.exists(c.filteredpath):
            os.makedirs(c.filteredpath) 
            
        outpath = c.filteredpath
        
        # Setup Counters
        count = Counter()
        read_count = [0]
            
        # Define filter function loop
        if c.log_fails: 
            # Set up log file
            print '\nSaving Failed reads to log file: fails.log'
            if os.path.isfile(os.path.join(outpath, 'fails.log')):
                print 'Overwriting previous Log file'
                logfile = open(os.path.join(outpath, 'fails.log'), 'wb')
            else:
                logfile = open(os.path.join(outpath, 'fails.log'), 'wb')
                   
        def make_pass_gen(recordgen):
            ''' An easier to read implimentation of a generator to only yield 
            records that pass all filter functions specified'''
            
            filterfuncs = self.filter_functions
            
            for rec in recordgen:
                for n, func in enumerate(filterfuncs):
                    # Run filter functions on record       
                    if func(rec) == False:
                        count.update([n])
                        read_count[0] += 1 
                        if c.log_fails:
                            logfile.write('%s\t%s\n' % (rec.id, n))
                        break # move on to next record and don't do else:
                else: # if runs through all filterfuncs, run this too. 
                    read_count[0] += 1 
                    yield rec # else yield record when all filters are passed        
                    
        # Single Record Cycler      
        RecCycler = Cycler(infiles=self.next_input_files, filepattern=False, 
                           inpath=self.next_input_path)    
    
        #===========================================================================
        # MainLoop
        #===========================================================================
          
        # timings
        toc = time.time()
        cum_t = 0
        outnames = []
         
        for recordgen in RecCycler.seqfilegen:
            
            self.current_file = RecCycler.curfilename
            print '\nProcessing {0}'.format(self.current_file)
               
            # Generator initiated 
            # Record only returned if it passes all filter functions in the list 
            passgen = make_pass_gen(recordgen)
                              
            # Construct file names
            name = RecCycler.curfilename.split('.')  
            pass_filename = '.'.join([name[0] + '-pass'] + name[1:]) 
            pass_file = os.path.join(outpath, pass_filename)
            name = '.'.join(name)
            outnames.append(pass_filename) 
            
            if name.endswith('.bgzf'):
                pass_filehdl = bgzf.BgzfWriter(pass_file)
            elif name.endswith('.fastq'):
                pass_filehdl = open(pass_file, 'wb')
            elif name.endswith('.gz'):
                pass_filehdl = gzip.open(pass_file, 'wb')
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
        for n in range(len(self.filter_functions)):
            print '%s\t\t%s' % (n,count[n])
    
        print '\nTotal No. Reads Processed:  %s' % read_count[0]
        print '\nTotal No. filtered:  %s (%.2f %%)' % (sum(count.values()), 
                                                        100 * (sum(count.values())/ 
                                                               float(read_count[0])))
        # Clean up
        if c.log_fails:     
            logfile.flush()
            logfile.close()
        
        total_t = time.time() - toc    
        print 'Processed all files in {0}'.format(time.strftime('%H:%M:%S', 
                                                            time.gmtime(total_t)))
        os.chdir(starting_dir)

        # Update internal Variables
        self.next_input_files = outnames
        self.next_input_path = outpath

    def process_MIDtag(self, max_edit_dist = None, outfile_postfix='-clean'):
        ''' Goes through files and corrects any errors in MIDtag and cutsite
        '''
        c = self.c
        
        if max_edit_dist is None:
            max_edit_dist = c.max_edit_dist
        

        # Setup Record Cycler
        RecCycler = Cycler(infiles=self.next_input_files, 
                           filepattern=False, inpath=self.next_input_path)
                
        # Make ouput directory if required
        if not os.path.exists(c.processedpath):
            os.makedirs(c.processedpath) 
            
        outpath = c.processedpath
      
        class Read_Corrector_Class():
        
            def __init__(self, tags,cutsite):
                ''' Class for a Generator to yield reads that pass a quality check or are corrected 
                
                Keeps track of number of skipped and corrected reads as class attributes.
                '''
                self.skipped_count = 0
                self.corrected_count = 0
                
                self.barMIDs = tags.keys()
                self.barMIDs.sort()
                
                self.cutsite = cutsite
                
                # Check length of MIDs
                lengths = [len(mid) for mid in self.barMIDs]
                assert all([x == lengths[0] for x in lengths]), 'Error: Lengths of MIDs are not equal'
                self.MIDlength = lengths[0]
        
            def clean_reads_gen(self, recgen):
                
                for rec in recgen:
                    
                    #===============================================================
                    # Analise MID tag section
                    #===============================================================
                    recMID = str(rec.seq[:self.MIDlength])
                    
                    if recMID not in self.barMIDs:
                        # Sequencing error in the tag. Work out nearest candidate.
                        distvec = np.array([ed.distance(recMID, tag) for tag in self.barMIDs]) 
                        
                        min_dist_candidates = []
                        distvec_min = distvec.min()
                        for elem in distvec:
                            if elem == distvec_min and elem <= max_edit_dist:
                                min_dist_candidates.append(self.barMIDs[elem])
    
                        if len(min_dist_candidates) > 1:
                            # Muliple candidates. True MID is Ambiguous 
                            self.skipped_count += 1
                            continue
                        elif len(min_dist_candidates) == 1:
                            # Correct the erroneous tag with the candidate.
                            # Letter annotations must be removed before editing rec.seq
                            temp_var = rec.letter_annotations
                            rec.letter_annotations = {}
                            # Change seq to mutableseq
                            rec.seq = rec.seq.tomutable()
                            rec.seq[:self.MIDlength] = min_dist_candidates[0]
                            rec.seq = rec.seq.toseq()
                            rec.letter_annotations.update(temp_var)
                            self.corrected_count += 1
                            
                    #===============================================================
                    # Analise Cut site section
                    #===============================================================
                    rec_cutsite = str(rec.seq[self.MIDlength: self.MIDlength + len(self.cutsite)])
                    if rec_cutsite != self.cutsite:
                        # Sequencing error in the cutsite. Correct if less than max_edit_dist
                        
                        cutsite_dist = ed.distance(rec_cutsite, self.cutsite)
    
                        if cutsite_dist <= max_edit_dist:
                            # Correct the erroneous cutsite with the actual cutsite.
                            # Letter annotations must be removed before editing rec.seq
                            temp_var = rec.letter_annotations
                            rec.letter_annotations = {}
                            # Change seq to mutableseq
                            rec.seq = rec.seq.tomutable()
                            rec.seq[self.MIDlength: self.MIDlength + len(self.cutsite)] = self.cutsite
                            rec.seq = rec.seq.toseq()
                            rec.letter_annotations.update(temp_var)
                            self.corrected_count += 1
                        else:
                            # Amount of error in cute site is too large to correct 
                            # May also be from contaminants. So read is skipped altogether 
                            self.skipped_count += 1
                            continue
                    # Note: rec only yielded if the MID tag and the cutsite are 
                    # sucessfuly cleaned         
                    yield rec
                
        #===========================================================================
        # MAIN LOOP
        #===========================================================================
        total_numwritten = 0
        total_numskipped = 0
        total_numcorrected = 0
        toc = time.time()
        cum_t = 0
        outnames = []
        
        if c.barcode_files_setup == 'global':
            tags = c.MIDtags 
        
        for seqfile in RecCycler.seqfilegen:
                                    
            if c.barcode_files_setup == 'individual':
                current_filename = RecCycler.curfilename.split('.')[0]                
                tags = c.MIDtags[current_filename.split('-')[0]]
                  
            # Make reads class for generator
            ReadCorrector = Read_Corrector_Class(tags,c.cutsite)

            # Set output filename and handle
            filename = RecCycler.curfilename
            filename = filename.split('.')
            filename[0] = filename[0] + outfile_postfix
            filename = '.'.join(filename)
            outnames.append(filename)
            
            outfile_path = os.path.join(outpath, filename)        
            output_filehdl = bgzf.BgzfWriter(outfile_path, mode='wb')
            # write records that are cleaned up
            numwritten = SeqIO.write(ReadCorrector.clean_reads_gen(seqfile), 
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
        # Update internal Variables
        self.next_input_files = outnames
        self.next_input_path = outpath

    def trim_reads(self, outfile='outfile.fasta', outpath = '', n = 1):
        ''' Trims off the MID tag of each read, as well as the last 'n' bases.
        Writes the trimed reads to one large fasta file for clustering'''
        
        c = self.c
        
        start_dir = os.getcwd() 
        
        RecCycler = Cycler(infiles=self.next_input_files, 
                           filepattern=False, inpath=self.next_input_path)
        count = 0
        outfile_part_list = []
    
        print ('Removing MID tags, triming reads and converting {0} files to'
               ' fasta format').format(RecCycler.numfiles)
    
        if c.barcode_files_setup == 'global':
            read_start_idx = len(c.cutsite) + len(c.MIDtags.keys()[0])
             
        # Generator to trim off MID tag and end of read.
        for seqfilegen in RecCycler.seqfilegen:
                       
            if c.barcode_files_setup == 'individual':
                fname = RecCycler.curfilename.split('.')[0].split('-')[0]
                lenMIDs = len(c.MIDtags[fname].keys()[0])
                read_start_idx = len(c.cutsite) + lenMIDs  
    
            read_gen = (rec[read_start_idx:-n] for rec in seqfilegen)
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
        with open(os.path.join(outpath, outfile), 'wb') as f:
            print 'Running "{0}" and saving to\n{1}'.format(cmd, os.path.join(outpath, outfile))
            call(cmd, stdout=f) 
        
        print 'Done, cleaning up temp files ....'
        for f in outfile_part_list:
            os.remove(f)   
        os.chdir(start_dir)

    def make_illumina_filter(self):
        ''' Returns filtering function based on illumina machine filter         
        N = was not pickup up by machine filter i.e. passed
        Y = was flagged by machine filter i.e. fail 
        '''
        # Define filter function
        def f(rec):
            ''' filter function for phred '''
            return rec.description.split()[1].split(':')[1] == 'N'
        return f  
        
    def make_phred_filter(self, value):
        ''' Returns a filtering function based on the mean phred of each read.
        
        If a read has a mean phred of less than the given value, the filter returns
         false.   
        '''
        # Define filter function
        def f(rec):
            ''' filter function for phred '''
            return np.array(rec.letter_annotations['phred_quality']).mean() > value
        return f  
        
    def make_propN_filter(self, value):
        ''' Returns a filter function based on the proportion of Ns in the read.
    
        If a read has a higher proportion of Ns than the given value, the filter 
        returns false.
        '''
        # Define filter function
        def f(rec):
            ''' Filter function for propN'''
            return float(rec.seq.count('N')) / len(rec.seq) < value
        return f
                
    def make_cutsite_filter(self, target_cutsite=None, max_edit_dist=None):
        ''' Returns a filter function based on the match of the read cutsite to the 
        target_cutsite given.
        
        Reads that differ in edit distance by more than mindist, cause filter to 
        return false 
        
        '''

        if target_cutsite is None:
            target_cutsite = self.c.cutsite

        if max_edit_dist is None:
            max_edit_dist = self.c.max_edit_dist

        cutsite_length = len(target_cutsite)
    
        # Define filterfunc
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
        
        return f
        
    def make_overhang_filter(self, target_cutsite=None, overhang=None, max_edit_dist=0):
        ''' Returns a filter function based on the overhang part of the cutsite. 
        
        The cut site should end with the specified overhang. Those that dont are likely 
        to be genetic contaminants which have been inedvertantly sequenced, and 
        therefore should be discarded. 
           
        Reads that mismatch in the overhang region by more than mindist, cause the 
        filter to return false. 
        '''   
    
        if target_cutsite is None:
            target_cutsite = self.c.cutsite
            overhang = target_cutsite[-2:]
        
        cutsite_length = len(target_cutsite)
        overhang_length = len(overhang)
        # Define filterfunc
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
            if cutsite.endswith(overhang):
                return True
            else:
                overhang_dist = ed.distance(cutsite[-overhang_length:], overhang)
                return overhang_dist <= max_edit_dist
            
        f.MIDlength = None
        f.target_file = None
            
        return f
        
#===============================================================================
# Individual Functions 
#===============================================================================

def filter_reads(infiles=None, filepattern='', inpath='', filterfunc=None, 
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

def filter_reads_pipeline(infiles=None, filepattern='', inpath='', filterfuncs=None, 
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
    if inpath:
        outpath = os.path.join(inpath, outdir) 
    else:        
        outpath = os.path.join(starting_dir, outdir) 
    
    if not os.path.isdir(outpath):
        if inpath:
            os.chdir(inpath)
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
    RecCycler = Cycler(infiles=infiles, filepattern=filepattern, inpath=inpath)    

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
    
def process_MIDtag(infiles=None, barcodes=None, filepattern=False, 
                   barcode_pattern=False, inpath='', barcode_path='',
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
                           inpath=barcode_path)
    # Setup Record Cycler
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, inpath=inpath)
    
    keys = [key[:6] for key in MIDdict.iterkeys()]
    keys.sort()
    
    cutsite = 'TGCAGG' 
    
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
 
def file2bgzf(infiles=None, filepattern=False, inpath='', SQLindex=True):
    ''' Convert given list of files from .gz or .fastq to .bgzf,
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
        
        f = smartopen(filename)
        
        # Checks for type of input
        if filename.endswith('.gz'):
            # Drop .gz and append .bgzf
            bgzfFileName = '.'.join(filename.split('.')[:-1]) + '.bgzf'  
        elif filename.endswith('.fastq'):
            # Append .bgzf
            bgzfFileName = filename + '.bgzf'

        print "Producing BGZF output from {0}...".format(filename)
        w = bgzf.BgzfWriter(bgzfFileName, 'wb')
        while True:
            data = f.read(65536)
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

def trim_reads(infiles=None, filepattern=False, inpath='', 
                outfile='outfile.fasta', outpath = '', n = 1):
    ''' Trims off the MID tag of each read, as well as the last 'n' bases.
    Writes the trimed reads to one large fasta file for clustering'''
    
    start_dir = os.getcwd() 
    
    RecCycler = Cycler(infiles=infiles, 
                       filepattern=filepattern, inpath=inpath)
    count = 0
    outfile_part_list = []

    print ('Removing MID tags, triming reads and converting {0} files to'
           ' fasta format').format(RecCycler.numfiles)

    # Generator to trim off MID tag and end of read.
    for seqfilegen in RecCycler.seqfilegen:
        
        read_gen = (rec[12:-n] for rec in seqfilegen)
        # File name 
        outfile_part = 'output_part' + str(count) + '.fasta'
        outfile_part_list.append(outfile_part)
        count += 1
        with open(os.path.join(outpath, outfile_part), 'wb') as f:
            write_count = SeqIO.write(read_gen, f, 'fasta')
            print 'Wrote {0} records to file\n{1}'.format(write_count, outfile_part)
    
    # Combine output parts into one big file
    os.chdir(outpath)
    cmd = ['cat'] + outfile_part_list 
    with open(os.path.join(outpath, outfile), 'wb') as f:
        print 'Running "{0}" and saving to\n{1}'.format(cmd, os.path.join(outpath, outfile))
        call(cmd, stdout=f) 
    
    print 'Done, cleaning up temp files ....'
    call('rm output*', shell=True)    
    os.chdir(start_dir)

if __name__ == '__main__':
    
#==============================================================================
    ''' RUNS SCRIPT FOR SMALL TEST SET '''
#===============================================================================

    starting_dir = os.getcwd()
    
    LANE = '6'
    
    # Path setup
    inpath = '/home/pgrad/musselle/ubuntu/workspace/popGen/testdata'
    barpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    raw_files = ['small_test_set.fastq']
    
    # Setup Filters
    outdir = 'filtered_reads'
    filter_functions = [setup_illumina_filter(), 
                        setup_propN_filter(0.1),
                        setup_phred_filter(25),
                        setup_cutsite_filter('TCGAGG', 2),
                        setup_overhang_filter('TCGAGG', 'GG', 0)]
    
    filter_reads_pipeline(infiles=raw_files, inpath=inpath, filterfuncs=filter_functions, 
                              outdir=outdir, log_fails=True)
    # Update names and path
    filtered_files = []
    for name in raw_files:
        temp = name.split('.')
        temp[0] = temp[0] + '-pass'
        temp = '.'.join(temp) 
        filtered_files.append(temp)
    filtered_inpath = os.path.join(inpath, outdir)
      
    cleaned_file_postfix = '-clean' 
    cleaned_outdir = '' # 'cleaned_data'
    barcode_pattern = '*[' + LANE + '].txt'

    process_MIDtag(infiles=filtered_files, barcodes =barcode_pattern,
                   barcode_pattern=True, 
                   inpath=filtered_inpath, barcode_path=barpath,
                   outfile_postfix=cleaned_file_postfix, outdir=cleaned_outdir, 
                   MIDtag_len = 6, max_edit_dist = 1, cutsite_len = 6)
