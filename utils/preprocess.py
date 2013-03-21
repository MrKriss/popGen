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

#from general_utilities import set_trace

from utils import smartopen, Cycler

class ConfigClass(object):
    pass

class Preprocessor(object):
    ''' Holds all functions in preprocessing pipeline. 
    
    Functions are useable individually, in which case the return value is in 
    the format:
                    out = (list of output files, output path)
    
    This output can then be fed to the next function of interest as: 
                out2 = func2(infiles=out[0], inpath=out[1])
    If infiles and inpath are omitted, the output from the last function is used by default.
    
    Other option is to use the workflow function to chain the outputs of function together
    
    
    '''
    
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
        self.next_input_path = c.data_inpath
        
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
                with open(os.path.join(c.barcode_inpath,name) + '.txt', 'rb') as f:
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
                with open(os.path.join(c.barcode_inpath,name) + '.txt', 'rb') as f:
                    for line in f:
                        elem = line.split()
                        global_tags[elem[0]] = elem[1]
            c.MIDtags = global_tags 
        else:
            raise Exception('Barcode file usage not specified.'
            'Set c.barcode_files_setup to "individual" or "global"')
                           
    def set_input_files(self, infiles, file_pattern, data_inpath):
        
        starting_dir = os.getcwd()
        
        if file_pattern:
            if data_inpath:
                os.chdir(data_inpath)
            else:
                os.chdir(self.c.paths.data_inpath)
            raw_files = glob.glob('infiles')
            raw_files.sort()
            self.c.next_input_files = raw_files
        else:
            if not type(infiles) == list:
                infiles = [infiles]
            self.c.next_input_files = infiles
        
        if os.getcwd() != starting_dir:
            os.chdir(starting_dir)
    
    def filter_reads_pipeline(self, infiles=None, inpath=None):
    
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
        if not os.path.exists(c.filtered_outpath):
            os.makedirs(c.filtered_outpath) 
            
        outpath = c.filtered_outpath
        
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
        
        # Setup Record Cycler 
        if infiles is None:
            infiles = self.next_input_files
        if inpath is None:
            inpath = self.next_input_path
 
        RecCycler = Cycler(infiles=infiles, filepattern=False, data_inpath=inpath)    
    
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
        
        return (outnames, outpath)

    def process_MIDtag(self, infiles=None, inpath=None, max_edit_dist=None, 
                       outfile_postfix='-clean'):
        ''' Goes through files and corrects any errors in MIDtag and cutsite
        '''
        c = self.c
        
        if max_edit_dist is None:
            max_edit_dist = c.max_edit_dist
        
        # Setup Record Cycler
        if infiles is None:
            infiles = self.next_input_files
        if inpath is None:
            inpath = self.next_input_path
        
        RecCycler = Cycler(infiles=infiles, filepattern=False, data_inpath=inpath)
                
        # Make ouput directory if required
        if not os.path.exists(c.tag_processed_outpath):
            os.makedirs(c.tag_processed_outpath) 
            
        outpath = c.tag_processed_outpath
      
        class Read_Corrector_Class():
        
            def __init__(self, tags, cutsite):
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
        
        return (outnames, outpath)

    def trim_reads(self, infiles=None, inpath=None, out_filename=None, outpath = None, n = 1):
        ''' Trims off the MID tag of each read, as well as the last 'n' bases.
        Writes the trimed reads to one large fasta file for clustering'''
        
        c = self.c
        
        if outpath is None:
            outpath = c.tag_processed_outpath
        if out_filename is None:
            out_filename = c.experiment_name + '_all_preprocessed.fasta'
        
        start_dir = os.getcwd() 

        # Setup Record Cycler        
        if infiles is None:
            infiles = self.next_input_files
        if inpath is None:
            inpath = self.next_input_path
        
        RecCycler = Cycler(infiles=infiles, filepattern=False, data_inpath=inpath)
        
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
        
        with open(os.path.join(outpath, out_filename), 'wb') as f:
            print 'Running "{0}" and saving to\n{1}'.format(cmd, os.path.join(outpath, out_filename))
            call(cmd, stdout=f) 
        
        print 'Done, cleaning up temp files ....'
        for f in outfile_part_list:
            os.remove(f)   
        os.chdir(start_dir)
        
        print '\nPreprocessing Complete.'
        
        self.next_input_files = out_filename
        self.next_input_path = outpath
        
        return (out_filename, outpath)

#    def split_by_tags(self, infiles=None, inpath=None, outpath=None, out_filename=None,
#                      report=True, savecounter=True):
#        ''' Split the file into separate files based on MID tags '''
#        
#        c = self.c
#        
#        if outpath is None:
#            outpath = c.tag_processed_outpath
#        i
#            print 'Database found with matching file name.'f out_filename is None:
#            out_filename = c.experiment_name
#             
#        start_dir = os.getcwd() 
#
#        # Setup Record Cycler        
#        if infiles is None:
#            infiles = self.next_input_files
#        if inpath is None:
#            inpath = self.next_input_path
#        
#        RecCycler = Cycler(infiles=infiles, filepattern=False, data_inpath=inpath)
#        
#        tag_counter = Counter
#        outfile_files_list = []
#    
#        if c.barcode_files_setup == 'individual':
#            # Process each file with its own list of barcodes
#            # tags_per_filename = {'filename' :  {'MIDtag' : 'individual' }}
#            
#        
#            for seqfilegen in RecCycler.seqfilegen:
#                
#                fname = RecCycler.curfilename.split('.')[0].split('-')[0]
#                lenMIDs = len(c.MIDtags[fname].keys()[0])
#                read_start_idx = len(c.cutsite) + lenMIDs
#                
#                
#        
#        
#        
#        
#        
#        
#        
#        elif c.barcode_files_setup == 'global':
#            # Process each file with reference to a global list of barcodes
#            # global_tags = {'MIDtag' : 'individual' }
#            
#            MID_length = len(c.MIDtags.keys()[0])
#            
#            print ('Spliting {0} files into a total of {1} files bases on MID tags'
#               '').format(RecCycler.numfiles, len(c.MIDtags))
#            
#            # Open Files for Writing for each tag  
#            for tag, individual in c.MIDtags.iteritems():
#                
#                fname = '-'.join(out_filename, tag, individual)
#                outfile_files_list.append(fname)
#                fvarname = 'f-' + tag
#                vars()[fvarname] = open(os.path.join(outpath, fname), 'w')
#    
#            for rec in RecCycler.recgen:
#        
#                tag = rec.seq[:MID_length].tostring()
#                
#                if tag not in c.MIDtags:
#                    raise Exception('MID tag not found in barcode library')
#                else:                   
#                    fvarname = 'f-' + tag                 
#                    SeqIO.write(rec, vars()[fvarname], 'fastq');
#                    tag_counter[tag] += 1
#
#            # Flush and Close Files for each tag  
#            for tag, individual in c.MIDtags.iterkeys():
#
#                fvarname = 'f-' + tag
#                vars()[fvarname].flush()
#                vars()[fvarname].close()
#
#            
#            print 'Wrote {0} records to file\n{1}'.format(write_count, outfile_part)
#                
#                    
#                    
#                    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#    
#                read_gen = (rec[read_start_idx:-n] for rec in seqfilegen)
#    
#    
#    
#    
#    
#    
#
#    
#    
#        if c.barcode_files_setup == 'global':
#            
#            read_start_idx = len(c.cutsite) + len(c.MIDtags.keys()[0])
#             
#        # Generator to trim off MID tag and end of read.
#        for seqfilegen in RecCycler.seqfilegen:
#                       
#            if c.barcode_files_setup == 'individual':
#                fname = RecCycler.curfilename.split('.')[0].split('-')[0]
#                lenMIDs = len(c.MIDtags[fname].keys()[0])
#                read_start_idx = len(c.cutsite) + lenMIDs  
#    
#
#        
#        
#        
#        
#        
#        
#        
#        
#        
#        
#        
#        
#       
#            




    def cleanup_files(self, *args):
        ''' Remove intermediate files that are not needed '''      
        # Choices of 'filtered', 'tag_processed', 'all'
        # Raw inputs are always unchanged, and output .fasta is always kept
             
        files2remove = []
          
        if 'filtered' in args or 'all' in args:
            pattern = os.path.join(self.c.filtered_outpath, 
                                   '*' + self.c.filtered_files_postfix) 
            files2remove.extend(glob.glob(pattern))
        elif 'tag_processed' in args or 'all' in args:
            pattern = os.path.join(self.c.tag_processed_outpath, 
                                   '*' + self.c.tag_processed_files_postfix) 
            files2remove.extend(glob.glob(pattern))
        
        # Remove files
        for f in files2remove:
            try:
                os.remove(f)   
            except OSError:
                pass   

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


if __name__ == '__main__':
    
#==============================================================================
    ''' RUNS SCRIPT FOR SMALL TEST SET '''
#===============================================================================

    starting_dir = os.getcwd()
    
    LANE = '6'
    
    # Path setup
    data_inpath = '/home/pgrad/musselle/ubuntu/workspace/popGen/testdata'
    barcode_inpath = '/space/musselle/datasets/gazellesAndZebras/barcodes'
    raw_files = ['small_test_set.fastq']
    
    # Setup Filters
    outdir = 'filtered_reads'
    filter_functions = [setup_illumina_filter(), 
                        setup_propN_filter(0.1),
                        setup_phred_filter(25),
                        setup_cutsite_filter('TCGAGG', 2),
                        setup_overhang_filter('TCGAGG', 'GG', 0)]
    
    filter_reads_pipeline(infiles=raw_files, data_inpath=data_inpath, filterfuncs=filter_functions, 
                              outdir=outdir, log_fails=True)
    # Update names and path
    filtered_files = []
    for name in raw_files:
        temp = name.split('.')
        temp[0] = temp[0] + '-pass'
        temp = '.'.join(temp) 
        filtered_files.append(temp)
    filtered_data_inpath = os.path.join(data_inpath, outdir)
      
    cleaned_file_postfix = '-clean' 
    cleaned_outdir = '' # 'cleaned_data'
    barcode_pattern = '*[' + LANE + '].txt'

    process_MIDtag(infiles=filtered_files, barcodes =barcode_pattern,
                   barcode_pattern=True, 
                   data_inpath=filtered_data_inpath, barcode_path=barcode_inpath,
                   outfile_postfix=cleaned_file_postfix, outdir=cleaned_outdir, 
                   MIDtag_len = 6, max_edit_dist = 1, cutsite_len = 6)
