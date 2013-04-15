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
import re

import numpy as np
from Bio import SeqIO, bgzf
import editdist as ed

from general_utilities import set_trace

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
    def __init__(self, c=None, db=None):
        
        # Setup config file
        if c:
            self.c = c
            # Set raw files for next input 
            self.next_input_files = c.raw_input_files
            self.next_input_path = c.data_inpath     
        else:
            self.c = ConfigClass()
            # Defaults for other parameters
            self.filter_functions = None 
        
        if db:
            self.db = db
        
    def get_data4file(self, filename, fields=['MIDtag']):
        ''' Return a list of records for the given filename with specified fields. 
        Each record is a tuple. 
        
        filename is a globed for without file extension.
        '''

        # Remove filename postfixes if present         
        filename_parts = filename.split('.')
        if filename_parts[0].endswith(self.c.tag_processed_files_postfix):
            filename_parts[0] = filename_parts[0][:-len(self.c.tag_processed_files_postfix)]
        if filename_parts[0].endswith(self.c.filtered_files_postfix):
            filename_parts[0] = filename_parts[0][:-len(self.c.filtered_files_postfix)]
        
        filename = '.'.join(filename_parts)
         
        # Extract tags that match filename of the datafile
        rows = self.db.select('''filename FROM datafiles WHERE filename GLOB ? ''', (filename,))
        assert rows and len(rows) > 0, 'No data files returned from filename querry'  
        rows = [tuple(elem) for elem in rows]
        
        # Check there are no repeats of filename glob
        names = zip(*rows)[0]
        assert len(names) == len(set(names)), 'Filenames in database are not unique'
        
        # Extract MIDtags that match datafiles for the globed filename
        rows = self.db.get_samples4datafile(filename, fields=fields)
        assert rows and len(rows) > 0, 'No MIDtags returned from filename querry'
        rows = [tuple(elem) for elem in rows]
        
        # Check there are no repeats of MIDtags
        tags = zip(*rows)[0]
        assert len(tags) == len(set(tags)), 'MIDtags are not unique for samples in files: {0}'.format(str(names))
        
        return rows

    def set_input_files(self, infiles=None, data_inpath=None, file_pattern=None):
        
        starting_dir = os.getcwd()
        
        if infiles is None:
            # Glob the file pattern 
            if file_pattern:
                if data_inpath:
                    os.chdir(data_inpath)
                files = glob.glob(file_pattern)
                files.sort()
                self.c.next_input_files = files
                self.c.next_input_path = os.getcwd()
            else:
                raise Exception('No input files passed and no file pattern defined')
        else:
            if type(infiles) is str:
                infiles = [infiles]
            if data_inpath:
                self.c.next_input_path = data_inpath
                   
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
        # Write the summary to a file 
        with open(os.path.join(outpath, "filter_summary_" + self.c.experiment_name +  ".log"), 'wb') as f:
            f.write("Filter parameters:\n") 
            f.write("------------------\n")
            f.write(self.c.filter_funtion_params + "\n")
            f.write("\n")
            f.write('Filter stats:\n')
            f.write("-------------------------------------------\n")
            f.write('Filter\t\t\tHits\n')    
            f.write("-------------------------------------------\n")
            for i, x in enumerate(self.filter_functions):
                percent = count[i] / float(sum(count.values())) 
                f.write('%s\t\t%s\t(%.2f %%)\n' % (x.__name__, count[i], percent * 100))
            f.write('\nTotal No. Reads Processed:  \t%s\n' % read_count[0])
            f.write('Total No. filtered:  \t\t%s (%.2f %%)\n' % (sum(count.values()), 
                                                        100 * (sum(count.values())/ 
                                                               float(read_count[0]))))
            f.write('Total No. Passed: \t\t%s (%.2f %%)\n' % (read_count[0] - sum(count.values()), 
                                                        100 * ((read_count[0] - sum(count.values()))/ 
                                                               float(read_count[0]))))
            
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
                self.MIDtag_corrected_count = 0
                self.cutsite_corrected_count = 0
                self.read_corrected_count = 0
                self.total_read_count = 0
                
                self.barMIDs = sorted(tags)
                self.cutsite = cutsite
                
                # Check length of MIDs
                lengths = [len(mid) for mid in self.barMIDs]
                assert all([x == lengths[0] for x in lengths]), 'Error: Lengths of MIDs are not equal'
                self.MIDlength = lengths[0]
        
            def clean_reads_gen(self, recgen):
                
                for rec in recgen:

                    self.total_read_count += 1
                    MID_corrected = False
                    
                    #===============================================================
                    # Analise MID tag section
                    #===============================================================
                    recMID = str(rec.seq[:self.MIDlength])

                    if recMID not in self.barMIDs:
                        # Sequencing error in the tag. Work out nearest candidate.
                        distvec = np.array([ed.distance(recMID, tag) for tag in self.barMIDs]) 
                        
                        distvec_min = distvec.min()
                        
                        if distvec_min > max_edit_dist:
                            self.skipped_count += 1
                            continue
                            
                        min_dist_idxs = np.argwhere(distvec==distvec.min())
                        if min_dist_idxs.size > 1:
                            # Muliple candidates. True MID is Ambiguous 
                            self.skipped_count += 1
                            continue
                    
                        elif min_dist_idxs.size == 1:
                            # Correct the erroneous tag with the candidate.
                            # Letter annotations must be removed before editing rec.seq
                            temp_var = rec.letter_annotations
                            rec.letter_annotations = {}
                            # Change seq to mutableseq
                            rec.seq = rec.seq.tomutable()
                            rec.seq[:self.MIDlength] = self.barMIDs[min_dist_idxs[0]]
                            rec.seq = rec.seq.toseq()
                            rec.letter_annotations.update(temp_var)
                            self.MIDtag_corrected_count += 1
                            self.read_corrected_count += 1
                            MID_corrected = True
                            
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
                            self.cutsite_corrected_count += 1
                            if not MID_corrected:
                                self.read_corrected_count += 1
                        else:
                            # Amount of error in cute site is too large to correct 
                            # May also be from contaminants. So read is skipped altogether 
                            self.skipped_count += 1
                            if MID_corrected:
                                # Roll back counters
                                self.MIDtag_corrected_count -= 1
                                self.read_corrected_count -= 1
                            continue
                    # Note: rec only yielded if the MID tag and the cutsite are 
                    # sucessfuly cleaned         
                    yield rec
                
        #===========================================================================
        # MAIN LOOP
        #===========================================================================
        total_reads = 0 
        total_written = 0 
        total_skipped =  0 
        total_corrected = 0 
        total_MIDcorrected = 0 
        total_cutsite_corrected = 0 
        
        toc = time.time()
        cum_t = 0
        outnames = []
        
        for seqfile in RecCycler.seqfilegen:
            
            tags = self.get_data4file(RecCycler.curfilename, fields=['MIDtag'])
            # tags is returned as a list of tuples for each record            
            tags = zip(*tags)[0]
            # tags is now a tuple of all the first elements in each record  
            
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
    
            # Increment Counters
            total_reads += ReadCorrector.total_read_count
            total_written += numwritten
            total_corrected += ReadCorrector.read_corrected_count
            total_MIDcorrected += ReadCorrector.MIDtag_corrected_count
            total_cutsite_corrected += ReadCorrector.cutsite_corrected_count
            total_skipped += ReadCorrector.skipped_count
            
            print ('{0} records written, of which ' 
            '{1} were corrected').format(numwritten, ReadCorrector.read_corrected_count)
            print '{0} records skipped'.format(ReadCorrector.skipped_count)
            loop_t = time.time() - toc - cum_t
            cum_t += loop_t
            print 'Finished {0} after {1}'.format(filename, 
                            time.strftime('%H:%M:%S', time.gmtime(loop_t)))
    
        print 'Total records written: {0}'.format(total_written)
        print 'Total records skipped: {0}'.format(total_skipped)
        print 'Total of {0} tags corrected.'.format(total_corrected)
        
        # Write the summary to a file 
        with open(os.path.join(outpath, "cleaner_summary_" + self.c.experiment_name + ".log"), 'wb') as f:
            f.write("Cleaner parameters:\n") 
            f.write("------------------\n")
            f.write("Maximum Edit Distance = " + str(max_edit_dist) + "\n")
            f.write("\n")
            f.write('Cleaner Stats:\n')
            f.write("Total reads: {0} \n".format(total_reads))
            f.write("-----------------------------------------------------\n")
            f.write("Count          Percentage      \n")
            f.write("-----------------------------------------------------\n")
            f.write("{0} \t({1:.2f}%)\tWritten to files\n".format(
                    total_written, (float(total_written)/total_reads) * 100))
            f.write("{0} \t({1:.2f}%)\tSkipped\n".format(
                    total_skipped, (float(total_skipped)/total_reads) * 100))
            f.write("{0} \t({1:.2f}%) \tMIDtag or cutsite corrected\n".format(
                    total_corrected, (float(total_corrected)/total_reads) * 100))
            f.write("\n")
            f.write("{0} \t({1:.2f}%) \tMIDtags corrected in total\n".format(
                    total_MIDcorrected, (float(total_MIDcorrected)/total_reads) * 100))
            f.write("{0} \t({1:.2f}%) \tcutsites corrected in total\n".format(
                    total_cutsite_corrected, (float(total_cutsite_corrected)/total_reads) * 100))
            
        total_t = time.time() - toc    
        print 'Processed all files in {0}\n'.format(time.strftime('%H:%M:%S', 
                                                            time.gmtime(total_t)))
        # Update internal Variables
        self.next_input_files = outnames
        self.next_input_path = outpath
        
        return (outnames, outpath)

    def split_by_tags(self, infiles=None, inpath=None, outpath=None, out_filename=None,
                      report=True, savecounter=True):
        ''' Split the file into separate files based on MID tags '''
        
        c = self.c
        
        if outpath is None:
            outpath = c.tag_split_outpath
        
        if out_filename is None:
            out_filename = c.experiment_name
             
        # Setup Record Cycler        
        if infiles is None:
            infiles = self.next_input_files
        if inpath is None:
            inpath = self.next_input_path
        
        RecCycler = Cycler(infiles=infiles, filepattern=False, data_inpath=inpath)
         
        print ('\nSpliting {0} file(s) based on MID tags'
               '').format(RecCycler.numfiles)
        
        outfiles_dict = {}
        
        first_run = 1
        
        
        # Running through all records in all passed files 
        for recordgen in RecCycler.seqfilegen:
            
            # Set / reset Counter
            tag_counter = Counter()
            
            dbtags = self.get_data4file(RecCycler.curfilename, fields=['MIDtag', 'description'])
            # tags is returned as a list of tuples for each record            
            MID_length = len(dbtags[0][0])
            # as tuple of descriptions then tuple of MIDtags
            tups = zip(*dbtags)
            # Cheack using MIDtags as keys would be unique
            assert len(set(tups[1])) == len(tups[1]), 'Duplicate MIDtags returned for file {0}'.format(RecCycler.curfilename) 

            # Convert to dictionary  {'MIDtag': 'description'} 
            dbtags = dict(dbtags)
            
            # Open Files for Writing for each tag  
            for tag, desc in dbtags.iteritems():
                
                fname = '-'.join([out_filename, tag, desc]) + '.bgzf'
                fnamevar = 'f_' + desc

                # Check that files don't already exist
                if first_run:
                    # If file already exists, overwrite it.
                    if os.path.isfile(os.path.join(outpath, fname)):
                        f = open(os.path.join(outpath, fname), 'w')
                        f.close()
                    
                vars()[fnamevar] = bgzf.open(os.path.join(outpath, fname), 'a')
                outfiles_dict[fnamevar] = fname
    
            first_run = 0
            
            for rec in recordgen:
        
                recMIDtag = rec.seq[:MID_length].tostring()
                
                if recMIDtag not in dbtags:
                    raise Exception('MID tag not found in database for file {0}'.format(RecCycler.curfilename))
                else:                   
                    fnamevar = 'f_' + dbtags[recMIDtag]           
                    SeqIO.write(rec, vars()[fnamevar], 'fastq');
                    tag_counter[recMIDtag] += 1

            # Flush and Close Files for each tag  
            for tag, desc in dbtags.iteritems():

                fnamevar = 'f_' + desc
                vars()[fnamevar].flush()
                vars()[fnamevar].close()
                
                # Update datafiles in database
                filename = outfiles_dict[fnamevar]
                self.db.add_datafile(filename, [desc], datafile_type='1sample')

            print 'Finished Splitting MIDtags for input file: {0}'.format(RecCycler.curfilename)
            
            # Update counts
            for tag, desc in dbtags.iteritems():
                  
                row = self.db.select('''read_count FROM samples WHERE description=? ''', (desc,))
                current_value = row[0]['read_count']
                if current_value is None:
                    current_value = 0
                    
                self.db.update('''samples SET read_count=? WHERE description=?''',
                                ( current_value + tag_counter[tag], desc))

        # Store file names 
        for outfile in outfiles_dict.itervalues():

            # Find sample description            
            fname = os.path.split(outfile)[1]
            if fname.endswith('.bgzf'):
                fname = fname[:-5]
            fname_parts = fname.split('-') 
            desc = fname_parts[-1]
            
            self.db.update('''samples SET read_file=? WHERE description=?''',
                                (outfile, desc))
            
        # Outputs return / update next inputs
        self.next_input_path = outpath
        self.next_input_files = outfiles_dict.values()
        
        return (outfiles_dict.values(), outpath)
    
    def split_by_subgroups(self, subgroups=None, infiles=None, inpath=None, outpath=None, out_filename=None,
                      report=True, savecounter=True):
        ''' Split the file into separate files based on MID tags '''
        
        if subgroups is None:
            # Dictionary of regular expressions to match sample discription
            subgroups = { 'zebra'  : '.*zebra.*',
                         'gazelle' : '.*gazelle.*'}
        
        # Compile regexes
        for k,v in subgroups.iteritems():
            subgroups[k] = re.compile(v)
        
        c = self.c
        
        if outpath is None:
            outpath = c.tag_split_outpath
        if out_filename is None:
            out_filename = c.experiment_name
        
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        
        # Setup Record Cycler        
        if infiles is None:
            infiles = self.next_input_files
        if inpath is None:
            inpath = self.next_input_path
        
        RecCycler = Cycler(infiles=infiles, filepattern=False, data_inpath=inpath)
         
        print ('\nSpliting {0} file(s) into zebras and gazelles'
               '').format(RecCycler.numfiles)
        
        outfiles_dict = {}
        
        first_run = 1
        
        for recordgen in RecCycler.seqfilegen:
            
            # Set / reset Counter
            tag_counter = Counter()
            
            dbtags = self.get_data4file(RecCycler.curfilename, fields=['MIDtag', 'description'])
            # tags is returned as a list of tuples for each record            
            MID_length = len(dbtags[0][0])
            # Convert to dictionary  {'MIDtag' : 'description' }
            dbtags = dict(dbtags)
            
            # Open Files for Writing for each subgroup  
            for group in subgroups.iterkeys():
                
                fname = '-'.join([out_filename, group]) + '.bgzf'
                fnamevar = 'f_' + group

                # Check that files don't already exist
                if first_run:
                    # If file already exists, overwrite it.
                    if os.path.isfile(os.path.join(outpath, fname)):
                        f = open(os.path.join(outpath, fname), 'w')
                        f.close()
                    
                vars()[fnamevar] = bgzf.open(os.path.join(outpath, fname), 'a')
                outfiles_dict[fnamevar] = fname
    
            first_run = 0
            
            for rec in recordgen:
        
                recMIDtag = rec.seq[:MID_length].tostring()
                
                if recMIDtag not in dbtags:
                    raise Exception('MID tag not found in database for file {0}'.format(RecCycler.curfilename))
                else:
                    # Get description
                    desc = dbtags[recMIDtag]
                    # Write to approprite file if it matches the regex
                    for group in subgroups.iterkeys():
                        if subgroups[group].match(desc):
                            
                            fnamevar = 'f_' + group                
                            SeqIO.write(rec, vars()[fnamevar], 'fastq');
                            tag_counter[recMIDtag] += 1
                            
            # Flush and Close Files for each tag  
            for group in subgroups.iterkeys():

                fnamevar = 'f_' + group
                vars()[fnamevar].flush()
                vars()[fnamevar].close()
                
                # Update datafiles in database
                filename = outfiles_dict[fnamevar]
                
                desc_list = filter(subgroups[group].match ,dbtags.values())
                
                self.db.add_datafile(filename, desc_list, datafile_type='1sample')

            print 'Finished Splitting reads for input file: {0}'.format(RecCycler.curfilename)
            
        # Outputs return / update next inputs
        self.next_input_path = outpath
        self.next_input_files = outfiles_dict.values()
        
        return (outfiles_dict.values(), outpath)
                        
    def trim_reads(self, infiles=None, inpath=None, out_filename=None, outpath=None, n=1, mode='grouped'):
        ''' Trims off the MID tag of each read, as well as the last 'n' bases.
        
        mode    = grouped
                Writes all trimmed reads from the list of infiles to one large 
                output fasta file for clustering. Returns one file.
        
                = separate
                Writes all trimmed reads from the list of infiles to separate 
                fasta files for clustering. Returns a list of files. 
        '''
        
        c = self.c
        
        if outpath is None:
            if mode == 'grouped':
                outpath = c.tag_processed_outpath
            elif mode == 'separate':
                outpath = c.tag_split_outpath
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
        outfile_list = []
    
        print ('\nRemoving MID tags, trimming reads and converting {0} files to'
               ' fasta format\n').format(RecCycler.numfiles)
    
        # Generator to trim off MID tag and end of read.
        for seqfilegen in RecCycler.seqfilegen:
                       
            tags = self.get_data4file(RecCycler.curfilename, fields=['MIDtag'])
            # tags is returned as a list of tuples for each record            
            tags = zip(*tags)[0]
            # tags is now a tuple of all the first elements in each record              
            lenMIDs = len(tags[0])
            read_start_idx = len(c.cutsite) + lenMIDs  

            read_gen = (rec[read_start_idx:-n] for rec in seqfilegen)
            # File name
            
            fname = RecCycler.curfilename
            if fname.endswith('.bgzf') or fname.endswith('.fastq') or fname.endswith('.fasta'):
                outfile_i = '.'.join(RecCycler.curfilename.split('.')[:-1])
                
            outfile_i = outfile_i + '.fasta'
            outfile_list.append(outfile_i)
            count += 1

            with open(os.path.join(outpath, outfile_i), 'wb') as f:
                write_count = SeqIO.write(read_gen, f, 'fasta')
                print 'Wrote {0} records to file\n{1}'.format(write_count, outfile_i)

        if mode == 'grouped':
            # Combine output parts into one big file
            if outpath:
                os.chdir(outpath)
            cmd = ['cat'] + outfile_list 
            with open(os.path.join(outpath, out_filename), 'wb') as f:
                print 'Running "{0}" and saving to\n{1}'.format(cmd, os.path.join(outpath, out_filename))
                call(cmd, stdout=f) 
            print 'Done, cleaning up temp files ....'
            for f in outfile_list:
                os.remove(f)   
            os.chdir(start_dir)

            print '\nPreprocessing Complete.'
            self.next_input_files = out_filename
            self.next_input_path = outpath
            return (out_filename, outpath)
            
        elif mode == 'separate':
            
            print '\nPreprocessing Complete.'
            self.next_input_files = outfile_list
            self.next_input_path = outpath
            return (outfile_list, outpath)
                
    def cleanup_files(self, *args):
        ''' Remove intermediate files that are not needed '''      
        # Choices of 'filtered', 'tag_processed', 'all'
        # Raw inputs are always unchanged, and output .fasta is always kept
             
        files2remove = []
          
        if 'filtered' in args or 'all' in args:
            pattern = os.path.join(self.c.filtered_outpath, 
                                   '*' + self.c.filtered_files_postfix + ".fastq.bgzf") 
            files2remove.extend(glob.glob(pattern))
        elif 'tag_processed' in args or 'all' in args:
            pattern = os.path.join(self.c.tag_processed_outpath, 
                                   '*' + self.c.filtered_files_postfix + 
                                   self.c.tag_processed_files_postfix + 
                                   ".fastq.bgzf") 
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
        def illumina_filter(rec):
            ''' filter function for phred '''
            return rec.description.split()[1].split(':')[1] == 'N'
        return illumina_filter  
        
    def make_phred_filter(self, value):
        ''' Returns a filtering function based on the mean phred of each read.
        
        If a read has a mean phred of less than the given value, the filter returns
         false.   
        '''
        # Define filter function
        def phred_filter(rec):
            ''' filter function for phred '''
            return np.array(rec.letter_annotations['phred_quality']).mean() > value
        return phred_filter  
        
    def make_propN_filter(self, value):
        ''' Returns a filter function based on the proportion of Ns in the read.
    
        If a read has a higher proportion of Ns than the given value, the filter 
        returns false.
        '''
        # Define filter function
        def propN_filter(rec):
            ''' Filter function for propN'''
            return float(rec.seq.count('N')) / len(rec.seq) < value
        return propN_filter
                
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
        def cutsite_filter(rec):
            ''' Filter function for cutsite '''
            
            fname = self.current_file
            
            if cutsite_filter.target_file is None or cutsite_filter.target_file != fname:
                cutsite_filter.target_file = fname
                tags = self.get_data4file(fname, fields=['MIDtag'])
                cutsite_filter.MIDlength =  len(tags[0][0])
            
            cutsite = rec.seq[cutsite_filter.MIDlength: cutsite_filter.MIDlength + cutsite_length].tostring()
            cutsite_dist = ed.distance(target_cutsite, cutsite)
            
            return cutsite_dist <= max_edit_dist
        
        cutsite_filter.target_file = None
        cutsite_filter.MIDlength = None
        
        return cutsite_filter
        
    def make_overhang_filter(self, target_cutsite=None, overhang=None, max_edit_dist=0):
        ''' Returns a filter function based on the overhang part of the cutsite. 
        
        The cut site should end with the specified overhang. Those that dont are likely 
        to be genetic contaminants which have been inadvertantly sequenced, and 
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
        def overhang_filter(rec):
            ''' Filter function for cutsite'''
            
            # Must calculate MIDlength, but this may vary between files
            fname = self.current_file
#            fname = self.current_file.split('.')[0].split('-')[0]
            
            if overhang_filter.target_file is None or overhang_filter.target_file != fname:
                overhang_filter.target_file = fname
                tags = self.get_data4file(fname, fields=['MIDtag'])
                overhang_filter.MIDlength =  len(tags[0][0])
            
            cutsite = rec.seq[overhang_filter.MIDlength: overhang_filter.MIDlength + cutsite_length].tostring()
            if cutsite.endswith(overhang):
                return True
            else:
                overhang_dist = ed.distance(cutsite[-overhang_length:], overhang)
                return overhang_dist <= max_edit_dist
            
        overhang_filter.MIDlength = None
        overhang_filter.target_file = None
                  
        return overhang_filter

        
#===============================================================================
# Individual Functions 
#===============================================================================


if __name__ == '__main__':
    pass

