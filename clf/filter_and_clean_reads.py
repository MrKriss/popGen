'''
Created on 3 Jul 2013

@author: musselle
'''
import os 
import sys 
import glob
import time

import numpy as np 
import matplotlib.pyplot as plt

import editdist

from utils.preprocess import Preprocessor, ConfigClass
from utils.fileIO import SeqRecCycler
from collections import Counter

from Bio import SeqIO, bgzf 
import gzip

import argparse



class RecordPreprocessor(object):
    
    def __init__(self, args):
        
        # Set counting variables 
        self.filterfail_counter = Counter()
        self.total_read_count = 0
        
        self.skipped_count = 0
        self.MIDtag_corrected_count = 0
        self.cutsite_corrected_count = 0
        self.read_corrected_count = 0
            
        # Load barcodes 
        if 'barcodes' in args:
            # input checks
            barcode_files = glob.glob(args.barcodes)
            assert barcode_files, "No barcode files found at destination"
            barcode_files.sort()

            # Store barcode dictionary            
            MIDs = []
            individuals = []
            for barcode_file in barcode_files:
                with open(barcode_file, 'rb') as f: 
                    for line in f:
                        parts = line.strip().split()
                        MIDs.append(parts[0])
                        individuals.append(parts[1])
                    diff = len(MIDs) - len(set(MIDs))
                    if diff > 0:
                        raise Exception('MID tags in barcode files are not unique.\n{0} duplicates found'.format(diff))
                    else:
                        self.barcode_dict = dict(zip(MIDs, individuals))
                        
        # Check length of MIDs
        lengths = [len(mid) for mid in self.barcode_dict.keys()]
        assert all([x == lengths[0] for x in lengths]), 'Error: Lengths of MIDs are not equal'
        self.MIDlength = lengths[0]

        # Setup filter functions 
        self.filter_functions = []
        
        if 'fi' in args:
            # Ilumina Machine filtering 
            self.filter_functions.append(self.set_illumina_filter())
        if 'fn' in args:
            # Filter by proportion of Ns
            self.filter_functions.append(self.set_propN_filter(args.fn))
        if 'fp' in args:
            # Filter by mean phred score
            self.filter_functions.append(self.set_phred_filter(args.fp))
        if 'fc' in args:
            # Filter by target Cutsites
            assert 'ed' in args, 'Edit distance must be specified to filter by cutsite'
            self.cutsites = args.fc
            
            # Check all lengths are equal
            cutsite_lengths = [len(cutsite) for cutsite in self.cutsites]
            assert all([x == cutsite_lengths[0] for x in cutsite_lengths]), 'Error: Lengths of cutsites are not equal'
            self.cutsite_length = cutsite_lengths[0]

            self.filter_functions.append(self.set_cutsite_filter(target_cutsites=self.cutsites, 
                                                                 max_edit_dist=args.ed, 
                                                                 midtag_length=self.MIDlength))
        if 'fo' in args:
            # Filter by target Cutsites
            assert 'ed' in args, 'Edit distance must be specified to filter by overhang'
            self.filter_functions.append(self.set_overhang_filter(target_cutsites=args.fc,
                                                                    overhang=args.fo, 
                                                                    max_edit_dist=args.ed, 
                                                                    midtag_length=self.MIDlength)) 
        # Setup Error correction 
        if 'ted' in args:
            self.error_corrected_dist = args.ted
                     
            
        
    def runfilters(self, rec):
        ''' Filter record based on criteria specified by a sequence of functions given
         in self.filter_functions.
        
        Each functions in self.filter_functions take in a sequence record object, 
        and returns a boolean.
        
        Only records that pass all filters are returned as True
          
        '''
        pass
    
    def write_summary_output(self, path):
        ''' Write the results of the filtering and cleaning to a summary output file'''
    
        filepath = os.path.join(path, "filtering_summary.log")
        
        if os.exists(filepath):
            var = raw_input("Summary file already exists. Overwrite? [y/n]")
            if var == 'y':
                file_handle = open(filepath, 'wb'):
            else:
                count = 0 
                while os.exists(filepath):
                    count += 1
                    filepath = os.path.join(path, "filtering_summary%.log" % str(count))
                
                
                file_handle = open(filepath, 'wb')
    
        # Write the summary to a file 
        with open(os.path.join(path, "filtering_summary_" + self.c.root_name + 
                               '_params-' + str(self.c.filterparam_id) + ".log"), 'wb') as f:
            f.write("Filter parameters:\n") 
            f.write("------------------\n")
            p = self.db.get_binary(col='params', target='filtering_parameterID',
                                    value = c.filterparam_id, table = 'filtering_parameters')
            f.write(str(p) + "\n")
            f.write("\n")
            f.write('Filter stats:\n')
            f.write("-------------------------------------------\n")
            f.write('Filter\t\t\tHits\n')    
            f.write("-------------------------------------------\n")
            for i, x in enumerate(self.filter_functions):
                percent = preprocessor.filterfail_counter[i] / float(sum(preprocessor.filterfail_counter.values())) 
                f.write('%s\t\t%s\t(%.2f %%)\n' % (x.__name__, preprocessor.filterfail_counter[i], percent * 100))
            f.write('\nTotal No. Reads Processed:  \t%s\n' % preprocessor.total_read_count)
            f.write('Total No. filtered:  \t\t%s (%.2f %%)\n' % (sum(preprocessor.filterfail_counter.values()), 
                                                        100 * (sum(preprocessor.filterfail_counter.values())/ 
                                                               float(preprocessor.total_read_count))))
            f.write('Total No. Passed: \t\t%s (%.2f %%)\n' % (preprocessor.total_read_count - sum(preprocessor.filterfail_counter.values()), 
                                                        100 * ((preprocessor.total_read_count - sum(preprocessor.filterfail_counter.values()))/ 
                                                               float(preprocessor.total_read_count))))
            self.db.update('filtering_parameters SET filtering_summary=? WHERE filtering_parameterId=?', 
                           (os.path.split(f.name)[1], c.filterparam_id))
            
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
    
    


    def make_processing_gen(self, recgen):
        ''' A generator to only yield records that pass all filter functions
         specified
         
         Plus if error_corrected_dist argument is set, the MIDtag and cutsite will 
         be error corrected up this number of edits. 
         '''
        
        filterfuncs = self.filter_functions
        
        MIDtagslist = self.barcode_dict.keys()
        
        goto_next_rec = False
        
        for rec in recgen:
        
            for n, func in enumerate(filterfuncs):
                #===============================================================
                # Run filter functions on record       
                #===============================================================
                if func(rec) == False:
                    self.filterfail_counter.update([n])
                    self.total_read_count += 1 
                    goto_next_rec = True
                    break # move on to next record and don't do else:
            
            if goto_next_rec:
                goto_next_rec = False
                continue
                
            else: # if runs through all filterfuncs and pased, run this too. 
                self.total_read_count += 1 
                    
                # RECORD HAS PASSED ALL FILTERS
            
                if hasattr('error_corrected_dist', self):
                    
                    MID_corrected = False
            
                    #===============================================================
                    # Analise MID tag section
                    #===============================================================
                    recMID = str(rec.seq[:self.MIDlength])
                    
                    if recMID not in self.barcode_dict:
                        # Sequencing error in the tag. Work out nearest candidate.
                        distvec = [editdist.distance(recMID, tag) for tag in MIDtagslist]
                        
                        distvec_min = min(distvec)
                        count_distvec_min = distvec.count(distvec_min)
                        
                        if distvec_min > self.error_corrected_dist:
                            self.skipped_count += 1
                            continue
                            
                        count_distvec_min = distvec.count(distvec_min)
                        if count_distvec_min > 1:
                            # Muliple candidates. True MID is Ambiguous 
                            self.skipped_count += 1
                            continue
                    
                        elif count_distvec_min == 1:
                            
                            correct_MIDtag = MIDtagslist[distvec.index(distvec_min)]
                            
                            # Correct the erroneous tag with the candidate.
                            # Letter annotations must be removed before editing rec.seq
                            temp_var = rec.letter_annotations
                            rec.letter_annotations = {}
                            # Change seq to mutableseq
                            rec.seq = rec.seq.tomutable()
                            rec.seq[:self.MIDlength] = correct_MIDtag
                            rec.seq = rec.seq.toseq()
                            rec.letter_annotations.update(temp_var)
                            self.MIDtag_corrected_count += 1
                            self.read_corrected_count += 1
                            MID_corrected = True
                            
                    #===============================================================
                    # Analise Cut site section
                    #===============================================================
                    # Must allow for multiple possible cutsites
                    rec_cutsite = str(rec.seq[self.MIDlength: self.MIDlength + len(self.cutsite_length)])
                    if rec_cutsite not in self.cutsites:
                        # Sequencing error in the cutsite. Correct if less than max_edit_dist
                        
                        cutsite_dists = [editdist.distance(rec_cutsite, cutsite) for cutsite in self.cutsites]
    
                        min_cutsite_dist = min(cutsite_dists)
                        
                        if min_cutsite_dist > self.error_corrected_dist:
                            # Amount of error in cut site is too large to correct 
                            # May also be from contaminants. So read is skipped altogether 
                            self.skipped_count += 1
                            if MID_corrected:
                                # Roll back counters
                                self.MIDtag_corrected_count -= 1
                                self.read_corrected_count -= 1
                            continue
                        
                        min_cutsite_dist_count = cutsite_dists.count(min_cutsite_dist)
                        
                        if cutsite_dists.count(min_cutsite_dist_count) > 1:
                            # Muliple candidates. True cutsite is Ambiguous 
                            self.skipped_count += 1
                            if MID_corrected:
                                # Roll back counters
                                self.MIDtag_corrected_count -= 1
                                self.read_corrected_count -= 1
                            continue
                            
                            
                        elif cutsite_dists.count(min_cutsite_dist_count) == 1:
                            
                            corrected_cutsite = self.cutsites[cutsite_dists.index(min_cutsite_dist)]
                            
                            # Correct the erroneous cutsite with the actual cutsite.
                            # Letter annotations must be removed before editing rec.seq
                            temp_var = rec.letter_annotations
                            rec.letter_annotations = {}
                            # Change seq to mutableseq
                            rec.seq = rec.seq.tomutable()
                            rec.seq[self.MIDlength: self.MIDlength + len(self.cutsite)] = corrected_cutsite
                            rec.seq = rec.seq.toseq()
                            rec.letter_annotations.update(temp_var)
                            self.cutsite_corrected_count += 1
                            if not MID_corrected:
                                self.read_corrected_count += 1
                        else:
                            # Amount of error in cut site is too large to correct 
                            # May also be from contaminants. So read is skipped altogether 
                            self.skipped_count += 1
                            if MID_corrected:
                                # Roll back counters
                                self.MIDtag_corrected_count -= 1
                                self.read_corrected_count -= 1
                            continue
        
                # Note: rec only yielded if the read passes all filters, and if specified,
                # only if the MIDtag and cutsites are cleaned sucessfully.   
                yield rec
                

    def set_illumina_filter(self):
        ''' Returns filtering function based on illumina machine filter         
        N = was not pickup up by machine filter i.e. passed
        Y = was flagged by machine filter i.e. fail 
        '''
        # Define filter function
        def illumina_filter(rec):
            ''' filter function for phred '''
            return rec.description.split()[1].split(':')[1] == 'N'
        return illumina_filter  
    
    def set_propN_filter(self, value):
        ''' Returns a filter function based on the proportion of Ns in the read.
    
        If a read has a higher proportion of Ns than the given value, the filter 
        returns false.
        '''
        # Define filter function
        def propN_filter(rec):
            ''' Filter function for propN'''
            return float(rec.seq.count('N')) / len(rec.seq) < value
        return propN_filter
    
    def set_phred_filter(self, value):
        ''' Returns a filtering function based on the mean phred of each read.
        
        If a read has a mean phred of less than the given value, the filter returns
         false.   
        '''
        # Define filter function
        def phred_filter(rec):
            ''' filter function for phred '''
            return np.array(rec.letter_annotations['phred_quality']).mean() > value
        return phred_filter  
    
    def set_cutsite_filter(self, target_cutsites=None, max_edit_dist=None, midtag_length=None):
        ''' Returns a filter function based on the match of the read cutsite to the 
        target_cutsite given.
        
        Reads that differ in edit distance by more than mindist, cause filter to 
        return false 
        
        '''
        cutsite_length = len(target_cutsites[0])
    
        # Define filterfunc
        def cutsite_filter(rec):
            ''' Filter function for cutsite '''
            
            cutsite = rec.seq[midtag_length: midtag_length + cutsite_length].tostring()
            
            for target_site in target_cutsites:
                cutsite_dist = editdist.distance(target_site, cutsite)
                if cutsite_dist <= max_edit_dist:
                    return True
                
            return False
        
        return cutsite_filter
        
    def set_overhang_filter(self, target_cutsites=None, overhang=None, midtag_length=None, max_edit_dist=0):
        ''' Returns a filter function based on the overhang part of the cutsite. 
        
        The cut site should end with the specified overhang. Those that dont are likely 
        to be genetic contaminants which have been inadvertantly sequenced, and 
        therefore should be discarded. 
           
        Reads that mismatch in the overhang region by more than mindist, cause the 
        filter to return false. 
        
        '''   
        overhang_patterns = [target_cutsites[-overhang:] for i in target_cutsites]
        cutsite_length = len(target_cutsites[0])
        
        # Define filterfunc
        def overhang_filter(rec):
            ''' Filter function for cutsite'''
            
            cutsite = rec.seq[midtag_length: midtag_length + cutsite_length].tostring()

            for i, pat in enumerate(overhang_patterns):
                
                dist = editdist.distance(target_cutsites[i], cutsite)
                if dist <= max_edit_dist:
                    if cutsite.endswith(pat):
                        return True
            
            return False
            
        return overhang_filter













#===============================================================================

    
parser = argparse.ArgumentParser(description='Filter and clean up FastQ files.')
parser.add_argument('-i', '--input', action='store', type=int, dest='input',
                    help='Input file(s) to process. (/path/filename) Will accept a glob')
parser.add_argument('-o', '--output_postfix', action='store', type=int, dest='output_postfix',
                    help='Output file postfix to use when writing to file.')
parser.add_argument('-p', '--output_path', action='store', type=int, dest='output_path', default='',
                    help='Output path to write reads to. Missing dirs will be created.')
parser.add_argument('-b', '--barcodes', action='store', type=int, dest='barcodes',
                    help='Barcodes accociated with input file(s). Will accept a glob')

parser.add_argument('--version', action='version', version='%(prog)s 1.0')


# Filter Parameters
parser.add_argument('-d', '--filterdefaults', action='store_true', dest='fd',
                    help='Use default values for all filters. n = 0.1, p = 20, l=True, c=TGCAGG, e=2, r=2, f=1')
parser.add_argument('-n', '--filterNs', action='store', dest='fn', type=float,
                    help='Threshold maximum for filtering by proportion of Ns in the read. Default = 0.1. Set to 0 to skip.')
parser.add_argument('-p', '--filterphred', action='store', dest='fp', type=int,
                    help='Threshold minimum for filtering by mean phred of the read. Default = 20. Set to 0 to skip.')
parser.add_argument('-l', '--filterillumina', action='store_true', dest='fi', default=False, 
                    help='Filter out reads that failed the Illumina machine filter. Default = True.')

parser.add_argument('-c', '--filtercutsite', action='append', dest='fc',
                    help='Filter reads that do not have one of the cutsites specified. May be used repeatedly for multiple cutsites.')
parser.add_argument('-e', '--filtercutsite_editdist', action='append', dest='ed',
                    help='Max edit distance allowed between target cutsite and read.')
parser.add_argument('-r', '--filteroverhang', action='store', dest='fo', type=int, 
                    help=('Number of bases in cutsite that make up the overhang. Reads are filtered out which' 
                    'have errors in the overhang of the cute site.'))

# parser.add_argument('-m', '--midtaglength', action='append', dest='ml',
#                     help='Length of MIDtag for all reads. Needed to filter by cutsite and overhang')

parser.add_argument('-f', '--corrected_editdist', action='append', dest='ted',
                    help='Max edit distance that is corrected between target MIDtag/cutsite and actual read.'
                    'If matched to more than one candidate barcode, the read is discarded due to abiguity of identity.')

# For displaying filter results to stderr 
parser.add_argument('-v', '--verbose', action='store_true', default=False, 
                    help='Whether to log the Sequence IDs of reads that fail the filtering')

# Cleaning parameters
parser.add_argument('-r', '--rootdir', action='store', help='Root directory where all ')
parser.add_argument('-n', '--name', action='store', type=int, help='Name of the experiemnt, and subdirectory where data is stored.')


# Parse args and set defaults 
args = parser.parse_args()

if 'fd' in args:
    default_args = []
    if 'fn' not in args:
        default_args.extend(['-n', '0.1'])
    if 'fp' not in args:
        default_args.extend(['-p', '20'])
    if 'fi' not in args:
        default_args.extend(['-l'])
    if 'fc' not in args:
        default_args.extend(['-c', 'TGCAGG', '-e', '2', ])
    if 'fo' not in args:
        default_args.extend(['-r', '2'])
    if 'ted' not in args:
        default_args.extend(['-f', '1'])
        
    args = parser.parse_args(default_args)




starting_dir = os.getcwd()
# Input checks 
if args.input == "-":    
    input_handle = sys.stdin
else:
    if type(args.input) is str:
        # Glob the file pattern 
        files = glob.glob(args.input)
        assert files, "No files returned for specified input"
        files.sort()
    else:
        raise Exception('Invalid entry for data_files.')
        

# Generator to cycle through files
reads_generator = SeqRecCycler(data_files=files)

# Initalise Class to Filter and cleanup reads
preprocessor = RecordPreprocessor(args)

# Output Path and file variables
if 'output_path' in args:
    outputnames = []
    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path) 

cum_t = 0
toc = time.time()


for recgen in reads_generator.seqfilegen:

    print >> sys.stderr, '\nProcessing {0}'.format(reads_generator.curfilename)
    
    # Initialise generator
    passes = preprocessor.make_processing_gen(recgen)
    
    for rec in passes:
        # These records have passed filter and been cleaned
        
        # DEfine output location
        # Output Path and file variables
        if 'output_postfix' in args:

            # Construct file names
            name = reads_generator.curfilename.split('.')  
            pass_filename = '.'.join([name[0] + args.output_postfix] + name[1:]) 
            pass_file = os.path.join(args.output_path, pass_filename)
            name = '.'.join(name)
            outputnames.append(pass_filename) 
            
            # Writes to the same filetype as the input
            if name.endswith('.bgzf'):
                pass_filehdl = bgzf.BgzfWriter(pass_file)
            elif name.endswith('.fastq'):
                pass_filehdl = open(pass_file, 'wb')
            elif name.endswith('.gz'):
                pass_filehdl = gzip.open(pass_file, 'wb')
            else:
                print >> sys.stderr, 'Input file format not supported: %s' % name
                sys.exit()
        
        else:
            # Output is written to std out
            pass_filehdl = sys.stdout
                     
        print >> sys.stderr, 'Writing passes to \n{0} ....'.format(pass_filename)
        numwritten = SeqIO.write(passes , pass_filehdl , 'fastq')
        
        if pass_filehdl == sys.stdout:
            pass_filehdl.flush()
        else:
            pass_filehdl.close()
        print >> sys.stderr, '{0} records Preprocessed'.format(numwritten)
        
        loop_t = time.time() - toc - cum_t
        cum_t += loop_t
        print >> sys.stderr, 'Finished file {0} after {1}'.format(reads_generator.curfilenum, 
                                    time.strftime('%H:%M:%S', time.gmtime(loop_t))) 
        
        
        if args.verbose:
            print >> sys.stderr, '\nFilter stats'
            print >> sys.stderr, '\nFilter No.\tHits'    
            for n in range(len(preprocessor.filter_functions)):
                print >> sys.stderr, '%s\t\t%s' % (n,preprocessor.filterfail_counter[n])
        
            print >> sys.stderr, '\nTotal No. Reads Processed:  %s' % preprocessor.total_read_count
            print >> sys.stderr, '\nTotal No. filtered:  %s (%.2f %%)' % (sum(preprocessor.filterfail_counter.values()), 
                                                            100 * (sum(preprocessor.filterfail_counter.values())/ 
                                                                    float(preprocessor.total_read_count)))
        
        preprocessor.write_summary_output(path)
        
        
        
        
        


        # Clean 
    

        
        # Write 


#===============================================================================















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
            if hasattr(self.c, 'starting_input_files'):
                self.next_input_files = c.starting_input_files
            if hasattr(self.c, 'starting_data_inpath') :
                self.next_input_path = c.starting_data_inpath   
        else:
            self.c = ConfigClass()
            # Defaults for other parameters
            self.filter_functions = None 
        
    
    def filter_reads_pipeline(self, infiles=None, inpath=None):
    
        ''' Filter reads based on criteria specified by a sequence of functions given
         in self.filter_functions.
        
        Individual functions in self.filter_functions must take in a sequence record object, 
        and return a boolean.
        
        differs from filter_reads in that this does not save records that fail a filter,
        only those that pass, but if log_fails = True, will log all reads that fail 
        any of the filter functions given, along with the filter number it failed at:
        logfile example:
        
        (rec.id) <TAB> (enumeration of filter failed from list given)  
        e.g.read_count
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
        with open(os.path.join(outpath, "filter_summary_" + self.c.root_name + 
                               '_params-' + str(self.c.filterparam_id) + ".log"), 'wb') as f:
            f.write("Filter parameters:\n") 
            f.write("------------------\n")
            p = self.db.get_binary(col='params', target='filtering_parameterID',
                                    value = c.filterparam_id, table = 'filtering_parameters')
            f.write(str(p) + "\n")
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
            self.db.update('filtering_parameters SET filtering_summary=? WHERE filtering_parameterId=?', 
                           (os.path.split(f.name)[1], c.filterparam_id))
            
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
            
            tags = self.get_data4file(RecCycler.curfilename, fields=['MIDtag', 'description'])
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
        with open(os.path.join(outpath, "cleaner_summary_" + self.c.root_name + 
                               '_params-' + str(self.c.filterparam_id) + ".log"), 'wb') as f:
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
            
            self.db.update('filtering_parameters SET cleaning_summary=? WHERE filtering_parameterId=?', 
                           (os.path.split(f.name)[1], c.filterparam_id))
            
            
        total_t = time.time() - toc    
        print 'Processed all files in {0}\n'.format(time.strftime('%H:%M:%S', 
                                                            time.gmtime(total_t)))

        # Update datafiles in database
        for fname in outnames:
            # find samples present in file 
            recs = self.get_data4file(fname, fields=['description'])
            # Store results in database 
            desc_list = [rec[0] for rec in recs]
            data_id = self.db.add_datafile(fname, desc_list, datafile_type='processed')
            
            # UPDATE table_name SET column1=value, column2=value,... WHERE some_column=some_value
            self.db.update('datafiles SET filtering_parameterId=? WHERE datafileId=?', (self.c.filterparam_id, data_id))
        
        # Update internal Variables
        self.next_input_files = outnames
        self.next_input_path = outpath
        
        return (outnames, outpath)



if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Filter and clean up FastQ files.')
    parser.add_argument('filename', action='store', help='Name of fastq input files. Can be a Glob')
    parser.add_argument('min', action='store', type=int, help='Minimum value of threshold for cluster size')
    parser.add_argument('max', action='store', type=int, help='Maximum value of threshold for cluster size')
    parser.add_argument('minreads', action='store', type=int, help='Minimum value of threshold for total reads in clusters of any size')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    
    
    
    c = ConfigClass()
    
    
    
    
    
    P = Preprocessor()
    
    
    
    
    P.

    
    
    