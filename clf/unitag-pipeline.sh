#!/bin/bash
#
# Example of how to run the various stages of the Salmon RAD-seq pipeline. 

# Summary specs of bens machine: 
#     198 GB RAM
#     32 core Intel(R) Xeon(R) CPU E5-2470 0 @ 2.30GHz


######################
# PARAMETERS DEFINES #
######################

# Unitag Reference creation
MIN_DEPTH=5
MAX_DEPTH=500

# Aligning stacks to reference and pstacks
STACK_MIN_DEPTH=2
CATALOGUE_MIN_DIST=0

# Num processors to run on.
NUM_THREADS=15


# File names and path locations 
RDFS_PATH="/home/musselle/salmon"     # home/musselle is my home on bens machine and the rdf storage is mounted with the name 'salmon'
PROJECT_ROOT="$RDFS_PATH/chris"       # Root directory for project on RDFS Storage. 
BIN_PATH="~/bin"                      # Where all python runscripts are located
STACKS_PATH="~bin/stacks/bin"         # Where all binaries for stacks are located. 


######################
# Data Preprocessing #
######################
# Assuming that stacks is installed and available on the system path. 
# The raw-data files in $PROJECT_ROOT/raw-data have already had the dummy barcodes added, 
# and are located under $PROJECT_ROOT/barcodes

echo "\nAbout to Run Preprocessing Steps"
# -p path to inputs directory 
# -b barcode file 
# -o
# -e Specify enzyme for cut site
# -q -c Filter based on phred quality (using sliding window) and discard reads with uncalled bases
# -r Correct RAD tags 
# --index_null Dummy barcode is in sequence header not in sequence itself. 
$BIN_PATH/process_radtags -p $PROJECT_ROOT/raw-data \                      
				-b $PROJECT_ROOT/barcodes/barcodes_only.txt  \
				-o $PROJECT_ROOT/processed-data/ \
				-e sbfI \     
				-q -c \      
				-r \          
				--index_null 
				
echo "\nPreprocessing Steps Complete"				 

#################################
# Unitag Reference Construction #
#################################

# Experiments with Multiple min max values (uncomment to run)
#------------------------------------------
# Uses GNU parallel to execute all combinations
#MIN_DEPTHS=2 5 10 15 20
#MAX_DEPTHS=500 
#parallel "$BIN_PATH/make_unitag.py -i $PROJECT_ROOT/processed-data/sample_GCAGGC.fq \
#                               -o $PROJECT_ROOT/processed-data/unitag/unitagref-m{1}-M{2}.fq \
#                               -m {1} -M {2}" ::: $MIN_DEPTHS ::: $MAX_DEPTHS
# Statistics for all runs are logged in unitag-logfile.log, though may be appended unordered.

echo "\nAbout to Construct Unitag Reference Sequence"
# -i input file to use as unitag ref 
# -o output file 
# -m min read depth threshold
# -M max read depth threshold
$BIN_PATH/make_unitag.py -i $PROJECT_ROOT/processed-data/sample_GCAGGC.fq \
            -o $PROJECT_ROOT/processed-data/unitag/unitagref-m$MIN_DEPTH-M$MAX_DEPTH.fq \
            -m $MIN_DEPTH -M $MAX_DEPTH

# Construct index using bowtie
# bowtie-build reads_file index_file
bowtie-build $PROJECT_ROOT/processed-data/unitag/unitagref-m$MIN_DEPTH-M$MAX_DEPTH.fq \
		     $PROJECT_ROOT/processed-data/unitag/unitag_idx-m$MIN_DEPTH-M$MAX_DEPTH
echo "\nBuilding Unitag Complete"

###########################
# Main Script: Unitag Ref #
###########################

# For each sample file, align reads to Unitag using bowtie, then run pstacks on the bowtie output
# -i Input path to processed reads 
# -o Output path to write files to
# -s Subset of samplet to run. run 'all' or select subset 'Alm' 'LoS' 'LoD' 'Tlt' 'UpS' or 'UpD' 
# -b barcodes file to map dummy barcodes to original file names 
# -q Starting index to use for MySQL database
# -x Index file 
# -m minimum number of matches to a locus to be considered a stack by pstacks.
# -p No. of threads to run with
# -k number of matches for bowtie to report per read. 
echo "\nAligning all samples to Unitag Reference and computing stacks"
$BIN_PATH/run_bowtie.py -i $PROJECT_ROOT/processed-data/ \       
						-s all \  
						-x $PROJECT_ROOT/processed-data/unitag/unitag_idx-m$MIN_DEPTH-M$MAX_DEPTH \
						-q 1 -m $STACK_MIN_DEPTH \
						-b $PROJECT_ROOT/barcodes/barcodes_filenames.txt \
						-o $PROJECT_ROOT/stacks/pstacks/ -p $NUM_THREADS -k 1
						--stackspath $STACKS_PATH
echo "\nFinished computing stacks"


echo "\nConstructing Global Catalogue"
# Construct Catalogue from all pstacks output
#   -i Input path to processed reads
#   -p  enable parallel execution with num_threads threads.
#   -b  MySQL ID of this batch.
#   -o  output path to write results.
#   -n  number of mismatches allowed between sample tags when generating the catalog.
$BIN_PATH/run_cstacks.py -i $PROJECT_ROOT/stacks/pstacks/* \ 
						 -o $PROJECT_ROOT/catalogues/unitag/ \
						 -b 200 \
						 -p $NUM_THREADS \
						 -n 0 \
						 --stackspath $STACKS_PATH 
						 
echo "\nMatching samples to Global Catalogue"















