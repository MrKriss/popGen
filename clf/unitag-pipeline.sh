#!/bin/bash
#
# Example of how to run the various stages of the Salmon RAD-seq pipeline. 

# Summary specs of bens machine: 
#     198 GB RAM
#     32 core Intel(R) Xeon(R) CPU E5-2470 0 @ 2.30GHz


######################
# PARAMETERS DEFINED #
######################

# Unitag Reference creation
MIN_DEPTH=5
MAX_DEPTH=500

# Aligning stacks to reference and pstacks
BOWTIE_MISMATCHES=3
STACK_MIN_DEPTH=2
CATALOGUE_MIN_DIST=0

# Batch Number
BATCH=200
BATCH2=400

# Num processors to run on.
NUM_THREADS=15

# File names and path locations 
RDFS_PATH="/home/musselle/salmon"     # home/musselle is my home on bens machine and the rdf storage is mounted with the name 'salmon'
PROJECT_ROOT="$RDFS_PATH/chris"       # Root directory for project on RDFS Storage. 
BIN_PATH="~/bin"                      # Where all python runscripts are located
STACKS_PATH="~bin/stacks/bin"         # Where all binaries for stacks are located. 
STACKS_OUTPUT="$PROJECT_ROOT/stacks/unitag"   # Where all output stacks files are written to by default 

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
						-o $STACKS_OUTPUT -p $NUM_THREADS -k 1
						--stackspath $STACKS_PATH
echo "\nFinished computing stacks"


echo "\nConstructing Global Catalogue"
# Construct Catalogue from all pstacks output
#   -i Input path to processed reads
#   -p  enable parallel execution with num_threads threads.
#   -b  MySQL ID of this batch.
#   -o  output path to write results.
#   -n  number of mismatches allowed between sample tags when generating the catalog.
$BIN_PATH/run_cstacks.py -i $PROJECT_ROOT/stacks/unitag/sample_* \ 
						 -o $STACKS_OUTPUT \
						 -b $BATCH \
						 -p $NUM_THREADS \
						 -n 0 \
						 --stackspath $STACKS_PATH 
						 
echo "\nMatching samples to Global Catalogue"
# Match all samples against the catalogue 
#  -s  Subset of samplet to run. run 'all' or select subset 'Alm' 'LoS' 'LoD' 'Tlt' 'UpS' or 'UpD' 
#  -b  barcodes file to map dummy barcodes to original file names
#  -P  Path to where all stacks output is located 
#  -c  File path to catalogue to use
#  -x  Starting Sql index/ batch number to use for output file
#  -p  Number of processors to run with
#  --stackspath  path to stacks binaries

$BIN_PATH/run_sstacks.py -s all \
						 -P $STACKS_OUTPUT \
						 -c $STACKS_OUTPUT/batch_$BATCH \
						 -o $STACKS_OUTPUT \
						 -b $PROJECT_ROOT/barcodes/barcodes_filenames.txt \
						 -x 400 
						 -p 15 \
						 --stackspath $STACKS_PATH 


###################################
# CALCULATE POPULATION STATISTICS #
###################################

# This section uses the populations program of stacks
 
# NOTE on population Map: the program requires a population map file defined as a tab delimited 
# list of individual vs the numerated group it belongs to e.g.
#    indiv_1    1
#    indiv_2    1
#    indiv_3    2
#    indiv_4    3
#    indiv_5    3
#
# This information in stored on the file:
#     $PROJECT_ROOT/barcodes/population_map.tsv
#
# With the mapping between this group number and the original subpopulations is stored in the file:
#     $PROJECT_ROOT/barcodes/population_map_dict.txt

# -b Batch number of catalogue to use for calculations 
# -M File path to where population map is


# Filtering parameters
#  -r    minimum percentage of individuals in a population required to process a locus for that population.
#  -p    minimum number of populations a locus must be present in to process a locus.
#  -m    specify a minimum stack depth required for individuals at a locus.
#  -a    specify a minimum minor allele frequency required before calculating Fst at a locus (0 < a < 0.5).
#  -f    specify a correction to be applied to Fst values: 'p_value', 'bonferroni_win', or 'bonferroni_gen'.
#  --p_value_cutoff [num] — required p-value to keep an Fst measurement (0.05 by default). Also used as base for Bonferroni correction.

#  --vcf  Output results in a vcf file.
 
$STACKS_PATH/populations -b $BATCH \
						 -M $PROJECT_ROOT/barcodes/population_map.tsv \
						 --vcf


#############################
# Load data into a database #
#############################

# Load in data from this pipeline into a MySQL database 

$STACKS_PATH/load_radtags.pl 








						 








