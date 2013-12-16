#!/bin/bash
#
# Example of how to run the various stages of the Salmon RAD-seq pipeline. 
#
########################################################################
# FULL PIPELINE: 
#       Preprocessing and Forming the Stacks.
#       Then Catalogue Construction, Matching and Populations Statistics
########################################################################

# Summary specs of bens machine: 
#     198 GB RAM
#     32 core Intel(R) Xeon(R) CPU E5-2470 0 @ 2.30GHz

######################
# PARAMETERS DEFINED #
######################

# Unitag Reference creation
MIN_DEPTH=5
MAX_DEPTH=500
UNITAG_REF="LoD06_"

# Aligning stacks to reference and pstacks
BOWTIE_MISMATCHES=3
STACK_MIN_DEPTH=5

# Unique Batch Number for this run 
BATCH_ID=1
BATCH_NAME="unitag"

# Aligning stacks to reference and pstacks
CATALOGUE_MIN_DIST=0

# Num processors to run on.
NUM_THREADS=15

# Population Analysis Parameters 
#-------------------------------

# Three sets of parameters are run inline with Floragenex's analysis: 
# Relaxed (RLX), Standard (STD), and Stringent (STR)

#  MIN_PERC     minimum percentage of individuals in a population required to process a locus for that population. 
#  MIN_DEPTH    specify a minimum stack depth required for individuals at a locus.
#  MIN_AF       specify a minimum minor allele frequency required before calculating Fst at a locus (0 < a < 0.5).
MIN_DEPTH_RLX=6
MIN_AF_RLX=0.1
MIN_PERC_RLX=0.60

MIN_DEPTH_STD=10
MIN_A_STDF=0.1
MIN_PERC_STD=0.80

MIN_DEPTH_STR=15
MIN_AF_STR=0.1
MIN_PERC_STR=0.90

# File names and path locations 
RDFS_PATH="/home/musselle/salmon"     # home/musselle is my home on bens machine and the rdf storage is mounted with the name 'salmon'
PROJECT_ROOT="$RDFS_PATH/chris"       # Root directory for project on RDFS Storage. 
BIN_PATH="$HOME/bin"                      # Where all python runscripts are located
STACKS_PATH="$HOME/bin/stacks/bin"         # Where all binaries for stacks are located. 
STACKS_OUTPUT="$PROJECT_ROOT/stacks/$BATCH_NAME"   # Where all output stacks files are written to by default 

DATABASE_NAME="$PROJECT_ROOT/salmon_database.db"

mkdir -p $STACKS_OUTPUT
mkdir -p "$STACKS_OUTPUT/rlx"
mkdir -p "$STACKS_OUTPUT/std"
mkdir -p "$STACKS_OUTPUT/str"
mkdir -p "$STACKS_OUTPUT/default"

######################
# Data Preprocessing #
######################

# Assuming that stacks is installed and available on the system path. 
# The raw-data files in $PROJECT_ROOT/raw-data have already had the dummy barcodes added, 
# and are located under $PROJECT_ROOT/barcodes

# Adding Dummy barcodes to Flourogenex Data
#-------------------------------------------

# Generate the barcode to use and store them in files:
#     barcodes_only.txt  -->  List of barcode one per line (for input to stacks preprocessing)
#     barcodes_filenames.txt  -->   List of barcode, file name pairs. One per line. 

# uncomment to run (barcodes already generated, will override if reran and require all subsequent analysis to be rerun)
# echo "\nAbout to Regenerate and re-append barcodes"
# $BIN_PATH/makebarcodes.py -i $PROJECT_ROOT/raw-data \
#                          -o $PROJECT_ROOT/barcodes/  \
# Add the generated barcodes to the raw reads files. 
#ls $PROJECT_ROOT/raw-data/* | parallel -j $NUM_THREADS "$BIN_PATH/add_bar2header.py -i {} \
#                                                                                    -o {} \
#                                                                                    -b $PROJECT_ROOT/barcodes/barcodes_filenames.txt"

echo "\nAbout to Run Preprocessing Steps"
# -p path to inputs directory 
# -b barcode file 
# -o Output path to save files to 
# -e Specify enzyme for cut site
# -q -c Filter based on phred quality (using sliding window) and discard reads with uncalled bases
# -r Correct RAD tags 
# --index_null Dummy barcode is in sequence header not in sequence itself. 

# -s The value for the sliding window average quality thresholding (default 10) 
# -w Set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).

$STACKS_PATH/process_radtags -p $PROJECT_ROOT/raw-data \
-b $PROJECT_ROOT/barcodes/barcodes_only.txt \
-o $PROJECT_ROOT/processed-data/ \
-e sbfI -q -c -s 15 -r --index_null
			
# rename stacks output files to contain original filename id 
$BIN_PATH/update_filenames.py -i $PROJECT_ROOT/processed-data/*.fq \
-b $PROJECT_ROOT/barcodes/barcodes_filenames.txt 

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
$BIN_PATH/make_unitag.py -i $PROJECT_ROOT/processed-data/$UNITAG_REF*.fq \
-o $PROJECT_ROOT/processed-data/unitag/unitagref-m$MIN_DEPTH-M$MAX_DEPTH.fq \
-m $MIN_DEPTH -M $MAX_DEPTH

# Construct index using bowtie
# bowtie-build reads_file index_file
bowtie-build $PROJECT_ROOT/processed-data/unitag/unitagref-m$MIN_DEPTH-M$MAX_DEPTH.fq \
$PROJECT_ROOT/processed-data/unitag/unitag_idx-m$MIN_DEPTH-M$MAX_DEPTH
echo "\nBuilding Unitag Complete"

###################
# STACK FORMATION #
###################

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
$BIN_PATH/run_bowtie.py -i $PROJECT_ROOT/processed-data/*.fq \
-x $PROJECT_ROOT/processed-data/unitag/unitag_idx-m$MIN_DEPTH-M$MAX_DEPTH \
-a $PROJECT_ROOT/alignments/ \
-q 1 -m $STACK_MIN_DEPTH \
-o $STACKS_OUTPUT -p $NUM_THREADS -k 1 \
--stackspath $STACKS_PATH
echo "\nFinished computing stacks"


#####################################
# CONSTRUCTING CATALOGUE & MATCHING #
#####################################

echo "\nConstructing Global Catalogue"
# Construct Catalogue from all pstacks output
#   -i Input path to processed reads
#   -p  enable parallel execution with num_threads threads.
#   -b  MySQL ID of this batch.
#   -o  output path to write results.
#   -n  number of mismatches allowed between sample tags when generating the catalog.
$BIN_PATH/run_cstacks.py -i $STACKS_OUTPUT/*snps.tsv \
-o $STACKS_OUTPUT \
-b $BATCH_ID \
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
-c $STACKS_OUTPUT/batch_$BATCH_ID \
-o $STACKS_OUTPUT \
-b $PROJECT_ROOT/barcodes/barcodes_filenames.txt \
-x 1 \
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

# Optional Filtering parameters
#  -r    minimum percentage of individuals in a population required to process a locus for that population.
#  -p    minimum number of populations a locus must be present in to process a locus.
#  -m    specify a minimum stack depth required for individuals at a locus.
#  -a    specify a minimum minor allele frequency required before calculating Fst at a locus (0 < a < 0.5).
#  -f    specify a correction to be applied to Fst values: 'p_value', 'bonferroni_win', or 'bonferroni_gen'.
#  --p_value_cutoff [num] — required p-value to keep an Fst measurement (0.05 by default). Also used as base for Bonferroni correction.

#  --vcf  Output results in a vcf file.

# To make population map file again, uncomment the following
#$BIN_PATH/make_popmap.py -i $PROJECT_ROOT/processed-data/*.fq \
#						 -o $PROJECT_ROOT/barcodes/population_map.tsv \
#						 -b $PROJECT_ROOT/barcodes/barcodes_filenames.txt


# Three sets of parameters are run for the RLX, STD and STR settings,
#  and the output files transfered to a suitable subdirectory (to save output summary files from being overwritten)

echo "\nAbout to Calculate Population Statistics for Relaxed Settings" 
$STACKS_PATH/populations -b $BATCH_ID \
-P $STACKS_OUTPUT \
-M $PROJECT_ROOT/barcodes/population_map.tsv \
--vcf \
-t $NUM_THREADS \
-m $MIN_DEPTH_RLX \
-a $MIN_AF_RLX \
-r $MIN_PERC_RLX

mv $STACKS_OUTPUT/batch_1.{fst*,hap*,pop*,sum*,vcf} $STACKS_OUTPUT/rlx/

echo "\nAbout to Calculate Population Statistics for Standard Settings" 
$STACKS_PATH/populations -b $BATCH_ID \
-P $STACKS_OUTPUT \
-M $PROJECT_ROOT/barcodes/population_map.tsv \
--vcf \
-t $NUM_THREADS \
-m $MIN_DEPTH_STD \
-a $MIN_AF_STD \
-r $MIN_PERC_STD

mv $STACKS_OUTPUT/batch_1.{fst*,hap*,pop*,sum*,vcf} $STACKS_OUTPUT/std/

echo "\nAbout to Calculate Population Statistics for Stringent Settings" 
$STACKS_PATH/populations -b $BATCH_ID \
-P $STACKS_OUTPUT \
-M $PROJECT_ROOT/barcodes/population_map.tsv \
--vcf \
-t $NUM_THREADS \
-m $MIN_DEPTH_STR \
-a $MIN_AF_STR \
-r $MIN_PERC_STR

mv $STACKS_OUTPUT/batch_1.{fst*,hap*,pop*,sum*,vcf} $STACKS_OUTPUT/str/

echo "\nAbout to Calculate Population Statistics for Default Settings" 
$STACKS_PATH/populations -b $BATCH_ID \
-P $STACKS_OUTPUT \
-M $PROJECT_ROOT/barcodes/population_map.tsv \
--vcf \
-t $NUM_THREADS \
-m $MIN_DEPTH_STR \
-a $MIN_AF_STR \
-r $MIN_PERC_STR

mv $STACKS_OUTPUT/batch_1.{fst*,hap*,pop*,sum*,vcf} $STACKS_OUTPUT/default/


#############################
# Load data into a database #
#############################

### EXPERIMENTAL, STILL NEED TO TEST AND DEBUG ###

# Load in data from this pipeline into a MySQL database 
#  -D  Database to load data into.
#  -p  Path to input files.
#  -b  Batch ID.
#  -M  if you have analyzed several populations, specify a population map.
#  -c  Load the catalog into the database.
#  -B  Load information into batch table.
#  -e  batch dEscription.
#  -d  perform a dry run. Do not actually load any data, just print what would be executed.
#  -W  only load file found on this white list.
#  -U  do not load stacks to unique_tags table to save database space.
#  -t  pipeline type (either 'map' or 'population'), load_radtags.pl will guess based on the presence/absence of progeny file types.

#echo "\nLoading all data into database: $DATABASE_NAME"
#$STACKS_PATH/load_radtags.pl -D $DATABASE_NAME \
#							 -p $STACKS_OUTPUT \
#							 -b $BATCH_ID \
#							 -M $PROJECT_ROOT/barcodes/population_map.tsv \
#							 -c -B -e $BATCH_NAME
#							 -d -t population

