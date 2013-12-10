#!/bin/bash
#
# Example of how to run the various stages of the Salmon RAD-seq pipeline. 
#
########################################################################
# PART II: Catalogue Construction, Matching and Populations Statistics #
########################################################################

# Summary specs of bens machine: 
#     198 GB RAM
#     32 core Intel(R) Xeon(R) CPU E5-2470 0 @ 2.30GHz


######################
# PARAMETERS DEFINED #
######################

# Population Parameters 
#  MIN_PERC     minimum percentage of individuals in a population required to process a locus for that population. 
#  MIN_DEPTH    specify a minimum stack depth required for individuals at a locus.
#  MIN_AF       specify a minimum minor allele frequency required before calculating Fst at a locus (0 < a < 0.5).
MIN_DEPTH=10
MIN_AF=0.1
MIN_PERC=80

# Aligning stacks to reference and pstacks
CATALOGUE_MIN_DIST=0

# Unique Batch Number for this run 
BATCH_ID=2
BATCH_NAME="unitag"

# Num processors to run on.
NUM_THREADS=15

# File names and path locations 
RDFS_PATH="/home/musselle/salmon"     # home/musselle is my home on bens machine and the rdf storage is mounted with the name 'salmon'
PROJECT_ROOT="$RDFS_PATH/chris"       # Root directory for project on RDFS Storage. 
BIN_PATH="$HOME/bin"                      # Where all python runscripts are located
STACKS_PATH="$HOME/bin/stacks/bin"         # Where all binaries for stacks are located. 
STACKS_OUTPUT="$PROJECT_ROOT/stacks/$BATCH_NAME"   # Where all output stacks files are written to by default 

DATABASE_NAME="$PROJECT_ROOT/salmon_database.db"

mkdir -p $STACKS_OUTPUT

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
$BIN_PATH/run_cstacks.py -i $STACKS_OUTPUT/sample_* \ 
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

echo "\nAbout to Calculate Population Statistics" 
$STACKS_PATH/populations -b $BATCH_ID \
                         -P $STACKS_OUTPUT \ 
                         -M $PROJECT_ROOT/barcodes/population_map.tsv \
                         --vcf \
                         -t $NUM_THREADS \
                         -m $MIN_DEPTH \
                         -a $MIN_AF \
                         -r $MIN_PERC

#############################
# Load data into a database #
#############################

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

