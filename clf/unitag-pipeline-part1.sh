#!/bin/bash
#
# Example of how to run the various stages of the Salmon RAD-seq pipeline. 
#
################################################
# PART I: Preprocessing and Forming the Stacks #
################################################

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