#!/bin/bash
#
# Example of how to run the various stages of the Salmon RAD-seq pipeline. 

# Summary specs of bens machine: 
#     198 GB RAM
#     32 core Intel(R) Xeon(R) CPU E5-2470 0 @ 2.30GHz


# File names and path locations 
RDFS_PATH="/home/musselle/salmon"     # home/musselle is my home on bens machine and the rdf storage is mounted with the name 'salmon'
PROJECT_ROOT="$RDFS_PATH/chris"       # Root directory for project on RDFS Storage. 
BIN_PATH="~/bin"                     # Where all binaries for stacks are located


outpath="output"                                 # Root directory to write output to.  
raw_data="testdata/testset_100k.fastq"           # Path after 'root' for raw data files. Can accept a glob to pass multiple files 
barcodes="testdata/barcodes/lane3_newzebra.txt"  # Path after 'root' for barcodes associated with raw input files
cdhit_path="lib/bin/"                            # Location of binaries for cd-hit-est clustering program. cd-hit-v4.6.1 is included in repository
database="$outpath/test100k.db"                     # Output name for central database.
fasta_filepath="$outpath/testset100k_precluster.fasta"  # Location to write extracted fasta data for input to cd-hit clustering program. 
cluster_filepath="$outpath/testset100k_cluster"    # Location and file name to write CDHIT clustering output to. 


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
#
echo "\nAbout to Construct Unitag Reference Sequence"
# -i input file to use as unitag ref 
# -o output file 
# -m min read depth threshold
# -M max read depth threshold
$BIN_PATH/make_unitag.py -i $PROJECT_ROOT/processed-data/sample_GCAGGC.fq \
            -o $PROJECT_ROOT/processed-data/unitag/unitagref-m5-M500.fq \
            -m 5 -M 500     

# Construct index using bowtie
# bowtie-build reads_file index_file
bowtie-build $PROJECT_ROOT/processed-data/unitag/unitagref-m5-M500 \
		     $PROJECT_ROOT/processed-data/unitag/unitag_idx-m5-M500
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
						-x $PROJECT_ROOT/processed-data/unitag/unitag_idx-m5-M500 \
						-q 1 -m 2 \
						-b $PROJECT_ROOT/barcodes/barcodes_filenames.txt \
						-o $PROJECT_ROOT/stacks/pstacks/ -p 10 -k 1 

# Construct Catalogue from all pstacks output





