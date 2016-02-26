#!/usr/bin/env bash

# Author Gregory S Mendez
#
# This script quickly fetches sequences based on def-lines from a text file and writes a new fasta file.
# It uses samtools faidx http://www.htslib.org to do the real work. This script just uses the text file
# file names to set the input fasta file and output file name, and parallelizes the process using xargs.
# This script differs from the similar get_seq.sh script in that it outputs both peptide and cds fasta files,
# as well as concatenating all the sequences down to one file per species.
#
# The script takes 4 arguments:
# 1) -i | --input_dir - The directory containing the text files with the def-lines of the sequences to fetch.
# 2) -o | --output_dir - The directory where the new fasta files should be written.
# 3) -f | --fasta - The directory containing the large fasta files to fetch the sequences from.
# 4) -t | --threads - How many threads to use.
#
# Usage: get_seq -i ~/GreenAlgae/hmmsearch_txt -f ~/GreenAlgae/translations -o ~/GreenAlgae/hmmsearch_out -t 24
#
# The input files should be names as follows:
# genename_taxon_id_with_underscores.txt
# The first part is the gene name/identifier. This must not contain underscores or periods.
# The second part is the species name/identifier. This can contain underscores or period.
# There must be a file extension at the end of the file name.
# 
# The contents of the file should contain a fasta def line on each line. The script parse_hmm_search.py generates
# the input text file expected by this script.

#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
  -i|--input_dir)
  INPUT="$2"
  shift # past argument
  ;;
  -o|--output_dir)
  OUT="$2"
  shift # past argument
  ;;
  -f|--fasta)
  FASTA="$2"
  shift # past argument
  ;;
  -t|--threads)
  THREADS="$2"
  shift # past argument
  ;;
  *)
        # unknown option
  ;;
esac
shift # past argument or value
done

# Function to parse the input file name and set variables for gene name, species, translation locations, and output file
function FIND_SPECIES_GENE()
{
        SPLIT_FILE=($(echo $1 | tr "\." "\n" | tr "_" "\n"))
        unset SPLIT_FILE[${#SPLIT_FILE[@]}-1]
        GENE=${SPLIT_FILE[0]}
        SPECIES_SPACED=${SPLIT_FILE[@]:1}
        SPECIES=${SPECIES_SPACED// /_}
        PEP_FILE=$SPECIES'.fa.transdecoder.pep'
        CDS_FILE=$SPECIES'.deg.fas'
        OUT_FILE=$GENE'_'$SPECIES'.faa'
        
}
# We need to export the function and variables so they can be used in the subshell
export -f FIND_SPECIES_GENE
export OUT
export FASTA

# We need to create our output directories
mkdir -p $OUT/TopHits/CDS $OUT/TopHits/PEP $OUT/OtherHits/PEP/tmp $OUT/OtherHits/CDS/tmp 2> /dev/null

# First lets do the Top Hits
cd $INPUT/TopHits
# send all the text files into a bash subshell (run X subshells at a time) and run samtools faidx (a tools to pull fasta sequences
# out of a fasta file based on the def lines and do it right quick) on the files, feeding each line from the file to samfiles one at a time.
find *.txt | xargs -n 1 -P $THREADS -I % bash -c 'FIND_SPECIES_GENE %; \
cat % | xargs samtools faidx $FASTA/$PEP_FILE > $OUT/TopHits/PEP/$OUT_FILE; \
cat % | xargs samtools faidx $FASTA/$CDS_FILE > $OUT/TopHits/CDS/$OUT_FILE;'
# Sort the output based on gene and rename files and def lines to just the species name
# First the peptide files
cd $OUT/TopHits/PEP
find *.faa | xargs -n 1 -P $THREADS -I % bash -c 'FIND_SPECIES_GENE %; \
sed -i "s,>.*,>$SPECIES,g" %; \
mkdir -p $OUT/TopHits/PEP/$GENE; \
mv % $OUT/TopHits/PEP/$GENE/$SPECIES".fas";'
# Next the cds files
cd $OUT/TopHits/CDS
find *.faa | xargs -n 1 -P $THREADS -I % bash -c 'FIND_SPECIES_GENE %; \
sed -i "s,>.*,>$SPECIES,g" %; \
mkdir -p $OUT/TopHits/CDS/$GENE; \
mv % $OUT/TopHits/CDS/$GENE/$SPECIES".fas";'

# Next lets do the Non-Top Hits
cd $INPUT/OtherHits
# send all the text files into a bash subshell (run X subshells at a time) and run samtools faidx (a tools to pull fasta sequences
# out of a fasta file based on the def lines and do it right quick) on the files, feeding each line from the file to samfiles one at a time.
find *.txt | xargs -n 1 -P $THREADS -I % bash -c 'FIND_SPECIES_GENE %; \
cat % | xargs samtools faidx $FASTA/$PEP_FILE > $OUT/OtherHits/PEP/tmp/$OUT_FILE; \
cat % | xargs samtools faidx $FASTA/$CDS_FILE > $OUT/OtherHits/CDS/tmp/$OUT_FILE;'
# Concatenate all the sequences from each species and append the gene name to the start of each def line
# First the peptide files
cd $OUT/OtherHits/PEP/tmp
find *.faa | xargs -n 1 -P $THREADS -I % bash -c 'FIND_SPECIES_GENE %; \
sed -i "s,>,>$SPECIES\___,g" %; \
cat % >> $OUT/OtherHits/PEP/$GENE".fasta";'
# Next the cds files
cd $OUT/OtherHits/CDS/tmp
find *.faa | xargs -n 1 -P $THREADS -I % bash -c 'FIND_SPECIES_GENE %; \
sed -i "s,>,>$SPECIES\___,g" %; \
cat % >> $OUT/OtherHits/CDS/$GENE".fasta";'