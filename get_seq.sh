#!/usr/bin/env bash

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
  -t|--trans)
  TRANS="$2"
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
        OUT_FILE=$GENE'_'$SPECIES'.faa'
}
# We need to export the function and variables so they can be used in the subshell
export -f FIND_SPECIES_GENE
export OUT
export TRANS
cd $INPUT
# send all the text files into a bash subshell (run X subshells at a time) and run samtools faidx (a tools to pull fasta sequences
# out of a fasta file based on the def lines and do it right quick) on the files, feeding each line from the file to samfiles one at a time.
ls *.txt | xargs -n 1 -P 22 -I % bash -c 'FIND_SPECIES_GENE %; cat % | xargs samtools faidx $TRANS/$PEP_FILE > $OUT/$OUT_FILE;'
