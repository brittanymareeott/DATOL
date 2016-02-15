#!/usr/bin/env bash
#
# cat_mafft_all.sh
#
# Author: Gregg Mendez
#
# This script creates alignments, using MAFFT, for separate fasta files in a directory in a highly parellel manner.
#
# Named variables. Every run needs the following defined:
# 1) -i | --input_dir - The directory containing the fasta files that need to be aligned.
# 2) -e | --file_extension - The file extension of the fasta input files.
# 3) -t | --threads - The number of threads to use for each mafft run. NOTE: -t X -r MUST NOT exceed the total number of threads available.
# 4) -o | --output_dir - The directory where you want the alignments saved.
#
# Example: cat_mafft_all.sh -i /home/mendezg/fasta -o /home/mendezg/alignments -e .fas -t 4

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
  -e|--extension)
  EXT="$2"
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

# function to round numbers
round() {
    printf "%.$2f" "$1"
}
# function for doing math
math() {
    echo "$*" | bc -l
}
# function to format seconds into hours, minutes, seconds
convertsecs() {
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "%02dh %02dm %02ds\n" $h $m $s
}


# make output directory
mkdir -p $OUT 2>/dev/null

# Set some variables to be used to keep track of progress during mafft runs
GENES=($(find $INPUT -type d -exec basename {} \;))
TOTAL=${#GENES[@]:1}
t0=$( date +%s )
# We need to export the function and variables so they can be used in the subshell
export OUT
export INPUT
export EXT
export TOTAL
export t0
export -f round
export -f math
export -f convertsecs

# find all the folders and launch a separate bash shell for each. In each shell concatenate the sequences then align them
printf "%s\n" "${GENES[@]:1}" | xargs -n 1 -P $THREADS -I %x bash -c 'cd $INPUT/%x; \
    cat *$EXT > $OUT/%x_combined.fasta; \
    mafft --quiet --thread 1 --auto $OUT/%x_combined.fasta > $OUT/%x.fas; \
    COUNT=$( find $OUT/*.fas -type f -size +1c -exec basename {} \; | wc -l ); \
    PROGRESS=$( math "$COUNT / $TOTAL *100" ); \
    tx=$( date +%s ); \
    DURATION=$( math "$tx - $t0" ); \
    REMAINING=$( math "(100-$PROGRESS)*($DURATION/$PROGRESS)" ); \
    printf "*************************** %s%% Completed ********* Estimated Time Remaining %s\r" "$(round $PROGRESS)" "$(convertsecs $(round $REMAINING))"'
printf "\n"