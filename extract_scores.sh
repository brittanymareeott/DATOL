#!/usr/bin/env bash

# This script loops through hmmsearch output files, saves the top bitscore to one file and saves all bitscores to another.
# Next it generates a histogram for each file using the python script histogram.py
# Output files are named for the query term of the hmmsearch.
# Output Files:
# 1) QUERY.txt = The list of all top hit bitscores
# 2) QUERY.all.txt = The list of all bitscores
# 3) QUERY.txt.pdf = The histogram of the top hits scores
# 4) QUERY.all.txt.pdf = The histogram of all the hit scores

# this script requires the script histogram.py written by me and packaged with this script.

# Example usage: extract_scores.sh -i /home/mendezg/hmmsearch_out -o scores

#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
  -i|--input)
  INPUT="$2"
  shift # past argument
  ;;
  -o|--output_directory)
  OUT="$2"
  shift # past argument
  ;;
  *)
        # unknown option
  ;;
esac
shift # past argument or value
done

cd $INPUT
mkdir $OUT
# Loop through all the files in the input directory ending in .out
for FILENAME in *.out
  do
# Loop through the lines of the file and save score of the top hit to a file
    while read LINE
      do
        if [[ ! $LINE =~ ^# ]]
          then
            ARR=($LINE)
           # echo ${ARR[2]} ${ARR[5]}
            echo ${ARR[5]} >> $OUT/${ARR[2]}".txt"
            break
        fi
      done < $FILENAME
# Loop through the lines of the file and save all the scores to a file
    while read LINE
      do
        if [[ ! $LINE =~ ^# ]]
          then
            ARR=($LINE)
           # echo ${ARR[2]} ${ARR[5]}
            echo ${ARR[5]} >> $OUT/${ARR[2]}".all.txt"
        fi
      done < $FILENAME
  done

# Loop through the score files and make histograms for them all
cd $OUT
for FILENAME in *.txt
  do
    histogram.py $FILENAME
  done
