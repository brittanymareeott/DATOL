#!/usr/bin/env bash
#
# subset_sorted.sh
#
# Author: Gregory Mendez
#
# This script takes lists of species and genes and loops through a directory
# like that produced by write_cds_pep.sh and and copies just those
# species and genes in the lists to a new directory with the same structure.

# This script takes four arguments:
# 1) -i | --input_dir - The input directory as discussed above.
# 2) -o | --output_dir - The directory in which to write the the new genes directories and fasta files.
# 3) -s | --species_list - A text file containing a species name on each line.
# 4) -g | --gene_list - A text file containing a gene on each line.
#
# Example Usage: subset_sorted.sh -i /home/mendezg/crawly_things/all_sorted -o /home/mendezg/crawly_things/subset_sorted -s crawliest_species.txt -g best_genes.txt
#

#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
  -i|--input_dir)
  INPUT="$2"
  shift # past argument
  ;;
  -s|--species_list)
  SPECIES_LIST="$2"
  shift # past argument
  ;;
  -g|--gene_list)
  GENE_LIST="$2"
  shift # past argument
  ;;
  -o|--out_dir)
  OUT="$2"
  shift # past argument
  ;;
  *)
        # unknown option
  ;;
esac
shift # past argument or value
done

#LOAD input files into arrays
SPECIES=($(cat $SPECIES_LIST))
GENES=($(cat $GENE_LIST))

# function to check contents of arrays
ARRAY_CONTAINS () {
  local ARRAY="$1[@]"
  local SEEKING=$2
  local IN=1
  for ELEMENT in "${!ARRAY}"; do
    if [[ $ELEMENT == $SEEKING ]]; then
      IN=0
      break
    fi
  done
  return $IN
}
# function to round numbers
round() {
    printf "%.$2f" "$1"
}
# # function for doing math
math() {
    echo "$*" | bc -l
}
# # function to format seconds into hours, minutes, seconds
convertsecs() {
    ((h=${1}/3600))
    ((m=(${1}%3600)/60))
    ((s=${1}%60))
    printf "%02dh %02dm %02ds\n" $h $m $s
}
# loop through input directory and copy genes and species based on whether they are in the species and gene lists
mkdir $OUT 2>/dev/null
cd $INPUT
TOTAL=${#GENES[@]}
COUNT=0
SECONDS=0
for i in *
    do
        if [ -d $i ]
            then 
                cd $i
                if [ $( ARRAY_CONTAINS GENES $i && echo yes || echo no ) == yes  ]
                    then
                        GENE=${PWD##*/}
                        mkdir -p $OUT/$GENE
                        for g in *.fas
                            do
                                ARRAY_CONTAINS SPECIES ${g/.fas/} && cp $g $OUT/$GENE #&& echo ${g/.fas/} copied || echo ${g/.fas/} not copied                           
                            done
                COUNT=$( math $COUNT + 1 )
                PROGRESS=$( math "$COUNT / $TOTAL *100" )
                REMAINING=$( math "(100-$PROGRESS)*($SECONDS/$PROGRESS)" )
                printf '*************************** %s%% Completed ********* Estimated Time Remaining %s\r' "$(round $PROGRESS)" "$(convertsecs $(round $REMAINING))"
                fi
                cd $INPUT
        fi
    done
printf '\n'