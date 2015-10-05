#!/usr/bin/env bash

# This script takes lists of species and genes and loops through a directory
# like that produced by cleanup_explode_sort.sh and and copies just those
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

# loop through input directory and copy genes and species based on whether they are in the species and gene lists
mkdir $OUT 2>/dev/null
cd $INPUT
for i in *
  do
    if [ -d $i ]
     then 
       cd $i
       if [ $( ARRAY_CONTAINS GENES $i && echo yes || echo no ) == yes  ]
       then
        printf "********************* COPYING $i\n"
        GENE=${PWD##*/}
        mkdir -p $OUT/$GENE
        for g in *.fas
        do
          ARRAY_CONTAINS SPECIES ${g/.fas/} && cp $g $OUT/$GENE && echo ${g/.fas/} copied || echo ${g/.fas/} not copied
        done
       else
        printf "********************* SKIPPING  $i\n"
       fi
     cd $INPUT
    fi
  done
