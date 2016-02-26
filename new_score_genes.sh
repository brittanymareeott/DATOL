#!/usr/bin/env bash

# new_score_genes.sh

# Author: Gregg Mendez

# This script generates the necessary HMMs and HMM bitscore cut offs necessary to run the loop 
# with new genes. It takes fasta files as input. A separate fasta file for each gene is required.
# Multiple sequences for each gene are required.
#
# The output will be a directory containing HMMs (and alignments) for each gene and 2 files:
# averages.txt and score_cutoffs.txt 
#
# The input files MUST be named as follows:
# fasta filename = [gene name].[extension]
# example file names: KOG2606.fasta, BUSCO2606.fa, 000975.fsa
# fasta description lines = >[species name]___[gene name]
# example description lines: 
# >7294285___KOG2606
# >Amphidinium_carterae_MMET_c20244_g2_i1|m.35950___9658 c20244_g2_i1|g.35950  ORF c20244_g2_i1|g.35950 c20244_g2_i1|m.35950 type:5prime_partial len:164 (+) c20244_g2_i1:2-493(+) | aligned:3-164 (164)

# This script takes 4 named variables:
# 1) -i | --input_dir - The directory containing the fasta files to be scored.
# 2) -e | --extension - The file extension of the fasta input files. (Do not include the dot)
# 3) -t | --threads - The number of threads to use. Each thread will allow for another gene cutoff score to be calculated in parellel.
# 4) -o | --output_dir - The directory where you want the cutoff scores written.

#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
  -i|--input_dir)
  IN="$2"
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

# Function to find the gene and species names from given description lines
function FIND_SPECIES()
{
    SPLIT_FILE=($(echo $1 | tr " " "\n" | sed 's,___,\n,g' ))
    SPECIES=$( echo ${SPLIT_FILE[0]} | sed -e 's,>,,' )
}
# Function to break multi-fasta file into individual fasta files named as needed
function SPLIT_MULTI_FASTA() {
    INPUT=$1
    OUTPUT=$2/tmp/$GENE
    while read LINE
        do
            if [[ $LINE =~ ^\> ]]
                then
                    FIND_SPECIES $LINE
                    echo $LINE >> $OUTPUT/$SPECIES".fsa"
                else
                    if [[ $LINE == *\* ]]
                        then
                            NEW_LINE=${LINE/\*/__STOP_CODON__}
                            echo $NEW_LINE >> $OUTPUT/$SPECIES".fsa"
                        else
                            echo $LINE >> $OUTPUT/$SPECIES".fsa"
                    fi
            fi
            sed -i 's,__STOP_CODON__,\*,g' $OUTPUT/$SPECIES".fsa"
        done < $INPUT
}

#################################
# WORK
#################################

# make a temporary directory
mkdir -p $OUT/tmp
mkdir -p $OUT/hmms
# move into the input directory
cd $IN

# Inside each subshell we are first breaking each multi-fasta file into individual 
# fasta files for each sequence, then creating a list of all those single sequence
# files. We then loop through that list and create an alignment for every combination
# of sequences missing one of the sequences. An hmm is built for each of those alignments
# and hmmsearch is run with the missing sequences vs the hmm that doesn't contain it.
# The bitscores from all of these searches is averaged and then divided by 2 to generate
# a cut off score for that gene. The averages are saved to averages.txt and the cutoff scores
# are saved to score_cutoffs.txt.
function CALC_SCORES() {
    GENE=${FILE/.$EXT/}
    mkdir $OUT/tmp/$GENE
    SPLIT_MULTI_FASTA $FILE $OUT
    cd $OUT/tmp/$GENE
    ALIGNMENT=(*.fsa)
    for j in ${ALIGNMENT[@]}
        do ARRAY_BOB=(${ALIGNMENT[@]/$j})
            cat ${ARRAY_BOB[@]} > allminus_${j/fsa/fasta}
            mafft --auto allminus_${j/fsa/fasta} > allminus_${j/fsa/fa}
            hmmbuild allminus_${j/fsa/hmm} allminus_${j/fsa/fa}
            hmmsearch allminus_${j/fsa/hmm} $j > ${j/fsa/out}
        done
    MATHS=0
    for h in *.out
        do SCORE=($(sed -n "15 p" $h))
            NUMBER=${SCORE[1]}
            MATHS=$(echo $NUMBER + $MATHS | bc -l)
            ITEMS=$(ls *.out | wc -l)
        done
    AVERAGE=$(echo "scale=3; $MATHS/$ITEMS" | bc -l)
    CUTOFF=$(echo "scale=3; $AVERAGE/2" | bc -l)
    echo $GENE $AVERAGE >> $OUT/average.txt
    echo $GENE $CUTOFF >> $OUT/score_cutoffs.txt
}
# make variables and funtions available to our subshells
export EXT
export OUT
export IN
export -f SPLIT_MULTI_FASTA
export -f FIND_SPECIES
export -f CALC_SCORES

# Feed everything ending in specified extension into its own sub shell to calculate scores.
# Then create alignments and build HMMs from the alignments.
find *.$EXT | xargs -n 1 -P $THREADS -I % bash -c 'FILE=% ; \
    CALC_SCORES ; \
    cd $IN ; \
    mafft --auto $FILE > $OUT/hmms/$GENE.fasta ; \
    hmmbuild $OUT/hmms/$GENE.hmm $OUT/hmms/$GENE.fasta'
# delete the tmp directory
rm -r $OUT/tmp

