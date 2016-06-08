#!/usr/bin/env bash
#
# pretree_loop.sh
#
# Author: Gregg Mendez
#
# This script takes the output from loop.sh, user generated lists of genes and species, and creates alignment files
# for each gene (both dna and aa) and a supermatrix file of all the genes.
#
# Software Dependencies:
# trimal 1.2 http://trimal.cgenomics.org
# MAFFT v7.215 http://mafft.cbrc.jp/alignment/software/
# Python 2.7
#   Bio
#   glob
#   argparse
#
# You need to generate a couple input data files before running this script:
# 1) The output directory from loop.sh, this should include a Working Directory which has directories with data from each step of the loop.sh file.
# 2) A text file listing the species names you want to use in the analysis
# 3) A text file listing the genes you want to use in the analysis
#
# Scripts called by this script:
# subset_sorted.sh
# cat_mafft_all.sh
# supermatrix.py
# nexus_to_phylip.py
#
# Named variables. Every run needs the following defined:
# 1) -i | --input_dir - The directory containing the Working directory from loop.sh.
# 2) -s | --species_list - A text file listing species; one species per line.
# 3) -g | --gene_list - A text file listing genes; one gene per line.
# 4) -t | --threads - How many threads to use.
# 5) -p | --paralogs - [OPTIONAL] A directory containing files specifying sequences that have been identified as paralogous. The directory should contain text files for each gene titled in the format outlier_taxa.GENE.txt containing a species name on each line.
#
# Example:
# pretree_loop.sh -t 24 -i ~/bio/data/critter -s ~/bio/data/critter/Working_Dir_Mon_Dec_7_161512_EST_2015/lists/species_list.txt -g ~/bio/data/critter/Working_Dir_Mon_Dec_7_161512_EST_2015/lists/gene_list.txt
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
    -t|--threads)
    THREADS="$2"
    shift # past argument
    ;;
    -p|--paralogs)
    PARALOGS="$2"
    shift # past argument
    ;;
    *)
    # unknown option
    ;;
esac
shift # past argument or value
done

# Log variables and run time.
echo Pre Tree LOOP run on $(date) >> $INPUT/log.txt
echo Input Directory = "$INPUT" >> $INPUT/log.txt
echo Gene List = "$GENE_LIST" >> $INPUT/log.txt
echo Species List = "$SPECIES_LIST" >> $INPUT/log.txt
echo Threads = "$THREADS" >> $INPUT/log.txt

# find working directory name
LOOP_NUMBER=$(find $INPUT -path "$INPUT/loop*" -prune | wc -l )
LOOP_DIR=$INPUT/"loop_"$LOOP_NUMBER"_out/"
WORKING=$INPUT/"loop_"$LOOP_NUMBER"_out/tmp"
export WORKING
echo Working Directory = "$WORKING" >> $INPUT/log.txt

# Move desired sequences specified by gene and species lists to a new directory
mkdir $WORKING/subset_sorted
SUBSET_SORTED_OUT=$WORKING/subset_sorted
SUBSET_SORTED_INPUT=$LOOP_DIR/sequences/TopHits
if [ -z $PARALOGS  ]; then
    printf "***************************   Copying CDS files of genes and species specified ***************************\n"
    echo "Starting subset_sorted.sh CDS run on $(date)" >> $INPUT/log.txt
    echo "subset_sorted.sh -i $SUBSET_SORTED_INPUT/CDS -o $SUBSET_SORTED_OUT/CDS -s $SPECIES_LIST -g $GENE_LIST" >> $INPUT/log.txt
    subset_sorted.sh -i $SUBSET_SORTED_INPUT/CDS -o $SUBSET_SORTED_OUT/CDS -s $SPECIES_LIST -g $GENE_LIST
    echo "Completed subset_sorted.sh CDS run on $(date)" >> $INPUT/log.txt
    printf "***************************   Copying PEP files of genes and species specified ***************************\n"
    echo "Starting subset_sorted.sh PEP run on $(date)" >> $INPUT/log.txt
    echo "subset_sorted.sh -i $SUBSET_SORTED_INPUT/PEP -o $SUBSET_SORTED_OUT/PEP -s $SPECIES_LIST -g $GENE_LIST" >> $INPUT/log.txt
    subset_sorted.sh -i $SUBSET_SORTED_INPUT/PEP -o $SUBSET_SORTED_OUT/PEP -s $SPECIES_LIST -g $GENE_LIST
    echo "Completed subset_sorted.sh PEP run on $(date)" >> $INPUT/log.txt
else
    printf "***************************   Copying CDS files of genes and species specified ***************************\n"
    echo "Starting subset_sorted.sh CDS run on $(date)" >> $INPUT/log.txt
    echo "subset_sorted.sh -i $SUBSET_SORTED_INPUT/CDS -o $SUBSET_SORTED_OUT/CDS -s $SPECIES_LIST -g $GENE_LIST -p $PARALOGS" >> $INPUT/log.txt
    subset_sorted.sh -i $SUBSET_SORTED_INPUT/CDS -o $SUBSET_SORTED_OUT/CDS -s $SPECIES_LIST -g $GENE_LIST -p $PARALOGS
    echo "Completed subset_sorted.sh CDS run on $(date)" >> $INPUT/log.txt
    printf "***************************   Copying PEP files of genes and species specified ***************************\n"
    echo "Starting subset_sorted.sh PEP run on $(date)" >> $INPUT/log.txt
    echo "subset_sorted.sh -i $SUBSET_SORTED_INPUT/PEP -o $SUBSET_SORTED_OUT/PEP -s $SPECIES_LIST -g $GENE_LIST -p $PARALOGS" >> $INPUT/log.txt
    subset_sorted.sh -i $SUBSET_SORTED_INPUT/PEP -o $SUBSET_SORTED_OUT/PEP -s $SPECIES_LIST -g $GENE_LIST -p $PARALOGS
    echo "Completed subset_sorted.sh PEP run on $(date)" >> $INPUT/log.txt
fi

# Make alignments using MAFFT
mkdir $WORKING/cat_mafft_all
CAT_MAFFT_ALL_OUT=$WORKING/cat_mafft_all
printf "***************************   Creating Alignments of CDS files ***************************\n"
echo "Starting cat_mafft_all.sh CDS run on $(date)" >> $INPUT/log.txt
echo "cat_mafft_all.sh -i $SUBSET_SORTED_OUT/CDS -o $CAT_MAFFT_ALL_OUT/CDS -e .fas -t $THREADS" >> $INPUT/log.txt
cat_mafft_all.sh -i $SUBSET_SORTED_OUT/CDS -o $CAT_MAFFT_ALL_OUT/CDS -e .fas -t $THREADS
echo "Completed cat_muscle_all.sh CDS run on $(date)" >> $INPUT/log.txt
printf "***************************   Creating Alignments of PEP files ***************************\n"
echo "Starting cat_mafft_all.sh PEP run on $(date)" >> $INPUT/log.txt
echo "cat_mafft_all.sh -i $SUBSET_SORTED_OUT/PEP -o $CAT_MAFFT_ALL_OUT/PEP -e .fas -t $THREADS" >> $INPUT/log.txt
cat_mafft_all.sh -i $SUBSET_SORTED_OUT/PEP -o $CAT_MAFFT_ALL_OUT/PEP -e .fas -t $THREADS
echo "Completed cat_mafft_all.sh PEP run on $(date)" >> $INPUT/log.txt

# Remove sequence IDs from trimmed alignments prior to creating supermatrix
mkdir -p $WORKING/supermatrix/CDS/tmp
cd $CAT_MAFFT_ALL_OUT/CDS
find *.fas | xargs -n 1 -P $THREADS -I % bash -c 'sed "s,\(>.*\)___.*,\1,g" % > "$WORKING/supermatrix/CDS/tmp/"%'
# Now the PEP files
mkdir -p $WORKING/supermatrix/PEP/tmp
cd $CAT_MAFFT_ALL_OUT/PEP
find *.fas | xargs -n 1 -P $THREADS -I % bash -c 'sed "s,\(>.*\)___.*,\1,g" % > "$WORKING/supermatrix/PEP/tmp/"%'

# Trim alignments of poorly aligned sections using TrimAL
# mkdir $WORKING/trimal
# TRIMAL_OUT=$WORKING/trimal
# mkdir $TRIMAL_OUT/CDS
# mkdir $TRIMAL_OUT/PEP
# export TRIMAL_OUT
# cd $CAT_MAFFT_ALL_OUT/CDS
# printf "***********   Starting trimal for CDS sequences on `date` ...\n"
# echo "Starting trimal CDS run on $(date)" >> $INPUT/log.txt
# echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1" >> $INPUT/log.txt
# find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1
# echo "Completed trimal CDS run on $(date)" >> $INPUT/log.txt
# cd $TRIMAL_OUT/CDS/
# # trimal seems to have a bug where it lists all nexus output as protein data, so here we change it to dna
# echo "Correcting trimal bug that lists all nexus data as protein data on $(date)" >> $INPUT/log.txt
# echo "find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %" >> $INPUT/log.txt
# find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %

# Trim alignments of poorly aligned sections using TrimAL
printf "***********   Starting trimal for CDS sequences on `date` ...\n"
cd $WORKING/supermatrix/CDS/tmp/
echo "Starting trimal CDS run on $(date)" >> $INPUT/log.txt
echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out %'.nex' -automated1" >> $INPUT/log.txt
find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out %'.nex' -automated1
# trimal seems to have a bug where it lists all nexus output as protein data, so here we change it to dna
echo "Correcting trimal bug that lists all nexus data as protein data on $(date)" >> $INPUT/log.txt
echo "find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %" >> $INPUT/log.txt
find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %
printf "***********   Starting trimal for PEP sequences on `date` ...\n"
cd $WORKING/supermatrix/PEP/tmp/
echo "Starting trimal PEP run on $(date)" >> $INPUT/log.txt
echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out %'.nex' -automated1" >> $INPUT/log.txt
find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out %'.nex' -automated1


# Merge alignments into a single supermatrix alignment and save as Nexus and Phylip
printf "***********   Making supermatrix for CDS sequences on `date` ...\n"
echo "Starting supermatrix.py CDS run on $(date)" >> $INPUT/log.txt
echo "supermatrix.py --in_dir $TRIMAL_OUT/CDS/ --out $WORKING/supermatrix/CDS/supermatrix_cds.nex" >> $INPUT/log.txt
supermatrix.py --in_dir $WORKING/supermatrix/CDS/tmp/ --out $WORKING/supermatrix/CDS/supermatrix_cds.nex

# Remove file paths from nexus file gene names
cd $WORKING/supermatrix/CDS/
sed -i "s,$WORKING/supermatrix/CDS/tmp/,,g" supermatrix_cds.nex
# Convert nexus file to phylip
echo "Starting nexus_to_phylip.py CDS run on $(date)" >> $INPUT/log.txt
echo "nexus_to_phylip.py --input supermatrix_cds.nex --out supermatrix_cds.phylip" >> $INPUT/log.txt
nexus_to_phylip.py --input supermatrix_cds.nex --out supermatrix_cds.phylip

# Repeat Trimal and supermatrix steps now with the PEP files
# cd $CAT_MAFFT_ALL_OUT/PEP
# printf "***********   Starting trimal for PEP sequences on `date` ...\n"
# echo "Starting trimal PEP run on $(date)" >> $INPUT/log.txt
# echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/PEP/%'.nex' -automated1" >> $INPUT/log.txt
# find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/PEP/%'.nex' -automated1
# echo "Completed trimal PEP run on $(date)" >> $INPUT/log.txt
# cd $TRIMAL_OUT/PEP/

printf "***********   Making supermatrix for PEP sequences on `date` ...\n"
echo "Starting supermatrix.py PEP run on $(date)" >> $INPUT/log.txt
echo "supermatrix.py --in_dir $WORKING/supermatrix/PEP/tmp/ --out $WORKING/supermatrix/PEP/supermatrix_pep.nex" >> $INPUT/log.txt
supermatrix.py --in_dir $WORKING/supermatrix/PEP/tmp/ --out $WORKING/supermatrix/PEP/supermatrix_pep.nex

# remove file paths from nexus file gene names
cd $WORKING/supermatrix/PEP/
sed -i "s,$WORKING/supermatrix/PEP/tmp/,,g" supermatrix_pep.nex

# Convert nexus to phylip
printf "***********   Converting supermatrix PEP nexus file to phylip on `date` ...\n"
echo "Starting nexus_to_phylip.py PEP run on $(date)" >> $INPUT/log.txt
echo "nexus_to_phylip.py --input supermatrix_pep.nex --out supermatrix_pep.phylip" >> $INPUT/log.txt
nexus_to_phylip.py --input supermatrix_pep.nex --out supermatrix_pep.phylip
