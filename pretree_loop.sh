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
#   BioPython
#   glob
#   argparse
#   AlignIO
#   Nexus
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
WORKING=$( find $INPUT -name 'Working_Dir_*' -exec basename {} \; )
echo Working Directory = "$WORKING" >> $INPUT/log.txt

# Move desired sequences specified by gene and species lists to a new directory
mkdir $INPUT/$WORKING/subset_sorted
SUBSET_SORTED_OUT=$INPUT/$WORKING/subset_sorted
SUBSET_SORTED_INPUT=$INPUT/$WORKING/write_cds_pep/TopHits
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

# Make alignments using MAFFT
mkdir $INPUT/$WORKING/cat_mafft_all
CAT_MAFFT_ALL_OUT=$INPUT/$WORKING/cat_mafft_all
printf "***************************   Creating Alignments of CDS files ***************************\n"
echo "Starting cat_mafft_all.sh CDS run on $(date)" >> $INPUT/log.txt
echo "cat_mafft_all.sh -i $SUBSET_SORTED_OUT/CDS -o $CAT_MAFFT_ALL_OUT/CDS -e .fas -t $THREADS" >> $INPUT/log.txt
cat_mafft_all.sh -i $SUBSET_SORTED_OUT/CDS -o $CAT_MAFFT_ALL_OUT/CDS -e .fas -t $THREADS
echo "Completed cat_mafft_all.sh CDS run on $(date)" >> $INPUT/log.txt
printf "***************************   Creating Alignments of PEP files ***************************\n"
echo "Starting cat_mafft_all.sh PEP run on $(date)" >> $INPUT/log.txt
echo "cat_mafft_all.sh -i $SUBSET_SORTED_OUT/PEP -o $CAT_MAFFT_ALL_OUT/PEP -e .fas -t $THREADS" >> $INPUT/log.txt
cat_mafft_all.sh -i $SUBSET_SORTED_OUT/PEP -o $CAT_MAFFT_ALL_OUT/PEP -e .fas -t $THREADS
echo "Completed cat_mafft_all.sh PEP run on $(date)" >> $INPUT/log.txt

# Trim alignments of poorly aligned sections using TrimAL
mkdir $INPUT/$WORKING/trimal
TRIMAL_OUT=$INPUT/$WORKING/trimal
mkdir $TRIMAL_OUT/CDS
mkdir $TRIMAL_OUT/PEP
export TRIMAL_OUT
cd $CAT_MAFFT_ALL_OUT/CDS
printf "***********   Starting trimal for CDS sequences on `date` ...\n"
echo "Starting trimal CDS run on $(date)" >> $INPUT/log.txt
echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1" >> $INPUT/log.txt
find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1
echo "Completed trimal CDS run on $(date)" >> $INPUT/log.txt
cd $TRIMAL_OUT/CDS/
# trimal seems to have a bug where it lists all nexus output as protein data, so here we change it to dna
echo "Correcting trial bug that lists all nexus data as protein data on $(date)" >> $INPUT/log.txt
echo "find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %" >> $INPUT/log.txt
find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %
# Merge alignments into a single supermatrix alignment and save as Nexus and Phylip
printf "***********   Making supermatrix for CDS sequences on `date` ...\n"
echo "Starting supermatrix.py CDS run on $(date)" >> $INPUT/log.txt
mkdir -p $INPUT/$WORKING/supermatrix/CDS
echo "supermatrix.py --in_dir $TRIMAL_OUT/CDS/ --out $INPUT/$WORKING/supermatrix/CDS/supermatrix_cds.nex" >> $INPUT/log.txt
supermatrix.py --in_dir $TRIMAL_OUT/CDS/ --out $INPUT/$WORKING/supermatrix/CDS/supermatrix_cds.nex
echo "Completed supermatrix.py CDS run on $(date)" >> $INPUT/log.txt
printf "***********   Converting supermatrix CDS nexus file to phylip on `date` ...\n"
cd $INPUT/$WORKING/supermatrix/CDS/
sed -i "s,$TRIMAL_OUT/CDS/,,g" supermatrix_cds.nex
echo "Starting nexus_to_phylip.py CDS run on $(date)" >> $INPUT/log.txt
echo "nexus_to_phylip.py --input supermatrix_cds.nex --out supermatrix_cds.phylip" >> $INPUT/log.txt
nexus_to_phylip.py --input supermatrix_cds.nex --out supermatrix_cds.phylip
echo "Completed nexus_to_phylip.py CDS run on $(date)" >> $INPUT/log.txt
 
# Repeat Trimal and supermatrix steps now with the PEP files
cd $CAT_MAFFT_ALL_OUT/PEP
printf "***********   Starting trimal for PEP sequences on `date` ...\n"
echo "Starting trimal PEP run on $(date)" >> $INPUT/log.txt
echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/PEP/%'.nex' -automated1" >> $INPUT/log.txt
find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/PEP/%'.nex' -automated1
echo "Completed trimal PEP run on $(date)" >> $INPUT/log.txt
cd $TRIMAL_OUT/PEP/
mkdir -p $INPUT/$WORKING/supermatrix/PEP
printf "***********   Making supermatrix for PEP sequences on `date` ...\n"
echo "Starting supermatrix.py PEP run on $(date)" >> $INPUT/log.txt
echo "supermatrix.py --in_dir $TRIMAL_OUT/PEP/ --out $INPUT/$WORKING/supermatrix/PEP/supermatrix_pep.nex" >> $INPUT/log.txt
supermatrix.py --in_dir $TRIMAL_OUT/PEP/ --out $INPUT/$WORKING/supermatrix/PEP/supermatrix_pep.nex
echo "Completed supermatrix.py PEP run on $(date)" >> $INPUT/log.txt
printf "***********   Converting supermatrix PEP nexus file to phylip on `date` ...\n"
cd $INPUT/$WORKING/supermatrix/PEP/
sed -i "s,$TRIMAL_OUT/PEP/,,g" supermatrix_pep.nex
echo "Starting nexus_to_phylip.py PEP run on $(date)" >> $INPUT/log.txt
echo "nexus_to_phylip.py --input supermatrix_pep.nex --out supermatrix_pep.phylip" >> $INPUT/log.txt
nexus_to_phylip.py --input supermatrix_pep.nex --out supermatrix_pep.phylip
echo "Completed nexus_to_phylip.py PEP run on $(date)" >> $INPUT/log.txt