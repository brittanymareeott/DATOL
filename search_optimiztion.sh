#!/usr/bin/env bash
#
# search_optimiztion.sh
#
# Author: Gregory Mendez
#
# This is the command script that takes the output from loop.sh, pretree_loop.sh, and raxml
# and builds new hmms based on that data to repeat search using optimized hmms and blast queries
# selected from the output of loop.sh.
#
# NOTE: Scripts called by this script use X11 to generate images. If using ssh to execute this on a remote 
# headless server use ssh -Y when connecting to the remote server to use your local X11 installation.
#
# Software Dependencies:
# RAxML 8.x
# MAFFT v7.215 http://mafft.cbrc.jp/alignment/software/
# HMMER 3.1 http://hmmer.org/
# Python 2.7
#   Bio
#   numpy
#   matplotlib
#   ete2
#
# Scripts called by this script:
# 1) dir_nexus_to_phylip.sh
# 2) nexus_to_phylip.py
# 4) prune_tree.py
# 5) long_branches.py
# 6) distance_matrix_zscore.py
# 7) rebuild_queries.sh
# 8) new_score_genes.sh
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
    -og|--out_groups)
    OUTGROUPS="$2"
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
echo Search optimization LOOP run on $(date) >> $INPUT/log.txt
echo Input Directory = "$INPUT" >> $INPUT/log.txt
echo Threads = "$THREADS" >> $INPUT/log.txt

# find working directory name
WORKING=$( find $INPUT -name 'Working_Dir_*' -exec basename {} \; )
echo Working Directory = "$WORKING" >> $INPUT/log.txt

# Generate Gene Trees
# make gene_trees directory
mkdir -p $INPUT/$WORKING/gene_trees/CDS
mkdir -p $INPUT/$WORKING/gene_trees/PEP
# Convert nexus files to phylip
printf "***************************   Converting Nexus files to Phylip ***************************\n"
echo "Starting dir_nexus_to_phylip.sh CDS run on $(date)" >> $INPUT/log.txt
echo "dir_nexus_to_phylip.sh -i $INPUT/$WORKING/trimal/CDS -o $INPUT/$WORKING/gene_trees/CDS -e nex -t $THREADS" >> $INPUT/log.txt
dir_nexus_to_phylip.sh -i $INPUT/$WORKING/trimal/CDS -o $INPUT/$WORKING/gene_trees/CDS -e nex -t $THREADS
echo "Starting dir_nexus_to_phylip.sh PEP run on $(date)" >> $INPUT/log.txt
echo "dir_nexus_to_phylip.sh -i $INPUT/$WORKING/trimal/PEP -o $INPUT/$WORKING/gene_trees/PEP -e nex -t $THREADS" >> $INPUT/log.txt
dir_nexus_to_phylip.sh -i $INPUT/$WORKING/trimal/PEP -o $INPUT/$WORKING/gene_trees/PEP -e nex -t $THREADS
# now make gene trees constrained by the best tree from the supermatrix tree
# prune the constraint tree so the species match the alignment file
printf "***************************   Creating custom pruned constraint tree for each gene ***************************\n"
echo "Starting custom constraint trees on $(date)" >> $INPUT/log.txt
cp $INPUT/$WORKING/supermatrix/CDS/RAxML_bestTree.loop_1_cds.tre $INPUT/$WORKING/gene_trees/CDS/constraint.tre
cd $INPUT/$WORKING/gene_trees/CDS/
echo "FILE=($(find $INPUT/$WORKING/gene_trees/CDS/*.phy -type f -exec basename {} \; | sed 's,.phy,,' ))" >> $INPUT/log.txt
FILE=($(find $INPUT/$WORKING/gene_trees/CDS/*.phy -type f -exec basename {} \; | sed 's,.phy,,' ))
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % prune_tree.py constraint.tre %".phy"
cp $INPUT/$WORKING/supermatrix/PEP/RAxML_bestTree.loop_1_pep.tre $INPUT/$WORKING/gene_trees/PEP/constraint.tre
cd $INPUT/$WORKING/gene_trees/PEP/
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % prune_tree.py constraint.tre %".phy"
# run raxml
# since raxml insists on at least 2 threads we need to set a variable that is 1/2 of our total threads rounding down
printf "***************************   Generating constrained gene trees ***************************\n"
echo "Starting RAxML constrained trees on $(date)" >> $INPUT/log.txt
RUNS=$( echo "scale=0;$THREADS/2" | bc -l )
cd $INPUT/$WORKING/gene_trees/CDS/
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f e -t %".constraint.tre" -m GTRGAMMA -s %".phy" -n %".constrained.tre"
cd $INPUT/$WORKING/gene_trees/PEP/
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f e -t %".constraint.tre" -m PROTGAMMAAUTO -s %".phy" -n %".constrained.tre"
# run raxml to make unconstrained gene trees
#raxml_genetrees.sh -i $INPUT/$WORKING/gene_trees/CDS -o $INPUT/$WORKING/gene_trees/CDS -e phy -t $THREADS -m GTRGAMMA
#raxml_genetrees.sh -i $INPUT/$WORKING/gene_trees/PEP -o $INPUT/$WORKING/gene_trees/PEP -e phy -t $THREADS -m PROTGAMMAAUTO
# run the long branch length analysis on the fully constrained trees
mkdir -p $INPUT/$WORKING/long_branches/CDS
mkdir -p $INPUT/$WORKING/long_branches/PEP
printf "***************************   Running long branch analysis on each gene tree ***************************\n"
echo "Starting long_branch.py analysis on $(date)" >> $INPUT/log.txt
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % long_branches.py --tree $INPUT/$WORKING/gene_trees/CDS/RAxML_result.%.constrained.tre --multi 7 --out_dir $INPUT/$WORKING/long_branches/CDS --outgroups $OUTGROUPS
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % long_branches.py --tree $INPUT/$WORKING/gene_trees/PEP/RAxML_result.%.constrained.tre --multi 7 --out_dir $INPUT/$WORKING/long_branches/PEP --outgroups $OUTGROUPS
# Generate raxml pairwise distance matrixes for each gene
printf "***************************   Generating distance matrixes for constrained trees ***************************\n"
echo "Starting RAxML distance matrix calculation on $(date)" >> $INPUT/log.txt
cd $INPUT/$WORKING/gene_trees/CDS
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f x -p 258755 -s %".phy" -n %".txt" -m GTRGAMMA -t %".constraint.tre"
cd $INPUT/$WORKING/gene_trees/PEP
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f x -p 258755 -s %".phy" -n %".txt" -m PROTGAMMAAUTO -t %".constraint.tre"
# Run analysis of distances to flag outliers for each sequence
printf "***************************   Running distance matrix analysis  ***************************\n"
echo "Starting distance_matrix_zscore.py on $(date)" >> $INPUT/log.txt
mkdir -p $INPUT/$WORKING/dist_m_zscore/CDS
mkdir -p $INPUT/$WORKING/dist_m_zscore/PEP
echo "distance_matrix_zscore.py --input $INPUT/$WORKING/gene_trees/CDS/ --score 30 --tree $INPUT/$WORKING/gene_trees/CDS/ --out $INPUT/$WORKING/dist_m_zscore/CDS/ --outgroups $OUTGROUPS" >> $INPUT/log.txt
distance_matrix_zscore.py --input $INPUT/$WORKING/gene_trees/CDS/ --score 30 --tree $INPUT/$WORKING/gene_trees/CDS/ --out $INPUT/$WORKING/dist_m_zscore/CDS/ --outgroups $OUTGROUPS
echo "distance_matrix_zscore.py --input $INPUT/$WORKING/gene_trees/PEP/ --score 30 --tree $INPUT/$WORKING/gene_trees/PEP/ --out $INPUT/$WORKING/dist_m_zscore/PEP/ --outgroups $OUTGROUPS" >> $INPUT/log.txt
distance_matrix_zscore.py --input $INPUT/$WORKING/gene_trees/PEP/ --score 30 --tree $INPUT/$WORKING/gene_trees/PEP/ --out $INPUT/$WORKING/dist_m_zscore/PEP/ --outgroups $OUTGROUPS
# Use the output from distance_matrix_z-score.py and long_branches.py to select species for new queries
printf "***************************   Building new query files ***************************\n"
echo "Starting rebuild_queries.sh on $(date)" >> $INPUT/log.txt
echo "rebuild_queries.sh -i $INPUT -t $THREADS" >> $INPUT/log.txt
rebuild_queries.sh -i $INPUT -t $THREADS
# Use new query sequences from rebuild_queries.sh to build hmms and generate bitscore cut off file
printf "***************************   Building new HMMs files and calculating bitscore cut offs   ***************************\n"
echo "Starting new_score_genes.sh on $(date)" >> $INPUT/log.txt
echo "new_score_genes.sh -i $INPUT/$WORKING/rebuilt_queries/query -o $INPUT/$WORKING/rebuilt_queries -e faa -t $THREADS" >> $INPUT/log.txt
new_score_genes.sh -i $INPUT/$WORKING/rebuilt_queries/query -o $INPUT/$WORKING/rebuilt_queries -e faa -t $THREADS
printf "***************************   DONE! Now run loop.sh again using the new queries.   ***************************\n"