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
# trimal 1.2 http://trimal.cgenomics.org
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
LOOP_NUMBER=$(find $INPUT -path "$INPUT/loop*" -prune | wc -l )
LOOP_DIR=$INPUT/"loop_"$LOOP_NUMBER"_out/"
WORKING=$INPUT/"loop_"$LOOP_NUMBER"_out/tmp"
REPORT=$INPUT/report
export LOOP_NUMBER
export INPUT
export WORKING
export REPORT
echo Working Directory = "$WORKING" >> $INPUT/log.txt

######################################
##      Generate Gene Trees
######################################

# make gene_trees directory
mkdir -p $WORKING/gene_trees/CDS
mkdir -p $WORKING/gene_trees/PEP

# Trim alignments of poorly aligned sections using TrimAL
mkdir $WORKING/gene_trees/trimal
TRIMAL_OUT=$WORKING/gene_trees/trimal
CAT_MAFFT_ALL_OUT=$WORKING/cat_mafft_all
mkdir $TRIMAL_OUT/CDS
mkdir $TRIMAL_OUT/PEP
export TRIMAL_OUT
cd $CAT_MAFFT_ALL_OUT/CDS
printf "***********   Starting trimal for CDS sequences on `date` ...\n"
echo "Starting trimal CDS run on $(date)" >> $INPUT/log.txt
echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1" >> $INPUT/log.txt
find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1
cd $TRIMAL_OUT/CDS/

# trimal seems to have a bug where it lists all nexus output as protein data, so here we change it to dna
echo "Correcting trimal bug that lists all nexus data as protein data on $(date)" >> $INPUT/log.txt
echo "find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %" >> $INPUT/log.txt
find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %

# Repeat Trimal steps now with the PEP files
cd $CAT_MAFFT_ALL_OUT/PEP
printf "***********   Starting trimal for PEP sequences on `date` ...\n"
echo "Starting trimal PEP run on $(date)" >> $INPUT/log.txt
echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/PEP/%'.nex' -automated1" >> $INPUT/log.txt
find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/PEP/%'.nex' -automated1
cd $TRIMAL_OUT/PEP/

# Convert nexus files to phylip
printf "***************************   Converting Nexus files to Phylip ***************************\n"
echo "Starting dir_nexus_to_phylip.sh CDS run on $(date)" >> $INPUT/log.txt
echo "dir_nexus_to_phylip.sh -i $TRIMAL_OUT/CDS -o $WORKING/gene_trees/CDS -e nex -t $THREADS" >> $INPUT/log.txt
dir_nexus_to_phylip.sh -i $TRIMAL_OUT/CDS -o $WORKING/gene_trees/CDS -e nex -t $THREADS
echo "Starting dir_nexus_to_phylip.sh PEP run on $(date)" >> $INPUT/log.txt
echo "dir_nexus_to_phylip.sh -i $TRIMAL_OUT/PEP -o $WORKING/gene_trees/PEP -e nex -t $THREADS" >> $INPUT/log.txt
dir_nexus_to_phylip.sh -i $TRIMAL_OUT/PEP -o $WORKING/gene_trees/PEP -e nex -t $THREADS

# now make gene trees constrained by the best tree from the supermatrix tree
# prune the constraint tree so the species match the alignment file
printf "***************************   Creating custom pruned constraint tree for each gene ***************************\n"
echo "Starting custom constraint trees on $(date)" >> $INPUT/log.txt
cp $WORKING/supermatrix/CDS/RAxML_bestTree.loop_"$LOOP_NUMBER"_cds.tre $WORKING/gene_trees/CDS/constraint.tre
cd $WORKING/gene_trees/CDS/
echo "FILE=($(find $WORKING/gene_trees/CDS/*.phy -type f -exec basename {} \; | sed 's,.phy,,' ))" >> $INPUT/log.txt
FILE=($(find $WORKING/gene_trees/CDS/*.phy -type f | sed 's#.*/##' | sed 's,.phy,,' ))
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % prune_tree.py constraint.tre %".phy"
cp $WORKING/supermatrix/PEP/RAxML_bestTree.loop_"$LOOP_NUMBER"_pep.tre $WORKING/gene_trees/PEP/constraint.tre
cd $WORKING/gene_trees/PEP/
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % prune_tree.py constraint.tre %".phy"

# run raxml
# since raxml insists on at least 2 threads we need to set a variable that is 1/2 of our total threads rounding down
printf "***************************   Generating constrained gene trees ***************************\n"
echo "Starting RAxML constrained trees on $(date)" >> $INPUT/log.txt
RUNS=$( echo "scale=0;$THREADS/2" | bc -l )
cd $WORKING/gene_trees/CDS/
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f e -t %".constraint.tre" -m GTRGAMMA -s %".phy" -n %".constrained.tre"
cd $WORKING/gene_trees/PEP/
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f e -t %".constraint.tre" -m PROTGAMMAAUTO -s %".phy" -n %".constrained.tre"

#####################################
#      Analyze Trees
#####################################

# Long Branch Analysis

# run the long branch length analysis on the fully constrained trees
mkdir -p $WORKING/long_branches/CDS
mkdir -p $WORKING/long_branches/PEP
printf "***************************   Running long branch analysis on each gene tree ***************************\n"
echo "Starting long_branch.py analysis on $(date)" >> $INPUT/log.txt
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % long_branches.py --tree $WORKING/gene_trees/CDS/RAxML_result.%.constrained.tre --multi 7 --out_dir $WORKING/long_branches/CDS --outgroups $OUTGROUPS
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $THREADS -I % long_branches.py --tree $WORKING/gene_trees/PEP/RAxML_result.%.constrained.tre --multi 7 --out_dir $WORKING/long_branches/PEP --outgroups $OUTGROUPS

# Generate Report for Long Branch Analysis
RFILE=($(find $REPORT/genes/*.html -type f | sed 's#.*/##'))
function REPORT_LONG {
    HTML_FILE=$1
    FILE_PATH=$INPUT/report/genes/$HTML_FILE
    GENE=${HTML_FILE/.html/}
    LONG_CDS=$WORKING/long_branches/CDS/"longbranch_taxa."$GENE".txt"
    LONG_PEP=$WORKING/long_branches/PEP/"longbranch_taxa."$GENE".txt"
    if [ -f $LONG_CDS ]
        then
            readarray SPECIES_LIST < $LONG_CDS
            for SPECIES in ${SPECIES_LIST[@]}
                do
                    ROW=\<\!--ROW_"$SPECIES"--\>
                    sed -i "s,$ROW,<td class='fail long_cds'>FAIL</td>$ROW," $FILE_PATH
                done
            sed -i "s,\(<td class='hmm'>[0-9]*</td>\)<\!--,\1<td class='long_cds'>PASS</td><\!--,g" $FILE_PATH
        else
            sed -i "s,<\!--ROW,<td class='long_cds'> </td><\!--ROW,g" $FILE_PATH
    fi
    if [ -f $LONG_PEP ]
        then
            readarray SPECIES_LIST < $LONG_PEP
            for SPECIES in ${SPECIES_LIST[@]}
                do
                    ROW=\<\!--ROW_"$SPECIES"--\>
                    sed -i "s,$ROW,<td class='fail long_pep'>FAIL</td>$ROW," $FILE_PATH
                done
            sed -i "s,\(long_cds'>[ A-Z]*</td>\)<\!--,\1<td class='long_pep'>PASS</td><\!--,g" $FILE_PATH
        else
            sed -i "s,<\!--ROW,<td class='long_pep'> </td><\!--ROW,g" $FILE_PATH
    fi
}
export -f REPORT_LONG
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'REPORT_LONG %'

# Add header column to gene report files and fill in zero counts in species and gene tables
cd $REPORT/genes
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % sed -i "s,\(<\!--THEAD-->\),<th>Long Branch CDS</th><th>Long branch PEP</th>\1," %

# Add Long Branch Analysis Figures
function LONG_FIGS {
    HTML_FILE=$1
    GENE=${HTML_FILE/.html/}
    RPATH=../figures/"loop_"$LOOP_NUMBER/long_branches
    HIST=$GENE".hist.pdf"
    TREE=$GENE".tre.pdf"
    if [ -f $LONG_FIG_PEP/$TREE ]
        then sed -i "s,\(<\!--FIGURES-->\),<h2>Loop $LOOP_NUMBER</h2>\n\t\1,
            s,\(<\!--FIGURES-->\),<h3>Long Branch CDS Histogram</h3>\n\t<img src=\"$RPATH/CDS/$HIST\" />\n\t\1,
            s,\(<\!--FIGURES-->\),<h3>Long Branch PEP Histogram</h3>\n\t<img src=\"$RPATH/PEP/$HIST\" />\n\t\1,
            s,\(<\!--FIGURES-->\),<h3>Long Branch CDS Tree</h3>\n\t<img src=\"$RPATH/CDS/$TREE\" />\n\t\1,
            s,\(<\!--FIGURES-->\),<h3>Long Branch PEP Tree</h3>\n\t<img src=\"$RPATH/PEP/$TREE\" />\n\t\1," $HTML_FILE
    fi
}
export -f LONG_FIGS
LONG_FIG_CDS=$REPORT/figures/"loop_"$LOOP_NUMBER/long_branches/CDS
LONG_FIG_PEP=$REPORT/figures/"loop_"$LOOP_NUMBER/long_branches/PEP
export LONG_FIG_PEP
mkdir -p $LONG_FIG_CDS $LONG_FIG_PEP
cp $WORKING/long_branches/CDS/*.pdf $LONG_FIG_CDS
cp $WORKING/long_branches/PEP/*.pdf $LONG_FIG_PEP
RFILE=($(find $REPORT/genes/*.html -type f | sed 's#.*/##'))
cd $REPORT/genes
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'LONG_FIGS %'

# Distance Matrix Analysis

# Generate raxml pairwise distance matrixes for each gene
printf "***************************   Generating distance matrixes for constrained trees ***************************\n"
echo "Starting RAxML distance matrix calculation on $(date)" >> $INPUT/log.txt
cd $WORKING/gene_trees/CDS
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f x -p 258755 -s %".phy" -n %".txt" -m GTRGAMMA -t %".constraint.tre"
cd $WORKING/gene_trees/PEP
printf "%s\n" "${FILE[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f x -p 258755 -s %".phy" -n %".txt" -m PROTGAMMAAUTO -t %".constraint.tre"
# Run analysis of distances to flag outliers for each sequence
printf "***************************   Running distance matrix analysis  ***************************\n"
echo "Starting distance_matrix_zscore.py on $(date)" >> $INPUT/log.txt
mkdir -p $WORKING/dist_m_zscore/CDS
mkdir -p $WORKING/dist_m_zscore/PEP
printf "***************************   CDS Files  ***************************\n"
echo "distance_matrix_zscore.py --input $WORKING/gene_trees/CDS/ --score 85 --tree $WORKING/gene_trees/CDS/ --out $WORKING/dist_m_zscore/CDS/ --outgroups $OUTGROUPS" >> $INPUT/log.txt
distance_matrix_zscore.py --input $WORKING/gene_trees/CDS/ --score 85 --tree $WORKING/gene_trees/CDS/ --out $WORKING/dist_m_zscore/CDS/ --outgroups $OUTGROUPS
printf "***************************   PEP Files  ***************************\n"
echo "distance_matrix_zscore.py --input $WORKING/gene_trees/PEP/ --score 90 --tree $WORKING/gene_trees/PEP/ --out $WORKING/dist_m_zscore/PEP/ --outgroups $OUTGROUPS" >> $INPUT/log.txt
distance_matrix_zscore.py --input $WORKING/gene_trees/PEP/ --score 90 --tree $WORKING/gene_trees/PEP/ --out $WORKING/dist_m_zscore/PEP/ --outgroups $OUTGROUPS

# Generate Report for Distance Matrix Analysis
RFILE=($(find $REPORT/genes/*.html -type f | sed 's#.*/##'))
function REPORT_DIST {
    HTML_FILE=$1
    FILE_PATH=$INPUT/report/genes/$HTML_FILE
    GENE=${HTML_FILE/.html/}
    DIST_CDS=$WORKING/dist_m_zscore/CDS/"outlier_taxa."$GENE".txt"
    DIST_PEP=$WORKING/dist_m_zscore/PEP/"outlier_taxa."$GENE".txt"
    if [ -f $DIST_CDS ]
        then
            readarray BAD_SPECIES_LIST < $DIST_CDS
            for SPECIES in ${BAD_SPECIES_LIST[@]}
                do
                    ROW=\<\!--ROW_"$SPECIES"--\>
                    sed -i "s,$ROW,<td class='fail dist_cds'>FAIL</td>$ROW," $FILE_PATH
                done
            sed -i "s,\(long_pep'>[ A-Z]*</td>\)<\!--,\1<td class='dist_cds'>PASS</td><\!--,g" $FILE_PATH
        else
            sed -i "s,<\!--ROW,<td class='dist_cds'>N/A</td><\!--ROW,g" $FILE_PATH
    fi
    if [ -f $DIST_PEP ]
        then
            readarray BAD_SPECIES_LIST < $DIST_PEP
            for SPECIES in ${BAD_SPECIES_LIST[@]}
                do
                    ROW=\<\!--ROW_"$SPECIES"--\>
                    sed -i "s,$ROW,<td class='fail dist_pep'>FAIL</td>$ROW," $FILE_PATH
                done
            sed -i "s,\(dist_cds'>[ A-Z]*</td>\)<\!--,\1<td class='dist_pep'>PASS</td><\!--,g" $FILE_PATH
        else
            sed -i "s,<\!--ROW,<td class='dist_pep'>N/A</td><\!--ROW,g" $FILE_PATH
    fi
}
export -f REPORT_DIST
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'REPORT_DIST %'

# Add header column to gene report files and fill in zero counts in species and gene tables
cd $REPORT/genes
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % sed -i "s,\(<\!--THEAD-->\),<th>Distance CDS</th><th>Distance PEP</th>\1," %

# Add Distance Matrix Figures
function DIST_FIGS {
    HTML_FILE=$1
    GENE=${HTML_FILE/.html/}
    RPATH=../figures/"loop_"$LOOP_NUMBER/dist_m_zscore
    HIST=$GENE".hist.pdf"
    TREE=$GENE".tre.pdf"
    if [ -f $DIST_FIG_PEP/$TREE ]
        then sed -i "s,\(<\!--FIGURES-->\),<h3>Distance Matrix CDS Tree</h3>\n\t<img src=\"$RPATH/CDS/$TREE\" />\n\t\1,
            s,\(<\!--FIGURES-->\),<h3>Distance Matrix PEP Tree</h3>\n\t<img src=\"$RPATH/PEP/$TREE\" />\n\t\1," $HTML_FILE
    fi
}
export -f DIST_FIGS
DIST_FIG_CDS=$REPORT/figures/"loop_"$LOOP_NUMBER/dist_m_zscore/CDS
DIST_FIG_PEP=$REPORT/figures/"loop_"$LOOP_NUMBER/dist_m_zscore/PEP
export DIST_FIG_PEP
mkdir -p $DIST_FIG_CDS $DIST_FIG_PEP
cp $WORKING/dist_m_zscore/CDS/*.pdf $DIST_FIG_CDS
cp $WORKING/dist_m_zscore/PEP/*.pdf $DIST_FIG_PEP
RFILE=($(find $REPORT/genes/*.html -type f | sed 's#.*/##'))
cd $REPORT/genes
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'DIST_FIGS %'

############
# inparalog Analysis
############
collapse_inparalogs.sh -i $INPUT -s $LOOP_DIR/lists/species.txt -g $LOOP_DIR/lists/genes.txt -t $THREADS -og $OUTGROUPS

##########################################################################
#     Generate new inputs for another round of the loop
##########################################################################

# Use the output from distance_matrix_z-score.py and long_branches.py to select species for new queries
printf "***************************   Building new query files ***************************\n"
echo "Starting rebuild_queries.sh on $(date)" >> $INPUT/log.txt
echo "rebuild_queries.sh -i $INPUT -t $THREADS" >> $INPUT/log.txt
rebuild_queries.sh -i $INPUT -t $THREADS
# Use new query sequences from rebuild_queries.sh to build hmms and generate bitscore cut off file
printf "***************************   Building new HMMs files and calculating bitscore cut offs   ***************************\n"
echo "Starting new_score_genes.sh on $(date)" >> $INPUT/log.txt
echo "new_score_genes.sh -i $LOOP_DIR/rebuilt_queries/query -o $LOOP_DIR/rebuilt_queries -e fas -t $THREADS" >> $INPUT/log.txt
new_score_genes.sh -i $LOOP_DIR/rebuilt_queries/query -o $LOOP_DIR/rebuilt_queries -e fas -t $THREADS
printf "***************************   DONE! Now run loop.sh again using the new queries.   ***************************\n"

printf "***************************   Creating Paralog Screened PEP fasta files   ***************************\n"
trim_dbs.sh -i $INPUT -s $LOOP_DIR/lists/species.txt -t $THREADS

printf "*************************************************************\n\n\t
\tAfter checking that the script has properly executed, you should start another\n
\trun of the main loop with the following command:\n\n
loop.sh -d $INPUT/CDS_degenerate_seqs -p $LOOP_DIR/screened_fasta -q $LOOP_DIR/rebuilt_queries -o $INPUT -t $THREADS
\n\n*************************************************************\n"
