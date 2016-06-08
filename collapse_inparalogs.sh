#!/usr/bin/env bash
#
# collapse_inparalogs.sh
#
# Author Gregory S Mendez

# Given lists of species and genes this script will generate a list of sequences
# that are inparalogs and sister to parologous sequences listed in another file.
# This will take a long time to run since it needs to find trees for every
# species-gene combo for which there are more than one ortholog found.


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
    -og|--outgroups)
    OUTGROUPS_FILE="$2"
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

# Find working directory name
LOOP_NUMBER=$(find $INPUT -path "$INPUT/loop*" -prune | wc -l )
LOOP_DIR=$INPUT/"loop_"$LOOP_NUMBER"_out"
WORKING=$INPUT/"loop_"$LOOP_NUMBER"_out/tmp"

# Make directory for inparalog analysis
mkdir -p $WORKING/inparalogs
INP_DIR=$WORKING/inparalogs
OTHER_DIR=$LOOP_DIR/sequences/OtherHits/CDS
PARALOGS_DIR=$WORKING/dist_m_zscore/CDS
REPORT=$INPUT/report

# set alignment directory depending on loop NUMBER
if [ $LOOP_NUMBER -eq 1 ]; then
    ALI_DIR=$WORKING/cat_mafft_all/CDS
else
    ALI_DIR=$WORKING/early_gene_trees
    TRIMAL_OUT=$WORKING/early_gene_trees/trimal
    mkdir -p $TRIMAL_OUT/CDS
    mkdir -p $ALI_DIR/trees/CDS
fi

export TRIMAL_OUT
export REPORT
export INP_DIR
export OTHER_DIR
export ALI_DIR
export PARALOGS_DIR
export OUTGROUPS_FILE
export WORKING
export INPUT
export LOOP_NUMBER

################################
## FUNCTIONS
################################
# Function split input on input character
# use SPLIT "0" "___" species_genus_genbank___weird_trinity_seq_id
function SPLIT()
{
        SPLIT_LIST=($(echo $3 | sed -e "s,$2,\n,"))
        echo ${SPLIT_LIST[$1]}
}
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

# function to generate list of paralogs names
GET_LABELS () {
    while read LINE
        do
            if [[ $LINE =~ ^\> ]]
                then
                TAXON=$(SPLIT 0 ___ ${LINE:1})
                SEQ_ID=$(SPLIT 1 ___ ${LINE})
                    if [ $( ARRAY_CONTAINS SPECIES $TAXON && echo yes || echo no ) == yes  ]
                        then
                        printf "%s___%s\n" "$TAXON" "$SEQ_ID" >> $TAXON"___labels.txt"
                    fi
            fi
        done < $OTHER_DIR/$1".fasta"
}
# Function to use paralog name lists to create fasta file
GET_SEQS () {
    TAXON=$(SPLIT 0 ___ $1)
    usearch -fastx_getseqs $2 -labels $1 -fastaout $TAXON"___subset.fa"
}
# Function to do the loops
EACH_GENE () {
    if [ -f $OTHER_DIR/$GENE".fasta" ]
        then
        SPECIES=($(cat $SPECIES_LIST))
        mkdir $INP_DIR/$GENE
        cd $INP_DIR/$GENE
        GET_LABELS $GENE
        for SPECIES_FILE in *___labels.txt
            do
                GET_SEQS $SPECIES_FILE $OTHER_DIR/$GENE".fasta"
            done
    fi
}
# Function to build constraint trees
STARTING_TREES () {
    if [ -d $GENE ]; then
        cd $INP_DIR/$GENE
        cp $ALI_DIR/$GENE".fas" $GENE".fas"
        raxmlHPC-PTHREADS-SSE3 -p 558962 -m GTRGAMMA -s $GENE".fas" -n $GENE".tre" -T 2
        echo "$GENE already found"
    else
        if [ $LOOP_NUMBER > 1 ]; then
            mkdir $INP_DIR/$GENE
            cd $INP_DIR/$GENE
            cp $ALI_DIR/$GENE".fas" $GENE".fas"
            raxmlHPC-PTHREADS-SSE3 -p 558962 -m GTRGAMMA -s $GENE".fas" -n $GENE".tre" -T 2
            cp "RAxML_bestTree."$GENE".tre" $ALI_DIR/trees/CDS/"RAxML_result."$GENE".constrained.tre"
            cd ../
            rm -r $INP_DIR/$GENE
        fi
    fi
}
# Function to add Other sequences to alignments
ADD_OTHERS () {
    if [ -d $GENE ]; then
        cd $INP_DIR/$GENE
        for SUB_FILE in *___subset.fa
            do
                SPECIES=$(SPLIT 0 ___ ${SUB_FILE})
                mafft --quiet --add $SUB_FILE --reorder $GENE".fas" > $SPECIES"___others.fas"
            done
    fi
}
# Function to add make trees with other hits
FIND_OTHERS_TREES () {
    if [ -d $GENE ]; then
        cd $INP_DIR/$GENE
        for SUB_FILE in *___others.fas
            do
                SPECIES=$(SPLIT 0 ___ ${SUB_FILE})
                raxmlHPC-PTHREADS-SSE3 -p 558962 -m GTRGAMMA -s $SUB_FILE -n $SPECIES".tre" -T 2 -r RAxML_bestTree.$GENE.tre
            done
    fi
}
FIND_PARALOGS () {
    if [ -d $GENE ]; then
        cd $INP_DIR/$GENE
        # echo $GENE
        mkdir out
        for SUB_FILE in *___labels.txt
            do
                SPECIES=$(SPLIT 0 ___ ${SUB_FILE})
                echo $SPECIES
                TREE_FILE="RAxML_bestTree."$SPECIES".tre"
                paralogs.py --tree $TREE_FILE --para $PARALOGS_DIR/"outlier_taxa."$GENE".txt" --others $SUB_FILE --out out --outgroups $OUTGROUPS_FILE
            done
    fi
}

export -f FIND_PARALOGS
export -f FIND_OTHERS_TREES
export -f ADD_OTHERS
export -f STARTING_TREES
export -f GET_LABELS
export -f GET_SEQS
export -f ARRAY_CONTAINS
export -f SPLIT
export -f EACH_GENE
export SPECIES_LIST

#LOAD input files into arrays; We'll have to do the species list later since bash can't export arrays into the environment
GENES=($(cat $GENE_LIST))

# First we get the other sequences ready
printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'GENE=%; \
    EACH_GENE'

# If after initial loop, run cat_mafft_all
if [ $LOOP_NUMBER > 1 ]; then
    # Move desired sequences specified by gene and species lists to a new directory
    mkdir $WORKING/early_subset_sorted
    SUBSET_SORTED_OUT=$WORKING/early_subset_sorted
    SUBSET_SORTED_INPUT=$LOOP_DIR/sequences/TopHits
    printf "***************************   Copying CDS files of genes and species specified ***************************\n"
    echo "Starting subset_sorted.sh CDS run on $(date)" >> $INPUT/log.txt
    echo "subset_sorted.sh -i $SUBSET_SORTED_INPUT/CDS -o $SUBSET_SORTED_OUT/CDS -s $SPECIES_LIST -g $GENE_LIST" >> $INPUT/log.txt
    subset_sorted.sh -i $SUBSET_SORTED_INPUT/CDS -o $SUBSET_SORTED_OUT/CDS -s $SPECIES_LIST -g $GENE_LIST
    echo "Completed subset_sorted.sh CDS run on $(date)" >> $INPUT/log.txt
    # then run cat_mafft_all
    cat_mafft_all.sh -i $SUBSET_SORTED_OUT/CDS -o $ALI_DIR -e .fas -t $THREADS
fi

# Now we generate our starting gene trees
cd $INP_DIR
RUNS=$( echo "scale=0;$THREADS/2" | bc -l )
printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $RUNS -I % bash -c 'GENE=%; \
    STARTING_TREES'

# If not on initial loop, run distance_matrix_zscore. This requires trimming our alignments, generating distance matrixes, and then doing the analysis.
if [ $LOOP_NUMBER > 1 ]; then
    # First Run TrimAL
    cd $ALI_DIR
    printf "***********   Starting trimal for CDS sequences on `date` ...\n"
    echo "Starting trimal CDS run on $(date)" >> $INPUT/log.txt
    echo "find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1" >> $INPUT/log.txt
    find *.fas | sed 's,.fas,,' | xargs -n 1 -P $THREADS -I % trimal -nexus -in %'.fas' -out $TRIMAL_OUT/CDS/%'.nex' -automated1
    cd $TRIMAL_OUT/CDS/
    # trimal seems to have a bug where it lists all nexus output as protein data, so here we change it to dna
    echo "Correcting trimal bug that lists all nexus data as protein data on $(date)" >> $INPUT/log.txt
    echo "find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %" >> $INPUT/log.txt
    find *.nex | xargs -n 1 -P $THREADS -I % sed -i 's,PROTEIN,DNA,' %
    # Convert nexus files to phylip
    printf "***************************   Converting Nexus files to Phylip ***************************\n"
    echo "Starting dir_nexus_to_phylip.sh CDS run on $(date)" >> $INPUT/log.txt
    echo "dir_nexus_to_phylip.sh -i $TRIMAL_OUT/CDS -o $ALI_DIR/trees/CDS -e nex -t $THREADS" >> $INPUT/log.txt
    dir_nexus_to_phylip.sh -i $TRIMAL_OUT/CDS -o $ALI_DIR/trees/CDS -e nex -t $THREADS
    # Generate raxml pairwise distance matrixes for each gene
    # We need to first copy the starting trees all to one directory.
    # We need to rename our trees First
    mkdir $ALI_DIR/trees/CDS
    printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'GENE=%; \
    find $INP_DIR -name "RAxML_bestTree."$GENE".tre"  -exec cp {} $ALI_DIR/trees/CDS/"RAxML_result."$GENE".constrained.tre" \;'
    printf "***************************   Generating distance matrixes for constrained trees ***************************\n"
    echo "Starting RAxML distance matrix calculation on $(date)" >> $INPUT/log.txt
    cd $ALI_DIR/trees/CDS
    printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $RUNS -I % raxmlHPC-PTHREADS-SSE3 -T 2 -f x -p 258755 -s %".phy" -n %".txt" -m GTRGAMMA -t "RAxML_result."%".constrained.tre"
    # Run analysis of distances to flag outliers for each sequence
    printf "***************************   Running distance matrix analysis  ***************************\n"
    echo "Starting distance_matrix_zscore.py on $(date)" >> $INPUT/log.txt
    mkdir -p $WORKING/dist_m_zscore/CDS
    printf "***************************   CDS Files  ***************************\n"
    echo "distance_matrix_zscore.py --input $ALI_DIR/trees/CDS  --score 85 --tree $ALI_DIR/trees/CDS --out $WORKING/dist_m_zscore/CDS/ --outgroups $OUTGROUPS_FILE" >> $INPUT/log.txt
    distance_matrix_zscore.py --input $ALI_DIR/trees/CDS/ --score 85 --tree $ALI_DIR/trees/CDS/ --out $WORKING/dist_m_zscore/CDS/ --outgroups $OUTGROUPS_FILE
fi

# We add the other sequences to our alignments
cd $INP_DIR
printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'GENE=%; \
    ADD_OTHERS'

# We find trees with the other sequences added
cd $INP_DIR
printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $RUNS -I % bash -c 'GENE=%; \
    FIND_OTHERS_TREES'

# Now we analyze those trees to collapse inparalogs
cd $INP_DIR
printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'GENE=%; \
    FIND_PARALOGS'

##################
# Generate Report
##################
function REPORT_PARA {
    HTML_FILE=$1
    FILE_PATH=$INPUT/report/genes/$HTML_FILE
    GENE=${HTML_FILE/.html/}
    FIGURES=$REPORT/figures
    if [ $LOOP_NUMBER > 1 ]; then
        # Generate Report for Distance Matrix Analysis
        DIST_CDS=$WORKING/dist_m_zscore/CDS/"outlier_taxa."$GENE".txt"
        if [ -f $DIST_CDS ]
            then
                readarray BAD_SPECIES_LIST < $DIST_CDS
                for SPECIES in ${BAD_SPECIES_LIST[@]}
                    do
                        ROW=\<\!--ROW_"$SPECIES"--\>
                        sed -i "s,$ROW,<td class='fail dist_cds'>FAIL</td>$ROW," $FILE_PATH
                    done
                sed -i "s,\(hmm'>[0-9]*</td>\)<\!--,\1<td class='dist_cds'>PASS</td><\!--,g" $FILE_PATH
            else
                sed -i "s,\(hmm'>[0-9]*</td>\)<\!--,\1<td class='dist_cds'>PASS</td><\!--,g" $FILE_PATH
        fi
        # Add header column to gene report files and fill in zero counts in species and gene tables
        sed -i "s,\(<\!--THEAD-->\),<th>Distance CDS</th>\1," $FILE_PATH
        # Add Distance Matrix Figures
        RPATH=../figures/"loop_"$LOOP_NUMBER/dist_m_zscore
        HIST=$GENE".hist.pdf"
        TREE=$GENE".tre.pdf"
        DIST_FIG_CDS=$REPORT/figures/"loop_"$LOOP_NUMBER/dist_m_zscore/CDS
        export DIST_FIG_PEP
        mkdir -p $DIST_FIG_CDS
        cp $WORKING/dist_m_zscore/CDS/*.pdf $DIST_FIG_CDS

        if [ -f $DIST_FIG_CDS/$TREE ]
            then sed -i "s,\(<\!--FIGURES-->\),<h3>Loop $LOOP_NUMBER Distance Matrix CDS Tree</h3>\n\t<img src=\"$RPATH/CDS/$TREE\" />\n\t\1," $FILE_PATH
        fi

    fi
    if [ -d $WORKING/inparalogs/$GENE ]
        then
            cd $WORKING/inparalogs/$GENE/out/
            mkdir -p $FIGURES/"loop_"$LOOP_NUMBER/inparalogs/$GENE
            for FILE in *___paralogs.txt
                do
                    SPECIES=${FILE/___paralogs.txt/}
                    PARA_FILE=$WORKING/inparalogs/$GENE/out/$SPECIES"___paralogs.txt"
                    LABELS_FILE=$WORKING/inparalogs/$GENE/$SPECIES"___labels.txt"
                    START_NUMBER=$(wc -l < $LABELS_FILE)
                    NUMBER=$(wc -l < $PARA_FILE)
                    NEW_TOTAL=$(($START_NUMBER - $NUMBER))
                    ROW=\<\!--ROW_"$SPECIES"--\>
                    sed -i "s,$ROW,<td class='paralog'>$NEW_TOTAL</td>$ROW," $FILE_PATH
                    # Figures
                    TREE=$SPECIES".paralog_tree.pdf"
                    cp $TREE $FIGURES/"loop_"$LOOP_NUMBER/inparalogs/$GENE/
                    RPATH=../figures/"loop_"$LOOP_NUMBER/inparalogs/$GENE
                    sed -i "s,\(<\!--FIGURES-->\),<h3>Loop $LOOP_NUMBER $SPECIES Paralog Gene Tree</h3>\n\t<img src=\"$RPATH/$TREE\" />\n\t\1," $FILE_PATH
                done
            if [ $LOOP_NUMBER == 1 ]; then
            sed -i "s,\(dist_pep'>[ A-Z]*</td>\)<\!--,\1<td class='paralog'>0</td><\!--,g" $FILE_PATH
            else
                sed -i "s,\(dist_cds'>[ A-Z]*</td>\)<\!--,\1<td class='paralog'>0</td><\!--,g" $FILE_PATH
            fi
        else
            if [ $LOOP_NUMBER == 1 ]; then
            sed -i "s,\(dist_pep'>[ A-Z/]*</td>\)<\!--,\1<td class='paralog'>N/A</td><\!--,g" $FILE_PATH
            fi
    fi
}

export -f REPORT_PARA
RFILE=($(find $REPORT/genes/*.html -type f | sed 's#.*/##'))
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'REPORT_PARA %'

# Add header column to gene report files and fill in zero counts in species and gene tables
cd $REPORT/genes
printf "%s\n" "${RFILE[@]}" | xargs -n 1 -P $THREADS -I % sed -i "s,\(<\!--THEAD-->\),<th>Other Hits</th>\1," %
