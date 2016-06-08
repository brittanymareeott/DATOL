#!/usr/bin/env bash
#
# rebuild_queries.sh
#
# Author: Gregory Mendez
#
# For each gene, this script takes lists of outlier sequences from long_branches.py and distance_matrix_zscore.py
# combines them into a master list used to specify which sequences will be used to generate new query files. After
# this, one should run new_score_genes.sh to generate HMMs and create a list of HMMsearch bitscore cut off scores.
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

# find working directory name
LOOP_NUMBER=$(find $INPUT -path "$INPUT/loop*" -prune | wc -l )
LOOP_DIR=$INPUT/"loop_"$LOOP_NUMBER"_out/"
WORKING=$INPUT/"loop_"$LOOP_NUMBER"_out/tmp"
echo Working Directory = "$WORKING" >> $INPUT/log.txt
mkdir -p $WORKING/rebuilt_queries/cat_lists
mkdir -p $LOOP_DIR/rebuilt_queries/query

function get_good() {
    cat $WORKING/long_branches/CDS/longbranch_taxa.$GENE.txt $WORKING/long_branches/PEP/longbranch_taxa.$GENE.txt $WORKING/dist_m_zscore/CDS/outlier_taxa.$GENE.txt $WORKING/dist_m_zscore/PEP/outlier_taxa.$GENE.txt > $WORKING/rebuilt_queries/cat_lists/$GENE.txt
    FILE=$WORKING/rebuilt_queries/cat_lists/$GENE.txt
    #for each gene do
    BAD_SPECIES=()
    while read LINE
        do
            SPECIES=$LINE
            if [[ " ${BAD_SPECIES[*]} " != *" $SPECIES "* ]]
                then
                    BAD_SPECIES+=("$SPECIES")
            fi
        done < $FILE
    SEQS=$LOOP_DIR/sequences/TopHits/PEP/$GENE
    AVAILABLE_SPECIES=($(find $SEQS -type f | sed 's#.*/##' ))
    GOOD_SPECIES=()
    for AV_SPECIES in ${AVAILABLE_SPECIES[@]}
        do
            if [[ " ${BAD_SPECIES[*]} " != *" ${AV_SPECIES/.fas/} "* ]]
                then
                    GOOD_SPECIES+=("$AV_SPECIES")
            fi
        done
    cd $SEQS
    if [[ ${#GOOD_SPECIES[@]} < 1 ]]; then
        printf "%s\n" "$GENE" >> $WORKING/rebuilt_queries/no_good_species.txt
    else
        cat ${GOOD_SPECIES[@]} > $LOOP_DIR/rebuilt_queries/query/$GENE.fas
        cd $LOOP_DIR/rebuilt_queries/query/
    # This is necessary to make the formatting of these files match that of the original input sequences
    # so the new_score_genes.sh script will parse the def lines correctly.
        sed -i "s,\(>[0-9A-Za-z_.|]*\),\1___$GENE,g" $GENE".fas"
    fi
}

FILE_DIR=$LOOP_DIR/tmp/subset_sorted/CDS
GENES=($(find $FILE_DIR -not -name 'CDS' -type d | sed 's#.*/##' ))
export -f get_good
export INPUT
export WORKING
export LOOP_DIR
printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'GENE=% ; \
    get_good'
