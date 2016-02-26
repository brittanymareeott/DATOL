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

WORKING=$( find $INPUT -name 'Working_Dir_*' -exec basename {} \; )
mkdir -p $INPUT/$WORKING/rebuilt_queries/cat_lists
mkdir -p $INPUT/$WORKING/rebuilt_queries/query

function get_good() {
    cat $INPUT/$WORKING/long_branches/CDS/longbranch_taxa.$GENE.txt $INPUT/$WORKING/long_branches/PEP/longbranch_taxa.$GENE.txt $INPUT/$WORKING/dist_m_zscore/CDS/outlier_taxa.$GENE.txt $INPUT/$WORKING/dist_m_zscore/PEP/outlier_taxa.$GENE.txt > $INPUT/$WORKING/rebuilt_queries/cat_lists/$GENE.txt
    FILE=$INPUT/$WORKING/rebuilt_queries/cat_lists/$GENE.txt
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
    SEQS=$INPUT/$WORKING/subset_sorted/PEP/$GENE
    AVAILABLE_SPECIES=($(find $SEQS -type f -exec basename {} \; ))
    GOOD_SPECIES=()
    for AV_SPECIES in ${AVAILABLE_SPECIES[@]}
        do
            if [[ " ${BAD_SPECIES[*]} " != *" ${AV_SPECIES/.fas/} "* ]]
                then
                    GOOD_SPECIES+=("$AV_SPECIES")
            fi
        done
    cd $SEQS
    cat ${GOOD_SPECIES[@]} > $INPUT/$WORKING/rebuilt_queries/query/$GENE.faa
    cd $INPUT/$WORKING/rebuilt_queries/query/
    # This is necessary to make the formatting of these files match that of the original input sequences
    # so the new_score_genes.sh script will parse the def lines correctly.
    sed -i "s,\(>[0-9A-Za-z_.|]*\),\1___$GENE,g" $GENE".faa"
}

FILE_DIR=$INPUT/$WORKING/subset_sorted/PEP
GENES=($(find $FILE_DIR -not -name 'PEP' -type d -exec basename {} \; ))
export -f get_good
export INPUT
export WORKING
printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'GENE=% ;\
    get_good'