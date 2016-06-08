#!/usr/bin/env bash
#
# final_check.sh
#
# Author: Gregg Mendez

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
REPORT=$INPUT/report
export REPORT
export WORKING
export LOOP_DIR
mkdir -p $LOOP_DIR/lists

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
export -f ARRAY_CONTAINS

# First we need to rerun collapse_inparalogs. This will also run
# distance_matrix_zscore.
# collapse_inparalogs.sh -i $INPUT -s $SPECIES_LIST -g $GENE_LIST -t $THREADS -og $OUTGROUPS_FILE

# Next we need to generate some counts of the paralog hits for each gene and see
# if any genes fall below the 80% cut off for number of species with good hits.

# Count number of species in species list and genes in gene lists
SPECIES_NUMBER=$(wc -l < $SPECIES_LIST)
GENE_NUMBER=$(wc -l < $GENE_LIST)
export SPECIES_NUMBER
export GENE_NUMBER
export SPECIES_LIST

#LOAD input files into arrays; We'll have to do the species list later since bash can't export arrays into the environment
GENES=($(cat $GENE_LIST | sort))

function EVAL_GENE {
    # Get a count on the number of species prior to collapse in paralogs
    SPECIES_PRESENT=$(find $WORKING/early_subset_sorted/CDS/$GENE/*.fas | wc -l)
    # Now get count of how many were flagged as paralogs by distance matrix
    if [ -f $WORKING/dist_m_zscore/CDS/outlier_taxa.$GENE".txt" ]; then
        DIST_PARA=$(wc -l < $WORKING/dist_m_zscore/CDS/outlier_taxa.$GENE".txt")
        PARA_SPECIES=($(cat $WORKING/dist_m_zscore/CDS/outlier_taxa.$GENE".txt" ))
        RATIO_DIST_PARA=$( echo "$DIST_PARA / $SPECIES_PRESENT" | bc -l )
        SPECIES_PRESENT=$(($SPECIES_PRESENT - $DIST_PARA))
    else
        RATIO_DIST_PARA=0
    fi
    # Number of hits. We start by assuming at least 1 for every species present
    HIT_COUNTER=$SPECIES_PRESENT

    # per species
    if [ -d $WORKING/inparalogs/$GENE ]; then
        SPECIES_LIST_ARRAY=($(cat $SPECIES_LIST | sort))
        for SPECIES in ${SPECIES_LIST_ARRAY[@]}; do
            if [ $( ARRAY_CONTAINS PARA_SPECIES $SPECIES && echo yes || echo no ) == no  ]; then
                if [ -f $WORKING/inparalogs/$GENE/out/$SPECIES"___paralogs.txt" ]; then
                    PARA_FILE=$WORKING/inparalogs/$GENE/out/$SPECIES"___paralogs.txt"
                    LABELS_FILE=$WORKING/inparalogs/$GENE/$SPECIES"___labels.txt"
                    START_NUMBER=$(wc -l < $LABELS_FILE)
                    NUMBER=$(wc -l < $PARA_FILE)
                    NEW_TOTAL=$(($START_NUMBER - $NUMBER))
                    # if our NEW_TOTAL is negative then it means the original
                    # sequence was marked as a paralog because all its sisters were
                    # paralogs.
                    if [ $NEW_TOTAL == "-1" ]; then
                        SPECIES_PRESENT=$(($SPECIES_PRESENT - 1))
                    fi
                    HIT_COUNTER=$(($HIT_COUNTER + $NEW_TOTAL))
                fi
            fi
        done
    fi

    # A Gene must pass all 3 following checks to be recommended
    PERC_REMAINING=$( echo "($SPECIES_PRESENT) / $SPECIES_NUMBER" | bc -l)
    # Gene must have more than 75% of the species present
    if (( $( echo "$PERC_REMAINING > 0.75" | bc -l ) )); then
        PERC_REMAINING_CHECK="PASS"
        if [ $SPECIES_PRESENT != "0" ]; then
            PERC_HITS=$( echo "$HIT_COUNTER / $SPECIES_PRESENT" | bc -l)
            # Gene have less than 1.1 hits per species
            if (( $( echo "$PERC_HITS < 1.1" | bc -l ) )); then
                PERC_HITS_CHECK="PASS"
                # Gene must have less than 15% of its top hits marked as paralogs
                if (( $( echo "$RATIO_DIST_PARA < 0.15" | bc -l ) )); then
                    RATIO_DIST_PARA_CHECK="PASS"
                    printf "%s\n" "$GENE" >> $LOOP_DIR/lists/revised_genes.txt
                else
                    RATIO_DIST_PARA_CHECK="FAIL"
                fi
            else
                PERC_HITS_CHECK="FAIL"
            fi
        else
            SPECIES_PRESENT_CHECK="FAIL"
            PERC_HITS=0
        fi
    else
        PERC_REMAINING_CHECK="FAIL"
    fi
    # Now to report results in html files
    FINAL_CHECK_TABLE="<table><thead><tr><th>Test</th><th>Value</th><th>Pass/Fail</th></tr></thead><tbody><tr><td>Percent Species Remaining (75%)</td><td>$PERC_REMAINING</td><td>$PERC_REMAINING_CHECK</td></tr><tr><td>Average Number Hits (1.1)</td><td>$PERC_HITS</td><td>$PERC_HITS_CHECK</td></tr><tr><td>Percent Top Hit Paralogous (15%)</td><td>$RATIO_DIST_PARA</td><td>$RATIO_DIST_PARA_CHECK</td></tr></tbody></table>"
    sed -i "s,\(<\!--SPECIES TABLE-->\),$FINAL_CHECK_TABLE\n\1," $REPORT/genes/$GENE".html"
}
export -f EVAL_GENE

printf "%s\n" "${GENES[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'GENE=%; \
    EVAL_GENE'

printf "*************************************************************\n\n\t
\tReview the gene reports and the revised gene list written to\n\t $LOOP_DIR/lists/revised_genes.txt. Make any changes you wish then run\n\t pretree_loop.sh as follows:\n\n
pretree_loop.sh -i $INPUT -s $LOOP_DIR/lists/species.txt $LOOP_DIR/lists/revised_genes.txt -p $WORKING/dist_m_zscore/CDS -t $THREADS

\n\n*************************************************************\n"
