#!/usr/bin/env bash
#
# trim_dbs.sh
#
# Author: Gregory Mendez
#
# This script removes sequences from fasta files that were found to be paraogs
# according to the distance matrix analysis and inparalog analysis.
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
    SPECIES_FILE="$2"
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
mkdir -p $WORKING/trim_dbs/
mkdir -p $LOOP_DIR/screened_fasta/
TRIM_DBS_DIR=$WORKING/trim_dbs
SCREENED_DIR=$LOOP_DIR/screened_fasta
DIST_CDS_DIR=$WORKING/dist_m_zscore/CDS
DIST_PEP_DIR=$WORKING/dist_m_zscore/PEP
INP_DIR=$WORKING/inparalogs
ORIGINAL_FASTA_DIR=$LOOP_DIR/PEP_dereplicated

export SCREENED_DIR
export TRIM_DBS_DIR
export DIST_CDS_DIR
export DIST_PEP_DIR
export INP_DIR
export SPECIES_LIST
export ORIGINAL_FASTA_DIR

# Something like this:
# usearch -fastx_getseqs Amphidinium_carterae_MMET.fasta -labels remove.txt -notmatched subset.fas

function TRIM_DBS() {
    # Combine PEP and CDS Dist Files with all the inparalog analysis files
    # into a non-redundant file
    find $INP_DIR -name $SPECIES"___paralogs.txt" -exec cat {} > $TRIM_DBS_DIR/$SPECIES"___inp_all.txt" \;
    cat $TRIM_DBS_DIR/$SPECIES"___inp_all.txt" $DIST_PEP_DIR/$SPECIES"_seqids.txt" $DIST_CDS_DIR/$SPECIES"_seqids.txt" | sort | uniq > $TRIM_DBS_DIR/$SPECIES"___remove_all.txt"

    # Now use that file to remove those sequences from the fasta file used
    # by the loop.
    for RFILE in $TRIM_DBS_DIR/*___remove_all.txt
        do
            usearch -fastx_getseqs $ORIGINAL_FASTA_DIR/$SPECIES.fasta -labels $TRIM_DBS_DIR/$SPECIES"___remove_all.txt" -notmatched $SCREENED_DIR/$SPECIES".fasta"
        done
}

SPECIES_LIST=($(cat $SPECIES_FILE))
export -f TRIM_DBS
printf "%s\n" "${SPECIES_LIST[@]}" | xargs -n 1 -P $THREADS -I % bash -c 'SPECIES=% ; \
    TRIM_DBS'
