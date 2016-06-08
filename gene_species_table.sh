#!/usr/bin/env bash

# This script takes as input a directory, like that created by cleanup_explode_sort.sh,
# with directories for each gene with individual fasta files within each gene directory
# named for the species the sequence came from. It outputs a table in .csv format that
# displays 0/1 to indicate the absence or presence of a gene (column) for a species
# (row). I generally import this into Excel, transpose the table, and add percentage
# calculations for each row and colum. Based on those data I determine which species and
# genes should move forward in the pipeline by removing genes and species with a lot of
# gaps.
# Example Input:
# Root:
#     Gene001:
#         Homo_sapiens.fas
#         Rattus_rattus.fas
#         Procyon_lotor.fas
#         Canis_domesticus.fas
#     Gene002:
#         Homo_sapiens.fas
#         Procyon_lotor.fas
#         Canis_domesticus.fas
#     Gene003:
#         Homo_sapiens.fas
#         Rattus_rattus.fas
#         Procyon_lotor.fas
#         Canis_domesticus.fas
#
# This script takes two arguments:
# 1) -i | --input_dir - The input directory as discussed above.
# 2) -o | --output_dir - The directory to write the final .csv file and tmp files. WARNING: This directory must not contain files ending in .out.
#
# Example Usage: gene_species_table.sh -i /home/mendezg/crawly_things/all_sorted -o /home/mendezg/crawly_things/gene_coverage_table
#
# Some would say this should have been done in python or perl. Some would say I should have used more arrays/hash tables.
# I find this solution oddly amusing though.

#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
    -i|--input_dir)
    INPUT="$2"
    shift # past argument
    ;;
    -o|--output_dir)
    OUT="$2"
    shift # past argument
    ;;
    *)
         # unknown option
    ;;
esac
shift # past argument or value
done

mkdir -p $OUT
mkdir -p $OUT/tmp/gene_species_table
rm $OUT/tmp/gene_species_table/* 2> /dev/null
rm $OUT/big_fat_table.csv 2> /dev/null
cd $INPUT
SPECIES=('Gene ID')
SPECIES+=( $(find . -name '*.fas' -type f | sed 's#.*/##' | sort | uniq | sed 's,.fas,,') )
(IFS=",$IFS"; printf '%s\n' "${SPECIES[*]}" >> $OUT/tmp/gene_species_table/species_names.txt)
for GENE in *
    do
        cd $GENE
        ARRAY=("$GENE")
        for TAXON in ${SPECIES[@]:1}
            do
                if [ -e $TAXON.fas ]
                    then
                        ARRAY+=('1')
                    else
                        ARRAY+=('0')
                fi
            done
        (IFS=",$IFS"; printf '%s\n' "${ARRAY[*]}" >> $OUT/tmp/gene_species_table/$GENE.out)
        cd $INPUT
    done
cat $OUT/tmp/gene_species_table/species_names.txt $OUT/tmp/gene_species_table/*.out > $OUT/gene_species_table.csv
rm $OUT/tmp/gene_species_table/*
