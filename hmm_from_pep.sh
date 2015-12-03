#!/usr/bin/env bash

# Named variables. Every run needs the following defined:
# 1) -c | --cutoff_file - The file defining cut off values for each gene
# 2) -hmm | --hmm_dir - The directory containing the HMMs for each gene
# 3) -i | --input_dir - The directory with the Peptide sequences produced by the previous step in the pipeline
# 4) -o | --output_dir - The directory to put the output
# Example:
# hmm_from_pep.sh -c /home/mendezg/cegma/cutoff.txt -h /home/mendezg/cegma/hmm -i /home/mendezg/cegma_dinos/pepsfromblast -o /home/mendezg/cegma_dinos/hmmsearch_out -pre KOG

#Common  locations:
#cegma dino cleanup
#cutoff_file=/home/mendezg/Dropbox/Documents/UMD/genomes/cegma_cleanup/cat_genes/cutoff.txt
#hmm_dir=/home/mendezg/Dropbox/Documents/UMD/genomes/cegma_cleanup/cat_genes/hmm_profiles

#cegma locations
#cutoff_file=$datol/data/profiles_cutoff.tbl
#hmm_dir=$datol/data/hmm_profiles

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
  -c|--cutoff_file)
  CUTOFF_FILE="$2"
  shift # past argument
  ;;
  -hmm|--hmm_dir)
  HMM_DIR="$2"
  shift # past argument
  ;;
  *)
        # unknown option
  ;;
esac
shift # past argument or value
done
echo hmm_from_pep.sh run on $(date) >> $OUT/log.txt
echo Input Directory = "$INPUT" >> $OUT/log.txt
echo Output Directory = "$OUT" >> $OUT/log.txt
echo Cut off file = "$CUTOFF_FILE" >> $OUT/log.txt
echo HMM Directory = "$HMM_DIR" >> $OUT/log.txt

function FIND_SPECIES_GENE()
{
        SPLIT_FILE=($(echo $1 | tr "\." "\n" | tr "_" "\n"))
        unset SPLIT_FILE[${#SPLIT_FILE[@]}-1]
        GENE=${SPLIT_FILE[0]}
        SPECIES_SPACED=${SPLIT_FILE[@]:1}
        SPECIES=${SPECIES_SPACED// /_}
        OUT_FILE=$GENE'_'$SPECIES'.out'
}

# Loop through fasta files
for i in $INPUT/*.faa
  do
# Set some variable:
# We find the GENE name in the fasta file with some regular expressions. Specifically we are search for our prefix followed by numbers.
# We set HMM to the gene name plus the .hmm extension.
# We set the output file name (OUT) to the current file name minus the query_, .faa., and _db_. This simplifys the output names.
# We set the cutoff score by finding the name in the cutoff file matching GENE and a space then we remove the gene name so only the number remains.
    HMM=$GENE.hmm
    FIND_SPECIES_GENE ${i##*/}
    CUTOFF=`sed -n /$GENE/p $CUTOFF_FILE | sed -e 's,'$GENE' ,,'`
# With all those variables set we run our hmmsearch
    hmmsearch --tblout $OUT/$OUT_FILE -T $CUTOFF $HMM_DIR/$HMM $i
done