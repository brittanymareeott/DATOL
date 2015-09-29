#!/usr/bin/env bash

# Named variables. Every run needs the following defined:
# 1) -c | --cutoff_file - The file defining cut off values for each gene
# 2) -h | --hmm_dir - The directory containing the HMMs for each gene
# 3) -in | --input_dir - The directory with the Peptide sequences produced by the previous step in the pipeline
# 4) -out | --output_dir - The directory to put the output
# 5) -pre | --prefix - The prefix of the gene names. Examples: KOG, CHLOG, BUSCO, NOG
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
  -in|--input_dir)
  IN_DIR="$2"
  shift # past argument
  ;;
  -out|--output_dir)
  OUT_DIR="$2"
  shift # past argument
  ;;
  -c|--cutoff_file)
  CUTOFF_FILE="$2"
  shift # past argument
  ;;
  -h|--hmm_dir)
  HMM_DIR="$2"
  shift # past argument
  ;;
  -pre|--prefix)
  PREFIX="$2"
  shift # past argument
  ;;
  *)
        # unknown option
  ;;
esac
shift # past argument or value
done
echo hmm_from_pep.sh run on $(date) >> $OUT_DIR/log.txt
echo Input Directory = "$IN_DIR" >> $OUT_DIR/log.txt
echo Output Directory = "$OUT_DIR" >> $OUT_DIR/log.txt
echo Cut off file = "$CUTOFF_FILE" >> $OUT_DIR/log.txt
echo HMM Directory = "$HMM_DIR" >> $OUT_DIR/log.txt
echo Gene Prefix = "$PREFIX" >> $OUT_DIR/log.txt

# Loop through fasta files
for i in $IN_DIR/query*.faa
  do
# Set some variable:
# We find the GENE name in the fasta file with some regular expressions. Specifically we are search for our prefix followed by numbers.
# We set HMM to the gene name plus the .hmm extension.
# We set the output file name (OUT) to the current file name minus the query_, .faa., and _db_. This simplifys the output names.
# We set the cutoff score by finding the name in the cutoff file matching GENE and a space then we remove the gene name so only the number remains.
    GENE=`echo "${i##*/}" | sed -e 's,query_\('$PREFIX'[0-9]*\).*,\1,'`
    HMM=$GENE.hmm
    OUT=`echo "${i##*/}" | sed -e 's,query_,,' -e 's,.faa,,g' -e 's,_db,,'`.out
    CUTOFF=`sed -n /$GENE/p $CUTOFF_FILE | sed -e 's,'$PREFIX'[0-9]* ,,'`
# With all those variables set we run our hmmsearch
    hmmsearch --tblout $OUT_DIR/$OUT -T $CUTOFF $HMM_DIR/$HMM $i
done
