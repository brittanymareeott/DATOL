#!/usr/bin/env bash

# Named variables. Every run needs the following defined:
# 1) -c | --cutoff_file - The file defining cut off values for each gene
# 2) -hmm | --hmm_dir - The directory containing the HMMs for each gene
# 3) -i | --input_dir - The directory with the Peptide sequences produced by the previous step in the pipeline
# 4) -o | --output_dir - The directory to put the output
# 5) -t | --threads - The number of simultaneous hmmsearch runs to run at a time
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
echo hmm_from_pep.sh run on $(date) >> $OUT/log.txt
echo Input Directory = "$INPUT" >> $OUT/log.txt
echo Output Directory = "$OUT" >> $OUT/log.txt
echo Cut off file = "$CUTOFF_FILE" >> $OUT/log.txt
echo HMM Directory = "$HMM_DIR" >> $OUT/log.txt
echo Threads = "$THREADS" >> $OUT/log.txt

function FIND_SPECIES_GENE()
{
        SPLIT_FILE=($(echo $1 | tr "\." "\n" | tr "_" "\n"))
        unset SPLIT_FILE[${#SPLIT_FILE[@]}-1]
        GENE=${SPLIT_FILE[0]}
        SPECIES_SPACED=${SPLIT_FILE[@]:1}
        SPECIES=${SPECIES_SPACED// /_}
        OUT_FILE=$GENE'_'$SPECIES'.out'
}
# Export our function and some variables so they are available to the subshell
export -f FIND_SPECIES_GENE
export THREADS
export OUT
export HMM_DIR
export CUTOFF_FILE
# launch hmmsearch runs, one for each thread
cd $INPUT
ls *.faa | xargs -n 1 -P $THREADS -I % bash -c 'FIND_SPECIES_GENE %; \
HMM=$GENE".hmm"; \
CUTOFF=$(sed -n /$GENE/p $CUTOFF_FILE | sed -e "s,$GENE ,,"); \
hmmsearch --tblout $OUT/$OUT_FILE -T $CUTOFF $HMM_DIR/$HMM % 1> $OUT/stdout.txt 2> $OUT/stderr.txt;'