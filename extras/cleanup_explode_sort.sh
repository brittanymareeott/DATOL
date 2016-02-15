#!/usr/bin/env bash

#DEPENDENCIES: This script calls the script from the GMOD Genome Grid project by Don Gilbert gilbertd@indiana.edu. 
# http://iubio.bio.indiana.edu/gmod/genogrid/
# split_multifasta.pl which can be downloaded from
# http://iubio.bio.indiana.edu/gmod/genogrid/scripts/split_multifasta.pl

# Named variables. Every run needs the following defined:
#1) -i | --input_dir - The directory containing the input transcript fasta files
#2) -e | --file_extension - The file extension of the input fasta files
#3) -o | --output_dir - The directory to save the output to.
#4) -p | --gene_prefix - The prefix for all gene names. Examples: KOG, CHLOG, BUSCO, NOG
#5) -l | --gene_name_length - The total length of gene names. Examples: KOG3479 = 7, BUSCO321979 = 11.

#Example: cleanup_explode_sort.sh -i /home/mendezg/DATOL/KOG_Seqs/TopHits -e _TRANSCRIPTS.fasta -o /home/mendezg/DATOL/sorted -p KOG -l 7

#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
  -i|--input_dir)
  INPUT="$2"
  shift # past argument
  ;;
  -e|--file_extension)
  EXT="$2"
  shift # past argument
  ;;
  -o|--output_dir)
  OUT="$2"
  shift # past argument
  ;;
  -p|--gene_prefix)
  PRE="$2"
  shift # past argument
  ;;
  -l|--gene_name_length)
  LENGTH="$2"
  shift # past argument
  ;;
  *)
        # unknown option
  ;;
esac
shift # past argument or value
done

mkdir $OUT
#Write all command line arguments to log file
echo cleanup_explode_sort.sh run on $(date) >> $OUT/log.txt
echo Input Directory = "$INPUT" >> $OUT/log.txt
echo file extension = "$EXT" >> $OUT/log.txt
echo Output Directory = "$OUT" >> $OUT/log.txt
echo Gene Prefix = "$PRE" >> $OUT/log.txt
echo Gene Name Length = "$LENGTH" >> $OUT/log.txt

cd $INPUT
for FILENAME in *$EXT
  do
    DIR=$PWD
    SPECIES=${FILENAME/$EXT/}
    printf "***********   $SPECIES ...\n"
    #Split each species fasta file into individual fasta files for each gene.
    mkdir -p $SPECIES
    split_multifasta.pl -i $FILENAME -o $SPECIES
    cd $SPECIES
    #Loop through the individual gene files and add the species name to the file name with a _1_ in front of each species name to be used later. Also change the file extension from .fsa to .fasta.
    for i in *.fsa
      do mv $i  $(echo $i | sed -e 's,\('$PRE'.*\)*.fsa,\1_1_'$SPECIES'.fsa,' -e 's,.fsa,.fasta,')
    done
    #Loop through the contents of each individual gene files and replace the gene name in the definition line with the species name.
    for j in *.fasta
      do sed -i 's,>'$PRE'.*,>'$SPECIES',' $j;
    done
    # Move all of the individual gene files for this species into the output directory.
    mv *.fasta $OUT
    # Move back to the input directory so we can continue looping through the other species files.
    cd $DIR
    rm -r $SPECIES
  done
# Move to the Output directory  
cd $OUT
# Loop through all fasta files and create a new directory for each gene.
for FILENAME in *.fasta; do
      DIR_NAME=${FILENAME:0:$LENGTH};
      mkdir -p $DIR_NAME
      mv -i $FILENAME $DIR_NAME/${FILENAME}
done
#Loop through all the new gene directories and rename the files and remove the leading gene name and _1_, leaving just the species name. Then change the file extension from .fasta to .fas.
for i in *
  do DIR=$PWD
    if [ -d $i ]
    then cd $i
      for j in *.fasta
        do mv $j $(echo $j | sed -e 's,.*_1_,,' -e 's,.fasta,.fas,') 2>/dev/null
      done
    cd $DIR
    fi
  done
