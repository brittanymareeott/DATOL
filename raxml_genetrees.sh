#!/usr/bin/env bash

# Named variables. Every run needs the following defined:
# 1) -e | --file_extension - The file ending that we should match to use as input alignment files to use for RAXML
# 2) -i | --input_dir - The directory containing the alignments for input into RAXML
# 3) -o | --output_dir - The directory to save the output from RAXML. We will save all files to a directory called genetrees within this output directory.
# 4) -t | --threads - The number of threads to use for each RAXML run.
#
# Example: raxml_genetrees.sh -e .aligned.trim.faa -i /home/mendezg/alignments/fasta -o /home/mendezg/my_cool_critters/ -t 8

#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
  -e|--file_extension)
  EXT="$2"
  shift # past argument
  ;;
  -i|--input_dir)
  INPUT="$2"
  shift # past argument
  ;;
  -o|--output_dir)
  OUT="$2"
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

echo Gene Tree run started on $(date) >> $OUT/log.txt
echo File Extension = "$EXT" >> $OUT/log.txt
echo Output Directory = "$OUT" >> $OUT/log.txt
echo INPUT Directory = "$INPUT" >> $OUT/log.txt
echo Threads = "$THREADS" >> $OUT/log.txt

mkdir $OUT/genetrees
for ALIGNMENT in $INPUT/*$EXT
do
 printf "****************************** RAxML for $ALIGNMENT starting now!"
 raxmlHPC-PTHREADS-SSE3 -f a -p 12345 -x 12345 -# 100 -m GTRGAMMA -s $ALIGNMENT -T $THREADS -n $ALIGNMENT.tre -w $OUT/genetrees
 done
