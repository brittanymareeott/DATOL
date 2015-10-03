#!/usr/bin/env bash
# Named variables. Every run needs the following defined:
# 1) -i | --input_dir - The directory containing the fasta files that need to be aligned.
# 2) -e | --file_extension - The file extension of the fasta input files.
# 3) -t | --threads - The number of threads to use for each mafft run. NOTE: -t X -r MUST NOT exceed the total number of threads available.
# 4) -r | --runs - The number of simultaneous MAFFT runs to run. NOTE: -t X -r MUST NOT exceed the total number of threads available.
#
# Example: big_blastp.sh -i /home/mendezg/fasta -e .fas -t 4 -r 5

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
  -r|--runs)
  RUNS="$2"
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
mkdir $INPUT/alignments
OUT=$INPUT/alignments
echo Input Directory = "$INPUT" >> $OUT/log.txt
echo Output Directory = "$OUT" >> $OUT/log.txt
echo File Extension = "$EXT" >> $OUT/log.txt
echo Threads = "$THREADS" >> $OUT/log.txt
echo Simultaneous MAFFT runs = "$RUNS" >> $OUT/log.txt
cd $INPUT
find *$EXT | xargs -n 1 -P $RUNS -I % sh -c "\
mafft --thread $THREADS --auto % > $OUT/%\
"
