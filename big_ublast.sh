#!/usr/bin/env bash
#
# big_ublast.sh
#
# Author: Gregg Mendez
#
# Software Dependencies:
# usearch v8.1 : http://drive5.com/usearch/

# Named variables. Every run needs the following defined:
# 1) -db | --database_dir - The directory containing the ublast databases of transcriptome ORFs. Each database should be named in the format species_genus.udb
# 2) -q | --query_dir - The directory containing the fasta files for each gene to use as query terms. Each gene should have a separate fasta file
# 3) -o | --output_dir - The directory with the Peptide sequences produced by the previous step in the pipeline
# 4) -t | --threads - The number of threads to use for the ublast. Each ublast search will be given 6 threads.
#
# Example: big_ublast.sh -db /home/mendezg/blast_dbs -q /home/mendezg/cegma/fasta -o /home/mendezg/my_cool_critters/big_ublast -t 20


#Code to handle the named variable inputs:
while [[ $# > 1 ]]
do
key="$1"

case $key in
  -db|--database_dir)
  DBS="$2"
  shift # past argument
  ;;
  -o|--output_dir)
  OUT="$2"
  shift # past argument
  ;;
  -q|--query_dir)
  QUERY="$2"
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

export DBS
export OUT
export QUERY
#Loop through all protein ublast databases
SPECIES=($(find $DBS -name "*.udb" | sed 's#.*/##' ))

# Divide THREADS by 6 to see how many ublast searches to run at a time
RUNS=$( echo "scale=0;$THREADS/6" | bc -l )
export RUNS
# Using the xargs command to multithread the ublast job.
printf "%s\n" "${SPECIES[@]}" | xargs -n 1 -P 1 -I %x bash -c 'FILE=%x;\
    SPECIES=${FILE/.udb/};\
    QUERIES=($(find $QUERY -type f -exec basename {} \; | sed -e "s,.fas,," )); \
    printf "***********   Starting ublast for $SPECIES on `date` ...\n";\
    printf "%s\n" "${QUERIES[@]}" | xargs -n 1 -P $RUNS -I % usearch -ublast $QUERY/%.fas -db $DBS/$SPECIES".udb" -evalue 1e-9 -threads 6 -userout $OUT/%"_"$SPECIES".txt" -userfields target 1> $OUT"/stdout.out" 2> $OUT"/stderr.out"'
