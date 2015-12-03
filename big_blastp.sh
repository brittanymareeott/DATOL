#!/usr/bin/env bash

# Named variables. Every run needs the following defined:
# 1) -db | --database_dir - The directory containing the blast databases of transcriptome ORFs. Each database should be named in the format species.fa.transdecoder.pep
# 2) -q | --query_dir - The directory containing the fasta files for each gene to use as query terms. Each gene should have a separate fasta file
# 3) -o | --output_dir - The directory with the Peptide sequences produced by the previous step in the pipeline
# 4) -t | --threads - The number of threads to use for the blastp search. Each blastp search will be given 1 thread, but -t will specify how many blastp runs to do at a time.
#
# Example: big_blastp.sh -db /home/mendezg/blast_dbs -q /home/mendezg/cegma/fasta -o /home/mendezg/my_cool_critters/big_blastp -t 20

# Directory of BLAST DBs of transcriptome ORFs:
# $datol/blast_dbs

# Location of original CEGMA fasta files:
# $datol/data/fasta2

#Location of BUSCO fasta files:
#$datol/BUSCO/ancestral_split


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
echo Database Directory = "$DBS" >> $OUT/log.txt
echo Output Directory = "$OUT" >> $OUT/log.txt
echo Query Directory = "$QUERY" >> $OUT/log.txt
echo Threads = "$THREADS" >> $OUT/log.txt


#Loop through all protein blast databases produced from transdecoder (based on them ending in .fa.transdecoder.pep.phr)
cd $QUERY
for i in $DBS/*.fa.transdecoder.pep.phr
  do
    DBI=${i##*/}
    DBSPECIES=${DBI/.phr/}
    printf "***********   Starting blastp for $DBSPECIES on `date` ...\n"
# Using the xargs command to multithread the blastp job. I pipe the output of the ls command on the query directory into blastp
# 5 jobs at a time and give each of those jobs 4 threads to work with, usinga total of 20 threads to run 5 blastp jobs at ounce.
    ls $QUERY | xargs -n 1 -P $THREADS -I % blastp -query % -db ${i/.phr/} -num_threads 1 -outfmt 5 -out $OUT/'query_'%'_db_'$DBSPECIES'.blastp.out'
  done
