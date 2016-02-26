#!/usr/bin/env bash
#
# dir_nexus_to_phylip.sh
#
# Author: Gregory Mendez
#
# Scripts called by this script:
# 1) nexus_to_phylip.py
#
# This script converts a directory full of nexus files to phylip.
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
  -o|--output_dir)
  OUT="$2"
  shift # past argument
  ;;
  -e|--extension)
  EXT="$2"
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

# get basename of nexus files and load them into an array
BASE=($(find $INPUT/*.$EXT -type f -exec basename {} \; | sed 's,.nex,,'))

# find all the folders and launch a separate bash shell for each. In each shell concatenate the sequences then align them
printf "%s\n" "${BASE[@]}" | xargs -n 1 -P $THREADS -I % nexus_to_phylip.py --input "$INPUT/"%".$EXT" --out "$OUT/"%".phy"