#!/usr/bin/env python

# Author: Gregory S Mendez

# This script will create a super matrix alignment file in nexus format from input alignments in nexus format

# Named variables. Every run needs the following defined:
# 1) --in_dir - The directory containing the nexus alignments that need to be merged.
# 2) --out - The full filepath and name you want for the output file.

from Bio.Nexus import Nexus
import argparse, glob

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script will create a super matrix alignment file from input alignments')
parser.add_argument('--in_dir', required=True, help='The input directory containing alignment files.')
parser.add_argument('--out', required=True, help='The filepath and filename of the output file.')
args = parser.parse_args() 

IN_DIR = args.in_dir
OUT = args.out
FILE_LIST = glob.glob('%s/*.nex' % IN_DIR)
NEXI =  [(FNAME, Nexus.Nexus(FNAME)) for FNAME in FILE_LIST]
COMBINED = Nexus.combine(NEXI)
COMBINED.write_nexus_data(filename=open('%s' % OUT, 'w'))