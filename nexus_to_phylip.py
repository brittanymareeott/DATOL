#!/usr/bin/env python

# Author: Gregory S Mendez

# This script will convert a nexus alignment to phylip.

# Named variables. Every run needs the following defined:
# 1) --input - The input file in nexus format.
# 2) --out - The filename and file path desired for the output phylip file.

from Bio import AlignIO
import argparse

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script will create a super matrix alignment file from input alignments')
parser.add_argument('--input', required=True, help='The input file in nexus format.')
parser.add_argument('--out', required=True, help='The filename of the output file.')
args = parser.parse_args() 

INPUT = args.input
OUT = args.out
input_handle = open( INPUT, "rU")
output_handle = open( OUT, "w")
 
ALIGNMENT = AlignIO.parse(input_handle, "nexus")
AlignIO.write(ALIGNMENT, output_handle, "phylip-relaxed")
 
output_handle.close()
input_handle.close()