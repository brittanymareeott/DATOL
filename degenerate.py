#!/usr/bin/env python
#
# Author: Gregory Mendez
#
# This script uses the degenerate_dna python library from https://github.com/carlosp420/degenerate-dna
#
# This script takes an input fasta file and converts the third codon position to 
# degenerate code that could produce the same amino acid for that codon.
#
# Example:
# degenerate.py --dna /home/mendezg/greenalgae/translations/cds/env10972015_1.cds

from Bio import SeqIO
from degenerate_dna import Degenera
import argparse

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script uses the degenerate_dna python library to substitute the third codon position of input DNA sequences with degenerate codes for the same codon.')
parser.add_argument('--dna', required=True, help='The full file path to the desired input DNA sequence in FASTA format.') 
args = parser.parse_args() 

DNA = args.dna
OUTPUT_FILE = ( "%s.deg.fas" % DNA.split("/")[-1].split(".")[0] )
INPUT_HANDLE = open(DNA, "rU")

FASTA_SEQUENCES = SeqIO.parse(INPUT_HANDLE,'fasta')
with open (OUTPUT_FILE, 'a') as FILE:
    for SeqRecord in FASTA_SEQUENCES:
        res = Degenera(dna=SeqRecord.seq, table=1, method='S')
        res.degenerate()
        DEG = res.degenerated
        FILE.write( ">%s\n%s\n" % ( SeqRecord.id, DEG ) )