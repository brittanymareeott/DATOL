#!/usr/bin/env python

from Bio.Blast import NCBIXML
from glob import glob
from Bio import SeqIO, Seq
import sys
import argparse
#from pyfaidx import Fasta

def ExtractPeps(BLAST_Results, OutDirectory):
	SPECIES_GENE = {}
  # Loop over the BLAST XML and retrieve the compIDs
	for BLASTout in glob('%s*out' % BLAST_Results):
		# Define species and gene ID
		SPECIES = '_'.join(BLASTout.split('/')[-1].split('.')[1].split('_')[2:])
		GENE = BLASTout.split('query_')[-1].split('.')[0]
		# Create a list of ORFs we want to retrieve from Transdecoder
		# Take the sequence ID from the BLAST XML and put it into a list
		for blast_record in NCBIXML.parse(open(BLASTout)):
			for Alignment in blast_record.alignments:
				ORF_ID = str(Alignment).split()[1]
				try:
					SPECIES_GENE[SPECIES][GENE].add(ORF_ID)
				except KeyError:
					try:
						SPECIES_GENE[SPECIES]
						SPECIES_GENE[SPECIES][GENE] = set()
						SPECIES_GENE[SPECIES][GENE].add(ORF_ID)
					except KeyError:
						SPECIES_GENE[SPECIES] = {}
						SPECIES_GENE[SPECIES][GENE] = set()
						SPECIES_GENE[SPECIES][GENE].add(ORF_ID)
				continue

	for SPECIES, GENE_ORFS in SPECIES_GENE.iteritems():
		for KOG, GENE_ORF in GENE_ORFS.iteritems():
			OutFasta = '%s/%s_%s.txt' % (OutDirectory,KOG, SPECIES)
			with open(OutFasta, 'a') as Out:
				for VALUE in list(GENE_ORF):
					Out.write('%s\n' % VALUE)
				
# pyfaidx version
#  	for SPECIES, GENE_ORFS in SPECIES_GENE.iteritems():
#  		for KOG, GENE_ORF in GENE_ORFS.iteritems():
# 			ORFs = '%s%s%s' % (ORF_Directory, SPECIES, '.fa.transdecoder.pep')
# 			GENES = Fasta(ORFs)
# 			OutFasta = '%s/%s_%s.faa' % (OutDirectory,KOG, SPECIES)
# 			with open(OutFasta, 'a') as Out:
# 				for KEY in list(GENE_ORF):
# 					Out.write('>%s-%s\n%s\n' % (SPECIES, GENES[KEY].name, GENES[KEY]))
					
# biopython version
# 	for SPECIES, GENE_ORFS in SPECIES_GENE.iteritems():
# 		# Where are the ORFs located - transdecoder
# 		ORFs = '%s%s%s' % (ORF_Directory, SPECIES, '.fa.transdecoder.pep')
# 		with open(ORFs, 'rU') as FastaFile:
# 			for record in SeqIO.parse(FastaFile, "fasta"):
# 				for KOG, GENE_ORF in GENE_ORFS.iteritems():
# 					for ORF in GENE_ORF:
# 						if record.name == ORF:
# 							OutFasta = '%s/%s_%s.faa' % (OutDirectory,KOG, SPECIES)
# 							with open(OutFasta, 'a') as Out:
# 								Out.write('>%s-%s\n%s\n' % (SPECIES, record.name, record.seq))

# Argument Parser
parser = argparse.ArgumentParser(description = 'This is a program...')
parser.add_argument('--blast', required=True, help='BLAST results XML') 
parser.add_argument('--outdir', required=True, help='output gets written here') 
args = parser.parse_args() 

ExtractPeps(args.blast, args.outdir)
