#!/usr/bin/env python
#
# Authors: Gregory S Mendez and Bastian Bentlage
#
# This script fetches sequences listed in hmmsearch results files and writes out a plain text
# file. The plain text file is intended for use by another script to write new fasta files
# with the full length sequences using the def-lines listed in the hmmsearch results and the
# large fasta file used to generate the hmmsearch output file.
#
# This script takes 2 arguments:
# 1) --hmm - The directory containing the hmmsearch output files
# 2) --outdir - The directory to write the text files to
#
# Usage: parse_hmm_search.py --hmm ~/GreenAlgae/hmmsearch_out/ --outdir ~/GreenAlgae/parse_hmm_txt/
#
# The script write_cds_pep.sh should be used after this script to write the fasta files using the text
# files generated by this script. 

from glob import glob
from Bio import SeqIO, Seq
import sys, argparse
import os

def ParseHMM(HMMs, OUT_DIR):
	SPECIES_DICT_TOP_HITS = {}
	SPECIES_DICT_OTHER_HITS = {}
	for HMMout in glob('%s/*.out' % HMMs):
		with open(HMMout, 'r') as HMMhits:
			HMM_SCORE = float(0)
			for Line in HMMhits:
				if '#' not in Line:
#					SPECIES = '-'.join(Line.split()[0].split('-')[:-1])
					SPECIES = '_'.join(HMMout.split('/')[-1].split('.')[0].split('_')[1:])
					GENE = Line.split()[2]
					NEW_HMM_SCORE = float(Line.split()[5])
					ORF_ID = Line.split()[0]
					if HMM_SCORE < NEW_HMM_SCORE:
						HMM_SCORE = NEW_HMM_SCORE
						try:
							SPECIES_DICT_TOP_HITS[SPECIES][GENE].add(ORF_ID)
						except KeyError:
							try:
								SPECIES_DICT_TOP_HITS[SPECIES]
								SPECIES_DICT_TOP_HITS[SPECIES][GENE] = set()
								SPECIES_DICT_TOP_HITS[SPECIES][GENE].add(ORF_ID)
							except KeyError:
								SPECIES_DICT_TOP_HITS[SPECIES] = {}
								SPECIES_DICT_TOP_HITS[SPECIES][GENE] = set()
								SPECIES_DICT_TOP_HITS[SPECIES][GENE].add(ORF_ID)
					elif HMM_SCORE * 0.8 < NEW_HMM_SCORE:
						try:
							SPECIES_DICT_OTHER_HITS[SPECIES][GENE].add(ORF_ID)
						except KeyError:
							try:
								SPECIES_DICT_OTHER_HITS[SPECIES]
								SPECIES_DICT_OTHER_HITS[SPECIES][GENE] = set()
								SPECIES_DICT_OTHER_HITS[SPECIES][GENE].add(ORF_ID)
							except KeyError:
								SPECIES_DICT_OTHER_HITS[SPECIES] = {}
								SPECIES_DICT_OTHER_HITS[SPECIES][GENE] = set()
								SPECIES_DICT_OTHER_HITS[SPECIES][GENE].add(ORF_ID)
					else:
						continue

# Output text files of def-lines
# First check if output subdirectories exist
	if not os.path.exists('%s/TopHits' % OUT_DIR):
		os.makedirs('%s/TopHits' % OUT_DIR)
	if not os.path.exists('%s/OtherHits' % OUT_DIR):
		os.makedirs('%s/OtherHits' % OUT_DIR)
# Top Hits
	for SPECIES, GENE_ORFS in SPECIES_DICT_TOP_HITS.iteritems():
		for KOG, GENE_ORF in GENE_ORFS.iteritems():
			OUT_TXT = '%s/TopHits/%s_%s.txt' % (OUT_DIR, KOG, SPECIES)
			with open(OUT_TXT, 'a') as OUT_TOP:
				for VALUE in list(GENE_ORF):
					OUT_TOP.write('%s\n' % VALUE)
# Other Hits
	for SPECIES, GENE_ORFS in SPECIES_DICT_OTHER_HITS.iteritems():
		for KOG, GENE_ORF in GENE_ORFS.iteritems():
			OUT_OTHER_TXT = '%s/OtherHits/%s_%s.txt' % (OUT_DIR, KOG, SPECIES)
			with open(OUT_OTHER_TXT, 'a') as OUT_OTHER:
				for VALUE in list(GENE_ORF):
					OUT_OTHER.write('%s\n' % VALUE)


# Biopython Method
# Top Hits
#	 for SPECIES, GENE_ORFS in SPECIES_DICT_TOP_HITS.iteritems():
#		 with open('%s/%s_TRANSCRIPTS.fasta' % (OUT_DIR, SPECIES), 'w') as OUT_TRANSCRIPTS:
#			 with open('%s/%s_PEP.fasta' % (OUT_DIR, SPECIES), 'w') as OutPeps:
#				 for GENE, ORF in GENE_ORFS.iteritems():
#					 PeptideFile = '%s/%s.fa.transdecoder.pep' % (Transdecoder, SPECIES)
#					 Transcriptome = '%s/%s.fa.transdecoder.cds' % (Transdecoder, SPECIES)
#					 with open(Transcriptome, 'r') as Transcripts:
#						 for Record in SeqIO.parse(Transcripts, "fasta"):
#							 if ORF in Record.name:
#								 OUT_TRANSCRIPTS.write('>%s___%s-%s\n%s\n' % (GENE, SPECIES, Record.name, Record.seq))
#					 with open(PeptideFile, 'r') as Peps:
#						 for Record in SeqIO.parse(Peps, "fasta"):
#							 if ORF in Record.name:
#								 OutPeps.write('>%s___%s\n%s\n' % (GENE, Record.name, Record.seq))
# Other Hits
#	for SPECIES, GENE_ORFS in SPECIES_DICT_OTHER_HITS.iteritems():
#		with open('GENE_Seqs/OtherHits/%s_TRANSCRIPTS.fasta' % (SPECIES), 'w') as OUT_TRANSCRIPTS:
#			with open('GENE_Seqs/OtherHits/%s_PEP.fasta' % (SPECIES), 'w') as OutPeps:
#				for GENE, ORFs in GENE_ORFS.iteritems():
#					for ORF in ORFs:
#						PeptideFile = '../cegma_assembled_outgroups/%s.fa.transdecoder.pep' % (SPECIES)
#						Transcriptome = '../cegma_assembled_outgroups/%s.fa.transdecoder.cds' % (SPECIES)
#						with open(Transcriptome, 'r') as Transcripts:
#							for Record in SeqIO.parse(Transcripts, "fasta"):
#								if ORF in Record.name:
#									OUT_TRANSCRIPTS.write('>%s___%s-%s\n%s\n' % (GENE, SPECIES, Record.name, Record.seq))
#						with open(PeptideFile, 'r') as Peps:
#							for Record in SeqIO.parse(Peps, "fasta"):
#								if ORF in Record.name:
#									OutPeps.write('>%s___%s\n%s\n' % (GENE, Record.name, Record.seq))

# Argument Parser
parser = argparse.ArgumentParser(description = '# This script fetches sequences listed in hmmsearch results files and writes out a plain text file. The plain text file is intended for use by another script to write new fasta files with the full length sequences using the def-lines listed in the hmmsearch results and the large fasta file used to generate the hmmsearch output file.')

parser.add_argument('--hmm', required=True, help='HMMsearch output dir') 
# Only needed for the pure python method. This method is slow.
#parser.add_argument('--TransDecoder', required=True, help='Transdecoder files') 
parser.add_argument('--outdir', required=True, help='output gets written here') 

args = parser.parse_args() 

ParseHMM(args.hmm, args.outdir)
