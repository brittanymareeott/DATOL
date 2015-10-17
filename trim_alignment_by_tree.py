#!/usr/bin/env python

# This script will compare an alignment file (in Fasta or Phylip format) to a newick tree file,
# and create a new alignment file containing only species also in the tree file.
#
# This script takes 2 arguments:
# This first argument must be the newick formatted tree file.
# The second argument must be the fasta or phylip alignment.
#
# Example usage:
# trim_alignment_by_tree.py constraint.tre big_alignment.fasta

from Bio import Phylo, SeqIO, Seq
import sys
import re


#Set tree to the first command line argument (0 is the script itself)
tree = Phylo.read(sys.argv[1], "newick")
INPUT = sys.argv[2]
SPECIES_LIST = []
TREE_LIST = []
GENE = sys.argv[2].split(".")[0]
OUTPUT = '%s.%s.treetrimmed.fasta' % (GENE, sys.argv[1])

#Check if input alignment is a fasta or phylip file
with open(INPUT, 'r') as FILE_CHECK:
	for LINE in FILE_CHECK:
		if ">" in LINE:
			FORMAT = str("FASTA")
			break
		else:
			FORMAT = str("PHYLIP")
			break

#FUNCTIONS
# Generate Dictionary of clade names (always name your clades before running this)
def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

# Name all unnamed clades
def tabulate_names(tree):
    names = {}
    for idx, clade in enumerate(tree.find_clades()):
        if clade.name:
            clade.name = clade.name
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names
# Name all the clades in the tree
tabulate_names(tree)
# Fill Dictionary with clade names
names = lookup_by_names(tree)

# Make list of all species in tree.
for TERMINAL in names['0'].get_terminals():
	TREE_LIST.append(TERMINAL.name)
#Make list of all species in alignment
with open(INPUT, 'r') as ALIGNMENT:
	for LINE in ALIGNMENT:
		if FORMAT is str("FASTA"):
			if ">" in LINE:
				SPECIES = re.sub(r'>([0-9A-Za-z_]+)\n', r'\1', LINE)
				SPECIES_LIST.append(SPECIES)
		else:
			if FORMAT is str("PHYLIP"):
				SPECIES = LINE.split()[0]
				SPECIES_LIST.append(SPECIES)
					
# Make list of species in alignment that aren't in tree
MISSING_LIST = []
for ITEM in SPECIES_LIST:
	if ITEM	not in TREE_LIST:
		MISSING_LIST.append(ITEM)
print TREE_LIST
print SPECIES_LIST
print MISSING_LIST
# Write a new alignment file in fasta format that only includes species from the tree file
for PRUNE in MISSING_LIST:
	SPECIES_LIST.remove(PRUNE)
with open(OUTPUT, 'w') as OUT:
	with open(INPUT, 'r') as ALIGNMENT:
		for Record in SeqIO.parse(ALIGNMENT, "fasta"):
			if Record.name not in MISSING_LIST:
				#sys.stdout.write('>%s\n%s\n' % ( Record.name, Record.seq))
				OUT.write('>%s\n%s\n' % ( Record.name, Record.seq))