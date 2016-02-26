#!/usr/bin/env python
#
# prune_tree.py
#
# Author: Gregory Mendez
#
# This script will compare an alignment file (in Fasta or Phylip format) to a newick tree file,
# and create a new tree file containing only species also in the alignment file.
#
# This script takes 2 arguments:
# This first argument must be the newick formatted tree file.
# The second argument must be the fasta or phylip alignment.
#
# Example usage:
# prune_tree.py constraint.tre big_alignment.fasta

from Bio import Phylo
import sys
import re


#Set tree to the first command line argument (0 is the script itself)
tree = Phylo.read(sys.argv[1], "newick")
INPUT = sys.argv[2]
SPECIES_LIST = []
TREE_LIST = []
GENE = sys.argv[2].split(".")[0]
OUTPUT = '%s.%s' % (GENE, sys.argv[1])

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
# Search tree by clade name
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
				if SPECIES in TREE_LIST:
					SPECIES_LIST.append(SPECIES)
		else:
			if FORMAT is str("PHYLIP"):
			    if re.match(r'\w', LINE):
				    SPECIES = LINE.split()[0]
				    if SPECIES in TREE_LIST:
					    SPECIES_LIST.append(SPECIES)
					
# Make list of species in tree that aren't in alignment
MISSING_LIST = []
for ITEM in TREE_LIST:
	if ITEM	not in SPECIES_LIST:
		MISSING_LIST.append(ITEM)

# Remove missing species from tree and write this edited tree out to a new file
for PRUNE in MISSING_LIST:
	tree.prune(PRUNE)
Phylo.write(tree, OUTPUT, 'newick')