#!/usr/bin/env python
#
# This script analyses a tree file in newick format and generates a list of taxa with
# unusually long branches based on the median branch length within the tree. It will 
# also output a file with all the branch lengths and a histogram of the branch lengths
# with the median length and the cutoff score indicated with lines.
# 
# The script takes 2 arguments.
# 1) Tree file in newick format
# 2) Branch length cutoff multiplier. I suggest something around 7.
#
# Usage: long_branches.py KOG0023.tre 7

from Bio import Phylo
import sys
import numpy
from numpy import median, absolute
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#Set tree to the first command line argument (0 is the script itself)
tree = Phylo.read(sys.argv[1], "newick")

#Set gene to the gene name of the input file.
GENE = sys.argv[1].split(".")[1]
#Set the output file name to the gene name.
OUTPUT = '%s.txt' % (GENE)
#Set the multiplier to the second argument:
MULTI = float(sys.argv[2])
#tree.root_with_outgroup('Vitrella_brassicaformis_CEGMA', 'Alveolata_sp_CCMP3155_CEGMA')
# Functions

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

# Search tree by clade name
def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

# Median Absolute Deviation Statistic (MAD)
def mad(data, axis=None):
	return median(absolute(data - median(data, axis)), axis)
	
# Cut off score
def cut(data, axis=None):
	return median(data, axis) + median(absolute(data - median(data, axis)), axis)*MULTI

# First name clades, generate list of branches, calculate some statistics
# Run the name unnamed clades function
tabulate_names(tree)
#Look up info by clade name:
names = lookup_by_names(tree)
#Count species in tree:
TOTAL_SPECIES = len(names['0'].get_terminals())

#Print list of clade names and branch lengths and fill an array with the branch lengths
array_branch_lengths = []
with open (OUTPUT, 'w') as all_branch_lengths:
	for clade in tree.find_clades():
		array_branch_lengths.append(clade.branch_length)
		all_branch_lengths.write('%s\t%s\n' % (clade.name, clade.branch_length))
	all_branch_lengths.write('MEDIAN\t%s\n' % (median(array_branch_lengths)))
	all_branch_lengths.write('MAD\t%s\n' % (mad(array_branch_lengths)))
	all_branch_lengths.write('CUT_OFF\t%s\n' % (cut(array_branch_lengths)))

# Generate Histograms of branch lengths with lines for media, and various MAD Cut Offs
X = array_branch_lengths
plt.hist(X, bins=20, color='c', label='Branch lengths')
plt.axvline(median(array_branch_lengths), color='b', linestyle='dashed', linewidth=2, label='Median')
plt.axvline(cut(array_branch_lengths), color='r', linestyle='dashed', linewidth=2, label='Cut off')
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
plt.savefig(str(GENE) +'.pdf')
plt.close()

#Find Species that are the terminal nodes of any long branches
BadSpecies = []
for clade in tree.find_clades():
	if len(names[clade.name].get_terminals()) < TOTAL_SPECIES * .5:
		if clade.branch_length > cut(array_branch_lengths):
			print clade.name
			print clade.branch_length
			print cut(array_branch_lengths)
			capture = clade.get_terminals()
			for species in capture:
				if species.name not in BadSpecies:
					BadSpecies.append(species.name)