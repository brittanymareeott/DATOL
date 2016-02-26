#!/usr/bin/env python
#
# long_branches.py
#
# Author: Gregory Mendez
#
# This script analyses a tree file in newick format and generates a list of taxa with
# unusually long branches based on the median branch length within the tree. It will 
# also output a file with all the branch lengths and a histogram of the branch lengths
# with the median length and the cutoff score indicated with lines.
# 
# The script takes 4 arguments.
# 1) --tree | Tree file in newick format
# 2) --multi | Branch length cutoff multiplier. I suggest something around 5-10.
# 3) --outgroups | A text document with your outgroups listed. The line should start
# with the word Outgroup1 followed by a list of all the species in the outgroup with
# everything separated by spaces. You can specify Outgroup2 and Outgroup3 on other lines
# as backup outgroups if no species from your outgroup are present.
# 4) --out_dir | The directory in which to save all output files.
#
# Usage: long_branches.py --tree ~/constrained_trees/KOG0023.tre --multi 7 --outgroups ~/clades/outgroups.txt --out_dir ~/long_branches

from Bio import Phylo
import sys, argparse
import numpy
from numpy import median, absolute
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script analyses a tree file in newick format and generates a list of taxa with unusually long branches based on the median branch length within the tree. It will  also output a file with all the branch lengths and a histogram of the branch lengths with the median length and the cutoff score indicated with lines.')

parser.add_argument('--tree', required=True, help='The Newick formatted tree file you want to analyze.') 
parser.add_argument('--multi', required=True, help='The multiplier you want to use to set the cut-off. We recommend 5-10.') 
parser.add_argument('--outgroups', required=True, help='A text document with your outgroups listed. The line should start with the word Outgroup1 followed by a list of all the species in the outgroup with everything separated by spaces. You can specify Outgroup2 and Outgroup3 on other lines as backup outgroups if no species from your outgroup are present.') 
parser.add_argument('--out_dir', required=True, help='The directory in which to save all output files.')
args = parser.parse_args() 

#Set tree to the first command line argument (0 is the script itself)
TREE = Phylo.read(args.tree, "newick")
T = Tree(args.tree)
GENE = args.tree.split(".")[1]
MULTI = float(args.multi)
OUTGROUP_FILE = open( args.outgroups, 'r' )
OUT_DIR = args.out_dir
OUTPUT = '%s/%s.txt' % (OUT_DIR, GENE)

##### Functions
#### PHYLO BASED FUNCTIONS:
# Name all unnamed clades
def tabulate_names(TREE):
    names = {}
    for idx, clade in enumerate(TREE.find_clades()):
        if clade.name:
            clade.name = clade.name
        else:
            clade.name = str(idx)
        names[clade.name] = clade
    return names
# Search tree by clade name
def lookup_by_names(TREE):
    names = {}
    for clade in TREE.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names
#### ETE TOOLKIT BASED FUNCTIONS AND STYLES:
ts = TreeStyle() 
RED = NodeStyle()
RED["size"] = 0
RED["vt_line_width"] = 1
RED["hz_line_width"] = 1
RED["vt_line_type"] = 1 # 0 solid, 1 dashed, 2 dotted
RED["hz_line_type"] = 1
RED["bgcolor"] = "#dedede"
nstyle = NodeStyle()
nstyle["size"] = 0
#### STATISTICS FUNCTIONS:
# Median Absolute Deviation Statistic (MAD)
def mad(data, axis=None):
	return median(absolute(data - median(data, axis)), axis)
# Cut off score
def cut(data, axis=None):
	return median(data, axis) + median(absolute(data - median(data, axis)), axis)*MULTI

##### END FUNCTIONS

# First name clades, generate list of branches, calculate some statistics
# Run the name unnamed clades function
tabulate_names(TREE)
#Look up info by clade name:
names = lookup_by_names(TREE)
#Count species in tree:
TOTAL_SPECIES = len(names['0'].get_terminals())
		
# Root the tree using the outgroup specified in the text file
# First loop through text file and save the taxa from the Outgroup line to a list
OUTGROUP1 = []
OUTGROUP2 = []
OUTGROUP3 = []
# Make lists of each outgroup
for LINE in OUTGROUP_FILE:
	if LINE.split()[0] == "Outgroup1":
		for OUT in LINE.split()[2:]:
			OUTGROUP1.append( OUT )
	if LINE.split()[0] == "Outgroup2":
		for OUT in LINE.split()[2:]:
			OUTGROUP2.append( OUT )
	if LINE.split()[0] == "Outgroup3":
		for OUT in LINE.split()[2:]:
			OUTGROUP3.append( OUT )
# Next check if our outgroup taxa are in the tree and create a new list of just species present.
TREE_LIST = []
# Make list of all species in tree.
for TERMINAL in names['0'].get_terminals():
	TREE_LIST.append(TERMINAL.name)
NEW_OUTGROUP = []
for OUT in OUTGROUP1:
	if OUT in TREE_LIST:
		NEW_OUTGROUP.append( OUT )
# Root tree using the Outgroup taxa that are present, and if no outgroup taxa are present use the midpoint method to root the tree.
if len( NEW_OUTGROUP ) > 1:
	ANCESTOR = T.get_common_ancestor( NEW_OUTGROUP )
	T.set_outgroup( ANCESTOR )
if len( NEW_OUTGROUP ) == 1:
	T.set_outgroup( NEW_OUTGROUP[0] )
if len( NEW_OUTGROUP ) < 1:
	for OUT in OUTGROUP2:
		if OUT in TREE_LIST:
			NEW_OUTGROUP.append( OUT )
	if len( NEW_OUTGROUP ) > 1:
		ANCESTOR = T.get_common_ancestor( NEW_OUTGROUP )
		T.set_outgroup( ANCESTOR )
	if len( NEW_OUTGROUP ) == 1:
		T.set_outgroup( NEW_OUTGROUP[0] )
	if len( NEW_OUTGROUP ) < 1:
		for OUT in OUTGROUP2:
			if OUT in TREE_LIST:
				NEW_OUTGROUP.append( OUT )
		if len( NEW_OUTGROUP ) > 1:
			ANCESTOR = T.get_common_ancestor( NEW_OUTGROUP )
			T.set_outgroup( ANCESTOR )
		if len( NEW_OUTGROUP ) == 1:
			T.set_outgroup( NEW_OUTGROUP[0] )
		if len( NEW_OUTGROUP ) < 1:
			print("%s: No outgroup taxa present. Rooting at midpoint instead. This may break a monophyletic group." % args.tree )
			R = T.get_midpoint_outgroup()
			T.set_outgroup(R)
#Print list of clade names and branch lengths and fill an array with the branch lengths
array_branch_lengths = []
with open (OUTPUT, 'w') as all_branch_lengths:
	for clade in TREE.find_clades():
		array_branch_lengths.append(clade.branch_length)
		all_branch_lengths.write('%s\t%s\n' % (clade.name, clade.branch_length))
	all_branch_lengths.write('MEDIAN\t%s\n' % (median(array_branch_lengths)))
	all_branch_lengths.write('MAD\t%s\n' % (mad(array_branch_lengths)))
	all_branch_lengths.write('CUT_OFF\t%s\n' % (cut(array_branch_lengths)))
ts.scale = 300 / max(array_branch_lengths)
# Generate Histograms of branch lengths with lines for median, and MAD Cut Off
X = array_branch_lengths
plt.hist(X, bins=20, color='c', label='Branch lengths')
plt.axvline(median(array_branch_lengths), color='b', linestyle='dashed', linewidth=2, label='Median')
plt.axvline(cut(array_branch_lengths), color='r', linestyle='dashed', linewidth=2, label='Cut off')
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
plt.savefig('%s/%s.hist.pdf' % (OUT_DIR, GENE))
plt.close()
			
#Find Species that are the terminal nodes of any long branches and write the species to a text document
BadSpecies = []
# for clade in TREE.find_clades():
#     if len(names[clade.name].get_terminals()) < TOTAL_SPECIES * .5:
#         if clade.branch_length > cut(array_branch_lengths):
#             capture = clade.get_terminals()
#             for species in capture:
#                 if species.name not in BadSpecies:
#                     BadSpecies.append(species.name)
# with open ( '%s/longbranch_taxa.%s.txt' % ( OUT_DIR, GENE ), 'w') as LONG_OUT:
# 	for OTU in BadSpecies:
# 		print>>LONG_OUT, OTU
		
# Write a new tree file with the long branches indicated and their clades indicated			
for CLADE in T.traverse():
	CLADE.set_style(nstyle)
	if len(CLADE) < TOTAL_SPECIES * .5:
		if CLADE.dist > cut(array_branch_lengths):
			CLADE.img_style = RED
			CAPTURE = CLADE.get_leaves()
			for SPECIES in CAPTURE:
			    if SPECIES.name not in BadSpecies:
			        BadSpecies.append(SPECIES.name)
T.render( '%s/%s.tre.pdf' % ( OUT_DIR, GENE ), tree_style=ts)
with open ( '%s/longbranch_taxa.%s.txt' % ( OUT_DIR, GENE ), 'w') as LONG_OUT:
	for OTU in BadSpecies:
		print>>LONG_OUT, OTU