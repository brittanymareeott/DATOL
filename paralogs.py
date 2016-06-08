#!/usr/bin/env python
#
# paralogs.py
#
# Author: Gregory Mendez
#
# This script takes a gene tree in which a single assembly has had all of its
# orthologs for a given gene added to the tree and identifies paralogs.
# The paralogs identified fall into two categories:
# 1) In-paralogs - genes duplicated after the last speciation event, so that the
# duplicate gene is only found in members of the same species. We mark the
# sequence with the longest branch, suggesting the shortest be kept as the
# homolog.
# 2) Out-paralogs - We can identify some out paralogs simply by seeing if The
# sequence being examined is sister to a sequence that has previously been
# marked with high confidence as an out-paralog by another script (probably The
# distance_matriz_zscore.py script.).

# This script takes 5 arguments:
# 1) --tree | The tree file (in Newick format) to be examined.
# 2) --others | File with sequence IDS of non-top hits.
# 2) --para | The previously identified paralogs file.
# 3) --out | The directory where you want output files written.
# 4) --outgroups | A text file defining the outgroups used to root the tree. If you expect the outgroup taxa not always to be present it is a good idea to provide multiple outgroups.


from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace
import sys, argparse, os

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script takes a gene tree in which a single assembly has had all of its orthologs for a given gene added to the tree and identifies paralogs.')
parser.add_argument('--tree', required=True, help='File with sequence Ids of non-top hits.')
parser.add_argument('--others', required=True, help='The tree file (in Newick format) to be examined.')
parser.add_argument('--para', required=True, help='The previously identified paralogs file.')
parser.add_argument('--out', required=True, help='The directory where you want output files written.')
parser.add_argument('--outgroups', required=True, help='A text document with your outgroups listed. The line should start with the word Outgroup1 followed by a list of all the species in the outgroup with everything separated by spaces. You can specify Outgroup2 and Outgroup3 on other lines as backup outgroups if no species from your outgroup are present.')
args = parser.parse_args()

###############
## Functions
###############

#### ETE TOOLKIT STYLES:
ts = TreeStyle()
ts.scale = 200
RED = NodeStyle()
RED["size"] = 0
RED["vt_line_width"] = 1
RED["hz_line_width"] = 1
RED["vt_line_type"] = 1 # 0 solid, 1 dashed, 2 dotted
RED["hz_line_type"] = 1
RED["bgcolor"] = "#c74d52"
BLUE = NodeStyle()
BLUE["size"] = 0
BLUE["vt_line_width"] = 1
BLUE["hz_line_width"] = 1
BLUE["vt_line_type"] = 1 # 0 solid, 1 dashed, 2 dotted
BLUE["hz_line_type"] = 1
BLUE["bgcolor"] = "#dedede"
YELLOW = NodeStyle()
YELLOW["size"] = 0
YELLOW["vt_line_width"] = 1
YELLOW["hz_line_width"] = 1
YELLOW["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
YELLOW["hz_line_type"] = 1
YELLOW["bgcolor"] = "#7ebcff"

nstyle = NodeStyle()
nstyle["size"] = 0

# Check if sequence of interest is sister to sequence/s that has already been# identified as a paralog.
def OUT_PARA(SEQ_CHECK, SIS):
    # print SIS
    # print PARA_LIST
    if type(SIS) is list:
        SIS_SET = set(SIS)
        if SIS_SET.issubset(PARA_LIST):
            WRITE_LIST.append(SEQ_CHECK)
            # print SEQ_CHECK
    else:
        if SIS in PARA_LIST:
            WRITE_LIST.append(SEQ_CHECK)
            # print('%s is paralog' % SEQ_CHECK)

# Paralog checking loop
def PARA_LOOP():
    if len( IDS ) > 0:
        for ORG in IDS:
            # print('checking %s:' % ORG)
            # First we get the names of the sister leaves
            INT = (T&ORG).get_sisters()
            INT_NAME = INT[0].name
            INT_NAMES=T.search_nodes(name=INT_NAME)[0]
            # Get the node name for the leaf so we can check distances
            ORG_NODE = T.search_nodes(name=ORG)[0]
            # Lets take the easiest case first. A single sister taxon
            if len( INT_NAMES ) == 1:
                for LEAF in INT_NAMES.iter_leaves():
                    # get the species name of the sister taxon
                    SISTER_BINOM = "_".join( [ LEAF.name.split('_')[0], LEAF.name.split('_')[1] ] )
                    SISTER_ASSEMBLY = LEAF.name.split('___')[0]
                    # if its the same then we have an inparalog
                    if BINOM == SISTER_BINOM:
                        # Now check if it is the same assembly
                        if SISTER_ASSEMBLY == ORG_ID:
                            # Mark the longer branch as an inparalog
                            DISTANCES = ( LEAF.dist, ORG_NODE.dist) # put distance in list
                            SORT_DIST=(sorted(DISTANCES)) # sort the list
                            LONGEST = SORT_DIST[1] # mark second in list as longest
                            if LONGEST == LEAF.dist:
                                LONG_NODE = LEAF
                                PRUNE_LIST.remove(LONG_NODE.name) # remove the inparalog from list
                                T.prune( PRUNE_LIST, preserve_branch_length=True ) # remove the long inparalog from the tree
                                # Remove this sequence from our loop since its fully investigated.
                                WRITE_LIST.append(LONG_NODE.name)
                                IDS.remove(LONG_NODE.name) # take the pruned sequence out of the list we are looping through.
                                ##### write long inparalog to text file of inparalogs
                                # print('single sister. long = %s' % LONG_NODE.name)
                                # print('%s removed from IDS' % ORG)
                            else:
                                LONG_NODE = ORG_NODE
                                PRUNE_LIST.remove(LONG_NODE.name) # remove the inparalog from list
                                T.prune( PRUNE_LIST, preserve_branch_length=True ) # remove the long inparalog from the tree
                                # Remove this sequence from our loop since its fully investigated.
                                WRITE_LIST.append(LONG_NODE.name)
                                IDS.remove(LONG_NODE.name) # take the pruned sequence out of the list we are looping through.
                                ##### write long inparalog to text file of inparalogs
                                # print('single sister. long = %s' % LONG_NODE.name)
                                # print('%s removed from IDS' % ORG)
                        else:
                            # So its an inparalog, but the other sequence is from another assembly of the same species. Let's prune the other assembly just so its easier to consider whether this is the best sequence from this assembly to keep.
                            PRUNE_LIST.remove(LEAF.name) # remove the inparalog from list
                            T.prune( PRUNE_LIST, preserve_branch_length=True ) # remove the inparaog from the tree
                            # print('trimming other assembly')
                    else:
                        # its not an inparalog, but lets still check if this sister has already been flagged as an outparalog, because then we want to flag this the same way.

                        OUT_PARA(ORG, LEAF.name)
                        # print('single out paralog. sister: %s' % LEAF.name)
                        # Remove this sequence from our loop since its fully investigated.
                        # print('%s removed from IDS' % ORG)
                        IDS.remove(ORG) # take the pruned sequence out of the list we are looping through.
            # Now we handle what happens when the sister is a clade
            else:
                # first we need to see if any leaves are from this assembly, if not then we're done with this sequence.
                SISTER_LEAVES = []
                for LEAF in INT_NAMES.iter_leaves():
                    SISTER_LEAVES.append(LEAF.name)
                # Are any of the sister species from the same species?
                if any( BINOM in "_".join( [ s.split('_')[0], s.split('_')[1] ] ) for s in SISTER_LEAVES ):
                    # if these are all the same species we have in paralogs, if not then not an inparalog
                    if all( BINOM in "_".join( [ s.split('_')[0], s.split('_')[1] ] ) for s in SISTER_LEAVES ):
                        # if all of these are from different assemblies then trim them all out of the tree so we can reevaluate.
                        if all( ORG_ID != s.split('___')[0] for s in SISTER_LEAVES ):
                            for SISTER in SISTER_LEAVES:
                                PRUNE_LIST.remove(SISTER)
                                T.prune( PRUNE_LIST, preserve_branch_length=True ) # remove the long inparalog from the tree
                            # print('all sisters are from another assmebly')
                    else:
                        # its not an inparalog, but lets still check if this sister has already been flagged as an outparalog, because then we want to flag this the same way.
                        # print('sister clade includes some other taxa.')
                        # print SISTER_LEAVES
                        SUB_SISTERS = []
                        for CLADE_SEQ in SISTER_LEAVES:
                            # print('%s' % CLADE_SEQ)
                            # print ORG_ID
                            if ORG_ID != CLADE_SEQ.split('___')[0]:
                                SUB_SISTERS.append(CLADE_SEQ)
                        IDS.remove(ORG) # take the pruned sequence out of the list we are looping through.
                        # print('removing %s from IDS' % ORG)
                        OUT_PARA(ORG, SUB_SISTERS)
                else:
                    # its not an inparalog, but lets still check if this sister has already been flagged as an outparalog, because then we want to flag this the same way.
                    OUT_PARA(ORG, SISTER_LEAVES)
                    IDS.remove(ORG) # take the pruned sequence out of the list we are looping through.
                    # print('sister clade includes no same taxa')
        PARA_LOOP()

# Print tree to stdout
# print T.get_ascii(show_internal=True)

TREE = args.tree
# Load Newick tree file
T = Tree(TREE)
# Get assembly name
ORG_ID = TREE.split('.')[1]
BINOM = "_".join( [ ORG_ID.split('_')[0], ORG_ID.split('_')[1] ] )
# Get Species list
SPECIES_LIST = []
for LEAF in T:
    SPECIES = LEAF.name
    if SPECIES not in SPECIES_LIST:
        SPECIES_LIST.append(SPECIES)

PRUNE_LIST = list( SPECIES_LIST ) # copy the species list

# Get a list of all the sequence names
IDS = [ S for S in SPECIES_LIST if ORG_ID == S.split('___')[0] ]
OTHERS = [ LINE.strip() for LINE in open(args.others) ]
TOP = [ x for x in IDS if x not in OTHERS ]

# Root the tree
# First loop through text file and save the taxa from the Outgroup line to a list
OUTGROUP1 = []
OUTGROUP2 = []
OUTGROUP3 = []
# Make lists of each outgroup
OUTGROUP_FILE = open( args.outgroups, 'r' )
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
# Root the tree using the outgroup specified in the text file
# Next check if our outgroup taxa are in the tree and create a new list of just species present.
TREE_LIST = {}
# Make list of all species in tree.
for LEAF in T:
    SPECIES = LEAF.name.split("___")[0]
    SEQID = LEAF.name.split("___")[1]
    TREE_LIST[SPECIES] = SEQID
NEW_OUTGROUP = []
for SPECIES, SEQID in TREE_LIST.iteritems():
    if SPECIES in OUTGROUP1:
        FULL_NAME = "___".join( [ SPECIES, SEQID ] )
        NEW_OUTGROUP.append( FULL_NAME )
# Root tree using the Outgroup taxa that are present, and if no outgroup taxa are present use the midpoint method to root the tree.
def ROOT():
    if len( NEW_OUTGROUP ) > 1:
        ANCESTOR = T.get_common_ancestor( NEW_OUTGROUP )
        T.set_outgroup( ANCESTOR )
    if len( NEW_OUTGROUP ) == 1:
        T.set_outgroup( NEW_OUTGROUP[0] )
    if len( NEW_OUTGROUP ) < 1:
        for SPECIES, SEQID in TREE_LIST.iteritems():
            if SPECIES in OUTGROUP2:
                FULL_NAME = "___".join( [ SPECIES, SEQID ] )
                NEW_OUTGROUP.append( FULL_NAME )
        if len( NEW_OUTGROUP ) > 1:
            ANCESTOR = T.get_common_ancestor( NEW_OUTGROUP )
            T.set_outgroup( ANCESTOR )
        if len( NEW_OUTGROUP ) == 1:
            T.set_outgroup( NEW_OUTGROUP[0] )
        if len( NEW_OUTGROUP ) < 1:
            for SPECIES, SEQID in TREE_LIST.iteritems():
                if SPECIES in OUTGROUP3:
                    FULL_NAME = "___".join( [ SPECIES, SEQID ] )
                    NEW_OUTGROUP.append( FULL_NAME )
            if len( NEW_OUTGROUP ) > 1:
                ANCESTOR = T.get_common_ancestor( NEW_OUTGROUP )
                T.set_outgroup( ANCESTOR )
            if len( NEW_OUTGROUP ) == 1:
                T.set_outgroup( NEW_OUTGROUP[0] )
            if len( NEW_OUTGROUP ) < 1:
                print("%s: No outgroup taxa present. Rooting at midpoint instead. This may break a monophyletic group." % TREE )
                R = T.get_midpoint_outgroup()
                T.set_outgroup(R)
ROOT()

# Name unnamed nodes
NUMBER = 0
for NODE in T.traverse():
    if NODE.name == '':
        NODE.name = NUMBER
        NUMBER += 1

# Load the contents of the Paralog files
try:
    FILE_LIST = [ LINE.strip() for LINE in open(args.para) ]
    PARA_LIST = set()
    for ITEM in SPECIES_LIST:
        if ITEM not in IDS:
            ITEM_BINOM = ITEM.split('___')[0]
            if ITEM_BINOM in FILE_LIST:
                PARA_LIST.add(ITEM)
    if BINOM in FILE_LIST:
        PARA_LIST.add(TOP[0])
except IOError:
    PARA_LIST = set()

WRITE_LIST = []

# Then we start our loop. This calls itself so it will continue looping until
# every ortholog has been investigated.
PARA_LOOP()

# Now we write out our paralogs file
# with open ( '%s/%s___paralogs.txt' % ( args.out, ORG_ID ), 'a') as PARALOGS:
#     for WRITE in WRITE_LIST:
#         PARALOGS.write("%s\n" % WRITE.split('___')[1])

# Re-Load Newick tree file
T = Tree(TREE)
ROOT()
IDS = [ S for S in SPECIES_LIST if ORG_ID == S.split('___')[0] ]
# Write a new tree file with the paralogs indicated and their clades indicated
for CLADE in T.traverse():
    CLADE.set_style(nstyle)
for LEAF in T:
    if LEAF.name in IDS:
        LEAF.img_style = YELLOW
    if LEAF.name in WRITE_LIST:
        LEAF.img_style = RED
    if LEAF.name in PARA_LIST:
        LEAF.img_style = BLUE

T.render( '%s/%s.paralog_tree.pdf' % ( args.out, ORG_ID ), tree_style=ts)
