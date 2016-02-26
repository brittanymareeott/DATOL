#!/usr/bin/env python
#
# distance_matrix_zscore.py
#
# Author: Gregory S Mendez

# This script analyzes a set of distance matrixes for different genes for genes with significantly different
# distances than the rest of the genes from the same species. It starts by converting distances into ratios 
# of distance divided by median distance in each tree to standardize distances across genes that may have 
# different substitution rates. Then it looks across genes to compare all of the same distance measurements
# for each gene and converts these ratios into modified z-scores. A modified z-score greater than 3.5 is 
# considered an outlier. Each sequence will have a set of modified z-scores from all the distance measures, so
# you set a cut-off for what percentage of those distances can be outliers before the sequence itself is 
# considered an outlier. Anything over 30% is probably an outlier, but you can set the cutoff however you like.
# If an entire library has a lot of contamination some contaminates will probably slip through this test,
# so it is wise to have another contamination test if you suspect large amounts of contaminating sequences.

# This script takes 4 arguments:
# 1) --input | The directory containing the distance matrix files with names formatted so they all start with RAxML_distances.
# 2) --score | The percent of outlier distances required to flag a sequence as an outlier.
# 3) --tree | The directory containing newick formatted tree files used to generate the distance matrixes (used for visualization only)
# 4) --out | The directory where you want output files written.
# 5) --outgroups | A text file defining the outgroups used to root the tree. If you expect the outgroup taxa not always to be present it is a good idea to provide multiple outgroups.

# Example:
# distance_matrix_zscore.py --input ~/data/munchies/gene_trees/distances/ --score 30 --tree ~/data/munchies/gene_trees/ - ~/data/munchies/outliers/

from __future__ import print_function
import sys, argparse, os
import numpy as np
from glob import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy
from time import time
from ete2 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script will convert a RAxML distance matrix to a modified z-score normalized distance matrix')
parser.add_argument('--input', required=True, help='The input RAxML distance matrix file.')
parser.add_argument('--score', required=True, help='The percent of pairwise distance outliers required to mark an enture gene an outlier.')
parser.add_argument('--tree', required=True, help='The directory containing newick formatted tree files used to generate the distance matrixes (used for visualization only)')
parser.add_argument('--out', required=True, help='The directory where you want output files written.')
parser.add_argument('--outgroups', required=True, help='A text document with your outgroups listed. The line should start with the word Outgroup1 followed by a list of all the species in the outgroup with everything separated by spaces. You can specify Outgroup2 and Outgroup3 on other lines as backup outgroups if no species from your outgroup are present.') 
args = parser.parse_args() 

# Functions to calculate Median absolute deviation, and Modified Z-score
def MAD(VALUES, axis=None):
	return np.median(np.absolute(VALUES - np.median(VALUES, axis)), axis)
def MODZ(VALUES, axis=None):
    if MAD(VALUES) == 0:
        return VALUES
    else:
        return (0.6745 * (VALUES - np.median(VALUES))) / MAD(VALUES)

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
    
def MODZ_ALL(DIST_DIR, SCORE, TREE, OUT_DIR, OUTGROUPS):
#######################################
#        Load all distances into a data matrix
#######################################
    OUTGROUP_FILE = open( OUTGROUPS, 'r' )
    t0 = time()
    print('************************** Loading data from distance files.')
    DATA = {}
    for DIST_FILE in glob('%s/RAxML_distances.*' % DIST_DIR):
        GENE = os.path.basename(DIST_FILE).split('.')[1]
        JUST_DATA = []
        with open(DIST_FILE, 'r') as INPUT:
            for LINE in INPUT:
                # Find distance by splitting the line on a space tab space, select last item in that list, and stripping out the line break
                DATUM = float(LINE.split(' \t ')[-1].strip('\n'))
                JUST_DATA.append(DATUM)
            JUST_MEDIAN = np.median(JUST_DATA)
        with open(DIST_FILE, 'r') as INPUT:
            for LINE in INPUT:
                # Find distance by splitting the line on a space tab space, select last item in that list, and stripping out the line break
                DATUM = float(LINE.split(' \t ')[-1].strip('\n'))
                RATIO = ( DATUM / JUST_MEDIAN )
#             # Starting from the innermost parenthetical we are first splitting on a space tab space, selecting the first item
#             # from that split, splitting that on a space, storing that as a string, sorting that list, then joining the two 
#             # terms with a triple underscore, and setting equal to the variable KEY
                KEY = "___".join(sorted(str(LINE.split(' \t ')[0]).split(' ')))
                #Fill Dictionary with the newly found KEY and DATUM
                try:
                    DATA[KEY][GENE] = RATIO
                except KeyError:
                    DATA[KEY] = { GENE:RATIO }
    t1 = time()
    print('Loading took %f seconds' %( t1 - t0 ))

#######################################
#        Generate Modified Z-Scores For Each Distance
#######################################

# Create a new matrix by copying the data and dividing each datum by 0.1 * median - 
# This centers the distribution on 10 for every distance distribution, which will allow 
# normalizing the data (by taking the log) without shifting the distribution into negative
# values which will cause problems for a Z-Score calculation.
    t0 = time()
    print("************************** Calculating Z-Scores")
    DATA_NORMALIZED = copy.deepcopy(DATA)
    for SPECIES___SPECIES, GENES_DATA in DATA_NORMALIZED.iteritems():
        DEMON = 0.1 * np.median(GENES_DATA.values())
        for GENE, DATUM in GENES_DATA.iteritems():
            DATUM_WEIRDED = DATUM / DEMON
            DATUM_NORMALIZED = np.log(DATUM_WEIRDED)
            GENES_DATA[GENE] = DATUM_NORMALIZED
    DATA_MODZ = copy.deepcopy(DATA_NORMALIZED)
    for SPECIES___SPECIES, GENES_DATA in DATA_MODZ.iteritems():
        MEDIAN = np.median(GENES_DATA.values())
        MAD_GENES_DATA = MAD(GENES_DATA.values())
        if MAD_GENES_DATA == 0:
            for GENE, DATUM in GENES_DATA.iteritems():
                if DATUM >= 2.31:
                    GENES_DATA[GENE] = DATUM + 1.2
        else:
            for GENE, DATUM in GENES_DATA.iteritems():
              DATUM_MODZ = (0.6745 * (DATUM - MEDIAN)) / MAD_GENES_DATA
              GENES_DATA[GENE] = DATUM_MODZ
    t1 = time()
    print('Calculating took %f seconds' %( t1 - t0 ))
    
#######################################
#        Collect Modified Z-Scores by Gene
#######################################

# So here I am generating a new data matrix where all the data is organized by gene, rather
# than by species-to-species distances. Each gene will have a list of species, and each species
# will have a list of distances. Lots of nested dictionaries, but it makes sense to me.

    GENE_DATA = {}
    for SPECIES___SPECIES, GENES_DATA in DATA_NORMALIZED.iteritems():
        SPECIES_1 = SPECIES___SPECIES.split('___')[0]
        SPECIES_2 = SPECIES___SPECIES.split('___')[1]
        for GENE, DATUM in GENES_DATA.iteritems():
            try:
                GENE_DATA[GENE][SPECIES_1][SPECIES___SPECIES] = DATUM
            except KeyError:
                try:
                    GENE_DATA[GENE][SPECIES_1] = {}
                    GENE_DATA[GENE][SPECIES_1][SPECIES___SPECIES] = DATUM
                except KeyError:
                    GENE_DATA[GENE] = {}
                    GENE_DATA[GENE][SPECIES_1] = {}
                    GENE_DATA[GENE][SPECIES_1][SPECIES___SPECIES] = DATUM
            try:
                GENE_DATA[GENE][SPECIES_2][SPECIES___SPECIES] = DATUM
            except KeyError:
                try:
                    GENE_DATA[GENE][SPECIES_2] = {}
                    GENE_DATA[GENE][SPECIES_2][SPECIES___SPECIES] = DATUM
                except KeyError:
                    GENE_DATA[GENE] = {}
                    GENE_DATA[GENE][SPECIES_2] = {}
                    GENE_DATA[GENE][SPECIES_2][SPECIES___SPECIES] = DATUM
                    
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
                                   
#######################################
#        COUNT OUTLIERS AND WRITE FILES
#######################################      
    for GENE, VALUES in GENE_DATA.iteritems():
        BADSPECIES = {}
        TOTAL_DIST = len(VALUES.keys())
        for SPECIES, DISTANCES in VALUES.iteritems():
            SPECIES_COUNT = 0
            for SPECIES___SPECIES, DATUM in DISTANCES.iteritems():
                if DATUM >= 3.5:
                    SPECIES_COUNT = SPECIES_COUNT + 1
            PERCENT_OUTLIERS = float( ( float(SPECIES_COUNT) / float(TOTAL_DIST) ) * 100 )
            if PERCENT_OUTLIERS > float(SCORE):
                if SPECIES not in BADSPECIES.keys():
                    BADSPECIES[SPECIES] = PERCENT_OUTLIERS
        with open ( '%s/outlier_taxa.%s.txt' % ( OUT_DIR, GENE ), 'w') as OUT_OUT:
            for TAXON in BADSPECIES.keys():
	            OUT_OUT.write("%s\n" % TAXON)
	            
    #######################################
    #        GENERATE TREES
    #######################################    
        # Root the tree using the outgroup specified in the text file
        # Next check if our outgroup taxa are in the tree and create a new list of just species present.
        TREE_LIST = []
        # Make list of all species in tree.
        T = Tree( "%s/RAxML_result.%s.constrained.tre" % ( TREE, GENE ) )
        for LEAF in T:
            TREE_LIST.append(LEAF.name)
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
                    print("%s: No outgroup taxa present. Rooting at midpoint instead. This may break a monophyletic group." % GENE )
                    R = T.get_midpoint_outgroup()
                    T.set_outgroup(R)            
        # Write a new tree file with the long branches indicated and their clades indicated	
        for CLADE in T.traverse():
	        CLADE.set_style(nstyle)		
        for LEAF in T:
            if LEAF.name in BADSPECIES.keys():
                LEAF.img_style = RED
                LEAF.add_face(TextFace("\t%.2f" % BADSPECIES[LEAF.name]), column=1, position = "branch-right")
        T.render( '%s/%s.tre.pdf' % ( OUT_DIR, GENE ), tree_style=ts)            
    
#######################################
#        Create Histograms
#######################################
#     print("************************** Starting histograms")
#     TOTAL = len(DATA.keys())
#     COUNT = 0
#     t0 = time()
#     for KEY in DATA:
# #        print KEY
#         X = DATA[KEY].values() 
# #        Y = DATA_NORMALIZED[KEY].values()
# #        F = DATA_WEIRD[KEY].values()
# #        plt.axvline(np.median(DATA_NORMALIZED[KEY].values()), color='b', linestyle='dashed', linewidth=2, label='Median')
# #        plt.axvline(np.median(DATA_NORMALIZED[KEY].values())+MAD(DATA_NORMALIZED[KEY].values()), color='r', linestyle='dashed', linewidth=2, label='MAD')
#         plt.axvline(3.5, color='g', linestyle='dashed', linewidth=1, label='3.5')
#         plt.axvline(-3.5, color='g', linestyle='dashed', linewidth=1)
# #        Z = MODZ(DATA_NORMALIZED[KEY].values())
#         Y = DATA_MODZ[KEY].values()
# #        W = MODZ(
# #     # Write the contents of the new dictionary to a file
# # #     with open(OUT_FILE, 'w') as OUT:
# # #         for key, value in MODZ_DATA.items():
# # #           OUT.write('%s\t%s\n' % (key, value))  
# # #    print MODZ_DATA["Alexandrium_margalefi___Ceratium_fusus"]
# # #    X = MODZ_DATA["Alexandrium_margalefi___Ceratium_fusus"].values() 
# # #     Y = DATA["Alexandrium_margalefi___Ceratium_fusus"].values() 
#         BINS = 50
#         plt.hist(X, bins=BINS, color='c', alpha=0.5, label='Raw')
# #        plt.hist(Z, bins=BINS, color='y', alpha=0.5, label='Mod Z-score')
#         plt.hist(Y, bins=BINS, color='r', alpha=0.5, label='MODZ-matrix')
# #        plt.hist(Z, bins=BINS, color='g', alpha=0.5, label='MODZ-fly')
#         plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
#         plt.xlabel("Value")
#         plt.ylabel("Frequency")
#         plt.savefig(str(KEY) +'.pdf')
#         plt.close()
#         tx = time()
#         ELAPSED_TIME = round(tx - t0 , 2)
#         COUNT = ( COUNT + 1 )
#         PROGRESS = round(( float(COUNT) / float(TOTAL) ) * 100 , 2)
#         PROGRESS_INT = int(( float(COUNT) / float(TOTAL) ) * 100 )
#         REMAINING = int(( 100 - PROGRESS ) * ( ELAPSED_TIME / PROGRESS) )
#         sys.stdout.write( '********* %s%% Completed ****** Estimated Time Remaining %s seconds *********\r' % ( PROGRESS_INT, REMAINING ))
#         sys.stdout.flush()
#     t1 = time()
#     print('***************************   Writing histograms took %f seconds ***************************' %( t1 - t0 ))
        
           
#Invoke the function that does all the work    
MODZ_ALL(args.input, args.score, args.tree, args.out, args.outgroups)