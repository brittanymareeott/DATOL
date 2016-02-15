#!/usr/bin/env python

# Author: Gregory S Mendez

# This script converts a RAxML distance matrix file into a modified z-score normalized distance matrix.
# Example Usages: mod_zscore.py --input RAxML_distances.KOG0003.txt

import sys, argparse
import numpy as np

# Argument Parser
parser = argparse.ArgumentParser(description = 'This script will convert a RAxML distance matrix to a modified z-score normalized distance matrix')
parser.add_argument('--input', required=True, help='The input RAxML distance matrix file.')
args = parser.parse_args() 

# Functions to calculate Median absolute deviation, and Modified Z-score
def MAD(VALUES, axis=None):
	return np.median(np.absolute(VALUES - np.median(VALUES, axis)), axis)
def MODZ(VALUES, axis=None):
    return (0.6745 * np.absolute(VALUES - np.median(VALUES))) / MAD(VALUES)


def FIND_MODZ(FILE):
    OUT_FILE = '%s.distances.txt' % (FILE.split('.')[1])
    # First load the file into a dictionary
    with open(FILE, 'r') as INPUT:
        DATA = {}
        for LINE in INPUT:
            # Find distance by splitting the line on a space tab space, select last item in that list, and stripping out the line break
            DATUM = float(LINE.split(' \t ')[-1].strip('\n'))
            # Starting from the innermost parenthetical we are first splitting on a space tab space, selecting the first item
            # from that split, splitting that on a space, storing that as a string, sorting that list, then joining the two 
            # terms with a tripple underscore, and setting equal to the variable KEY
            KEY = "___".join(sorted(str(LINE.split(' \t ')[0]).split(' ')))
            #Fill Dictionary with the newly found KEY and DATUM
            DATA[KEY] = DATUM
    # Create a new dictionary by merging the key values from the original dictionary and the newly created modified z-scores
    MODZ_DATA = dict(zip(DATA.keys(), MODZ(DATA.values())))
    # Write the contents of the new dictionary to a file
    with open(OUT_FILE, 'w') as OUT:
        for key, value in MODZ_DATA.items():
          OUT.write('%s\t%s\n' % (key, value))  
          
#invoke the function that does all the work    
FIND_MODZ(args.input)