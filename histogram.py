#!/usr/bin/env python

# This is a simple little script that makes a histogram from data files given to it.
# Data should be for matted as a list of numbers, one number per line in the file.
# Usage: histogram.py inputfile.txt

import sys
import matplotlib.pyplot as plt

INPUT=sys.argv[1]
TABLE=[]
with open (INPUT, "rU") as FILE:
  for LINE in FILE:
    TABLE.append(float(LINE.rstrip('\n')))
X = TABLE[0:]

plt.hist(X, bins=20)
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.savefig(str(INPUT) +'.pdf')
plt.close()
