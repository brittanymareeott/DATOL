#!/usr/bin/env python

# This is a simple little script that makes a histogram from 2 data files given to it.
# Data should be formatted as a list of numbers, one number per line in the file, or
# as a space or tab separated file with a new entry on each line and all the numbers to
# be graphed in the same column. A third argument must be provided that indicates which
# column the data is in. Use 0 if there are no columns, just data.
#
# Usage: histogram.py inputfile_1.txtinputfile_2.txt 0

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_1 = sys.argv[1]
INPUT_2 = sys.argv[2]
COLUMN = int(sys.argv[3])
TABLE=[]

with open (INPUT_1, "rU") as FILE:
  for LINE in FILE:
    TABLE.append(float(LINE.split()[COLUMN]))
X = TABLE[0:]

TABLE=[]
with open (INPUT_2, "rU") as FILE:
  for LINE in FILE:
    TABLE.append(float(LINE.split()[COLUMN]))
Y = TABLE[0:]

BINS = 30

plt.hist(X, BINS, alpha=0.5, label='Data 1')
plt.hist(Y, BINS, alpha=0.5, label='Data 2')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.savefig(str(INPUT_1) +'.pdf')
plt.close()