#!/usr/bin/env python

#import os
import sys
import re
#import logging
#import logging.handlers

from Bio import SeqIO

args = sys.argv

try:
    reads_file = args[1]
except:
    sys.exit('usage: reads_coordinate.py reads_file')

s_desc = [s.description for s in SeqIO.parse(open(reads_file), 'fasta')]

pm = re.compile('xy=(\d*_\d*)')
xy = [pm.search(d).group(1) for d in s_desc]
coord = [map(int, f.split('_')) for f in xy]
x = [i[0] for i in coord]
y = [i[1] for i in coord]

print >> sys.stderr, 'Coordinates, I have'

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

plt.scatter(x, y, s=2)
plt.xlim(0, 4000)
plt.ylim(0, 4000)
print >> sys.stderr, 'Saving the figure'

imtype = 'pdf'
plt.savefig('coord.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                orientation='landscape', papertype=None, format=imtype,\
                transparent=False)

