#!/usr/bin/env python
import sys
import matplotlib
import matplotlib.pyplot as plt
from Bio import SeqIO

reads_file = sys.argv[1]

data = [len(s) for s in SeqIO.parse(open(reads_file), 'fasta')]

plt.hist(data)

plt.show()
