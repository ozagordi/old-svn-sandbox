#!/usr/bin/env python
'''This takes a far file in input, and removes
the reads that have more than a given percentage
of internal gaps (default=10%) or a stretch of
consecutive gaps longer than 4 nucleotides'''

import sys
from Bio import SeqIO
args = sys.argv
try:
    threshold = args[2]
except IndexError:
    threshold = 0.1
assert threshold < 1.0
inseqs = list(SeqIO.parse(open(args[1]), 'fasta'))
print >> sys.stderr, 'Found', len(inseqs), 'reads'
goodseqs = (s for s in inseqs if float(s.seq.tostring().count('-'))/len(s) < threshold)
goods = SeqIO.write(goodseqs, args[1], 'fasta')
print >> sys.stderr, 'Saving', goods, 'reads, in percent' , (100.0*int(goods))/len(inseqs)
