#!/usr/bin/env python
# modified by OZ

from Bio import SeqIO
import sys
import os
import random


args = sys.argv
del args[0]
#print args
seqs = []

for arg in args:
    if(os.path.isfile(arg)):
        f=open(arg)
        for s in SeqIO.parse(f,'fasta'):
            seqs.append(s)
        f.close()

try:
    size = int(args[1])
    shuffled = random.sample(seqs, size)
except:
    random.shuffle(seqs)
    shuffled = seqs

#    shuffled = random.sample(seqs, len(seqs))

SeqIO.write(shuffled, sys.stdout, 'fasta')
