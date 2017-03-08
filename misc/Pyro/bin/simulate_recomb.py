#!/usr/bin/env python
import sys
import random
import numpy as np
import textwrap
from Bio import SeqIO

args = sys.argv
file_clones = args[1]
rec_rate = 1.0
err_rate = 0.05
n = 1000
bases = set(['A', 'C', 'G', 'T', '', '-'])
clones = [s for s in SeqIO.parse(file_clones, 'fasta')]
s_len = max([len(c) for c in clones])

no_rec = int(n*(1-rec_rate))

oh = open('sim.fasta', 'w')

for i in range(no_rec):
    # flat distribution of reads over clones
    cs = random.choice(clones)
    slist = list(cs.seq.tostring())
    # uniform error rate
    n_err = np.random.binomial(len(slist), err_rate)
    err_pos = random.sample(xrange(len(slist)), n_err)
    for i in err_pos:
        slist[i] = random.choice(list(bases-set(slist[i])))
    mods = ''.join(slist)
    oh.write('>read_from_%s\n' % cs.id.split('|')[0])
    wt = textwrap.fill(mods, 80)
    oh.write(wt + '\n')

rec_eq = 0
for i in range(n-no_rec):
    break_point = int(random.gauss(s_len/2, int(s_len/10)))
    sa, sb = random.sample(clones, 2)
    rec_eq += (sa.id == sb.id)
    stw = sa.seq.tostring()[:break_point] + sb.seq.tostring()[break_point:]
    slist = list(stw)
    n_err = np.random.binomial(len(slist), err_rate)
    err_pos = random.sample(xrange(len(slist)), n_err)
    # modify slist
    for i in err_pos:
        slist[i] = random.choice(list(bases-set(slist[i])))
    mods = ''.join(slist).replace('-', '')
    oh.write('>%s_recomb_%s\n' % (sa.id.split('|')[0], sb.id.split('|')[0]))
    wt = textwrap.fill(mods, 80)
    oh.write(wt + '\n')
assert rec_eq == 0
