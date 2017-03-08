#!/usr/bin/env python

import sys
import os
import re
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
#import pythonlib
from pythonlib import Alignment
from Bio import SeqIO
Verbose = False


def find_best_split(seq):
    '''
    '''
    import heapq
    import operator
    l = len(seq)
    low_lim = int(0.25 * l)
    up_lim = int(0.76 * l)
    step = int(0.25 * l)
    ref_genome = 'all_clones.fas'
    best_score = 0
    best_split = 0
    best_gaps = 0
    
    for split in range(low_lim, up_lim, step):
        s1 = seq.seq[:low_lim]
        s2 = seq.seq[low_lim:]
        h1 = open('tmp1.fas', 'w')
        h1.write('>%s_1\n' % seq.id.split('#')[0])
        h1.write(s1.tostring() + '\n')
        h1.close()
        
        h2 = open('tmp2.fas', 'w')
        h2.write('>%s_1\n' % seq.id.split('#')[0])
        h2.write(s2.tostring() + '\n')
        h2.close()
        
        Alignment.needle_align('tmp1.fas', ref_genome, 'tmp1.needle')
        alset_1 = Alignment.alignfile2set(['tmp1.needle'], 'split_1', 6.0, 3.0)
        os.unlink('tmp1.needle')
        
        Alignment.needle_align('tmp2.fas', ref_genome, 'tmp2.needle')
        alset_2 = Alignment.alignfile2set(['tmp2.needle'], 'split_2', 6.0, 3.0)
        os.unlink('tmp2.needle')
        
        k1 = alset_1.keys()[0]
        l1 =  [(s[0], s[1].score) for s in alset_1[k1].iteritems()]
        best_1 = heapq.nlargest(2, iter(l1), operator.itemgetter(1))

        k2 = alset_2.keys()[0]
        l2 =  [(s[0], s[1].score) for s in alset_2[k2].iteritems()]
        best_2 = heapq.nlargest(2, iter(l2), operator.itemgetter(1))
        
        if best_1[0][1] + best_2[0][1] >= best_score:
            best_score = best_1[0][1] + best_2[0][1]
            best_split = split
            clones = best_1[0][0], best_2[0][0]
            alset_1[k1][clones[0]].summary()
            alset_2[k2][clones[1]].summary()
            best_gaps = alset_1[k1][clones[0]].int_gaps + alset_2[k2][clones[1]].int_gaps
            al_start_1, al_stop_1 = alset_1[k1][clones[0]].start, alset_1[k1][clones[0]].stop
            al_start_2, al_stop_2 = alset_2[k2][clones[1]].start, alset_2[k2][clones[1]].stop

    del alset_1
    del alset_2
    return best_score, best_split, best_gaps, clones, al_start_1, al_stop_1, al_start_2, al_stop_2


def main():
    '''
    '''
    
    h = open('outliers.fas')
    seq_dict = SeqIO.to_dict(SeqIO.parse(h, 'fasta'))
    count = 0
    print len(seq_dict), 'outliers found'
    for seq in seq_dict:
        res = find_best_split(seq_dict[seq])
        if res[2] <= 5 and abs(res[6] - res[5]) < 20:
            print >> sys.stderr, seq_dict[seq].id, res[2], res[3], res[4], res[5], res[6], res[7]
            count += 1
    print count, 'matched the requirements'
if __name__ == '__main__':
    main()
