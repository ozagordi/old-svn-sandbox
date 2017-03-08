#!/usr/bin/env python

import sys
import os

homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib
from Bio import SeqIO


import time
import optparse
import random, math

ref_genome = os.path.join(homedir, 'References/HIV-HXB2.fasta')

d2i = {
    'A': 0, 'a': 0,
    'C': 1, 'c': 1,
    'G': 2, 'g': 2,
    'T': 3, 't': 3,
    '-': 4,
    'N': 5, 'n': 5
    }

f_code = {}
f_code['A']    = 'A'
f_code['C']    = 'C'
f_code['G']    = 'G'
f_code['T']    = 'T'
f_code['AG']   = 'R'
f_code['CT']   = 'Y'
f_code['GT']   = 'K'
f_code['AC']   = 'M'
f_code['CG']   = 'S'
f_code['AT']   = 'W'
f_code['CGT']  = 'B'
f_code['AGT']  = 'D'
f_code['ACT']  = 'H'
f_code['ACG']  = 'V'
f_code['ACGT'] = 'N'

i2d = ['A', 'C', 'G', 'T', '-', 'N']
B = len(i2d)

thresh = 0.0 #1./(B-1)

def normalize(l):
    "Return a normalized probability vector."
    s = float(sum(l))
    if s == 0:
        return l
    else:
        return [ p/s for p in l ]


def shannon_entropy(l):
    from math import log
    """Return the Shannon entropy of random variable with probability
    vector l."""
    return sum([-p*log(p, 2) for p in l if p > 0])


def far_consensus(in_file, cons_file):
    """
    Parse reads from a file with aligned reads in fasta format
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna
    import random
    
    handle = open(in_file)
    format = 'fasta'
    print >> sys.stderr, 'Parsing aligned reads from file', in_file
    alignment = list(SeqIO.parse(handle, format))
    al_len = 0
    mstart_arr = []
    mod_reads = {}
    gen_length = len(alignment[0].seq)
    consensus = [ [1] * B for i in range(gen_length+1) ]
    n_reads = len(alignment)
    
    for s in alignment:
        id = s.id
        al_len += 1
        
        print >> sys.stderr,  '\x1B[1A\x1B[2K', '%4.1f percent done' % (100. * (al_len+1)/n_reads)
        # print >> sys.stderr, '\x1B[1A\x1B[2K parsing sequence', al_len
        
        assert gen_length == len(s.seq)
        
        if id in mod_reads:
            sys.exit('reads should all have different names')
            
        mod_reads[id] = []
        # sequences go from 1 to gen_length
        mst = s.seq.tostring().rstrip('-')
        
        mls = list(mst)
        mstop = len(mls)
        flanking = True
        mod_reads[id].append('X')
        
        j = 0
        for c in mls:
            j += 1
            
            if c != '-' and flanking:
                mstart = j
                flanking = False
            
            if flanking == True:
                mod_reads[id].append('X')
            else:
                mod_reads[id].append(c)
                base = d2i[c]
                consensus[j][base] += 1
    
        mstart_arr.append(mstart)
    
    # all sequences parsed
        
    
    cseq = []
    j = 1
    
    tmp = []
    for pos in consensus[1:]:
        tmax = max(pos)
        tot_bases = sum(pos)
        
        # this can be optimized
        if pos[4] < 0.8*tot_bases:
            b = []
            for c in range(len(pos)):
                if pos[c] == tmax:
                    b.append(c)
            to_append = i2d[random.choice(b)]
            cseq.append(to_append)
            if to_append != '-':
                tmp.append(pos)
        else:
            pass
#    assert len(cseq) == len(tmp)
    nor_cons = [ normalize(c) for c in tmp ]
    final_cons = ''.join(cseq).replace('-', '')
    
    print 'Length of the consensus', len(final_cons)
    
    cons_str = ''.join(final_cons)
    cons = SeqRecord(Seq(str(cons_str), generic_dna), id=in_file,
                     description='')
    cons_rec = [cons]
    
    ohandle = open(cons_file, "w")
    SeqIO.write(cons_rec, ohandle, "fasta")
    handle.close()
    
    return nor_cons



def plot_entropy(ent_1, ent_2, kld, hotspots, filename):
    """plot
    entropy o both samples and KL distance"""
    
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
        import scipy.stats as stats
    except:
        print 'exit, could not import matplotlib'
        sys.exit()
    """
    w, h = plt.figaspect(2.)
    fig = plt.Figure(figsize=(w,h))
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    """
    ax = plt.subplot(311)
    
    plt.hist(ent_1, 200, normed=1, facecolor = 'r', label='patient 1', alpha=0.5, log=True)
    plt.hist(ent_2, 200, normed=1, facecolor = 'g', label='patient 2', alpha=0.5, log=True)
    
    plt.ylabel('entropy')
    plt.legend(loc='upper right')
    plt.xticks(size='x-small')
    plt.yticks('')
    ax.set_xscale('log')
#    ax.set_xlim(0, 0.15)
    
    ax = plt.subplot(312)
    
    ent_2 = [-e for e in ent_2]
    m1 = stats.scoreatpercentile(ent_1, 98)
    m2 = stats.scoreatpercentile(ent_2, 2)
    
    m1 = np.mean(ent_1)
    m2 = np.mean(ent_2)
    print m1
    print m2
    

    plt.plot(ent_1, 'r', label='patient 1', alpha=0.7)
    plt.plot(ent_2, 'g', label='patient 2', alpha=0.7)
    plt.ylabel('entropy')
    plt.xticks('')
    plt.yticks([-1, 0, 1], ['1', '0', '1'], size='x-small')
    
    plt.axhspan(m2, m1, facecolor='0.5', alpha=0.5)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlim(0, 1300)


    ax = plt.subplot(313)

    try:
        for h in hotspots:
            plt.annotate('*', xy=(h-15, 1.1*kld[h-7]), xycoords='data', textcoords='data', color='magenta')
    except:
        pass
    
    plt.semilogy(kld)
    plt.ylabel('KL divergence')
    plt.xticks(size='x-small')
    plt.yticks([1E-4, 1E-2, 1E0, 1E2], size ='x-small')
    ax.set_xlim(0, 1300)
    
    """
    
    ax = plt.subplot(414)
    diff = [ e[0]+e[1] for e in zip(ent_1, ent_2) ]
    
    plt.plot(diff)
    # plt.xlabel('')
    plt.xlabel('sequence')
    plt.ylabel('$\Delta$ entropy')
    plt.xticks(size='x-small')
    plt.yticks([-1, 0, 1], ['1', '0', '1'], size='x-small')
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlim(0, 1300)


    """
    plt.subplots_adjust(hspace = 0.2)
    # plt.subplots_adjust(wspace = 0.4)
    
    imtype = 'pdf'
    plt.savefig('%s.%s' % (filename, imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()


def align_consensus(cons_file_1, cons_file_2):
    """ Align consensus.faa to each other and to the HIV reference
    """
    from pythonlib import EmbossStandalone
    from pythonlib.MarkxIO import Markx10Iterator
    
    needle_exe = 'needle'
    
    out_file = 'map_cons.needle'
    EmbossStandalone.needle(needle_exe, cons_file_1, cons_file_2,
                            out= out_file, aglobal3='False')
    
    alignment = Markx10Iterator(open(out_file)).next()
    
    cons_1 = alignment.get_seq_by_num(0).tostring().upper()
    cons_2 = alignment.get_seq_by_num(1).tostring().upper()
    
    map = []
    for c in zip(cons_1, cons_2):
        assert not ( c[0] == c[1] and c[1] == '-' ), 'not two gaps'
        if c[0] == '-':
            map.append(1)
        elif c[1] == '-':
            map.append(2)
        else:
            map.append(0)
    
    return map


def KL_distance(d_1, d_2):
    """ Computes Kullback Leibler distance on two distributions
    """
    from math import log
    if sum(d_1) != 1.:
        d_1 = normalize(d_1)
    
    if sum(d_2) != 1.:
        d_2 = normalize(d_2)
    
    assert len(d_1) == len(d_2), 'same number of entries, please'
    
    return sum([c[0]*log(float(c[0])/c[1], 2)  for c in zip(d_1, d_2)])


def main():
    """ What does the main do?
    """
    import pickle
    
    in_file_1 = sys.argv[1]
    in_file_2 = sys.argv[2]
    pick_1 = 'c1.pck'
    pick_2 = 'c2.pck'
    
    cons_1_name = 'consensus_1.faa'
    try:
        h1 = open(pick_1)
        consensus_1 = pickle.load(h1)
        h1.close()
    except:
        consensus_1 = far_consensus(in_file_1, cons_1_name)
        h1 = open(pick_1, 'w')
        pickle.dump(consensus_1, h1)
        h1.close()
    
    tmp_ent_array_1 = [ shannon_entropy(p) for p in consensus_1 ]
    
    cons_2_name = 'consensus_2.faa'
    try:
        h2 = open(pick_2)
        consensus_2 = pickle.load(h2)
        h2.close()
    except:
        consensus_2 = far_consensus(in_file_2, cons_2_name)
        h2 = open(pick_2, 'w')
        pickle.dump(consensus_2, h2)
        h2.close()
    
    tmp_ent_array_2 = [ shannon_entropy(p) for p in consensus_2 ]
    
    al_map = align_consensus('consensus_1.faa', 'consensus_2.faa')
    """
    print len(al_map), len(consensus_1), len(consensus_2)
    
    from Bio import SeqIO
    s1 = list(SeqIO.parse(open('consensus_1.faa'), 'fasta'))[0]
    s2 = list(SeqIO.parse(open('consensus_2.faa'), 'fasta'))[0]
    
    j1 = 0
    j2 = 0
    
    s1_l = list(s1.seq.tostring().replace('-', ''))
    s2_l = list(s2.seq.tostring().replace('-', ''))
    
    for c in al_map:
        if c == 0:
            print s1_l[j1], s2_l[j2]
            j1 += 1
            j2 += 1
        elif c == 1:
            print '-', s2_l[j2]
            #j1 += 1
            j2 += 1
        elif c == 2:
            print  s1_l[j1], '-'
            j1 += 1
           # j2 += 1
           """
    j1 = 0
    j2 = 0
    
    kld = []
    
    ent_array_1 = []
    ent_array_2 = []
    
    for c in al_map:
        if c == 0:
            kl = 0.5 *( KL_distance(consensus_1[j1], consensus_2[j2]) + KL_distance(consensus_2[j2], consensus_1[j1]) )
            kld.append(kl)
            ent_array_1.append(tmp_ent_array_1[j1])
            ent_array_2.append(tmp_ent_array_2[j2])
            j1 += 1
            j2 += 1
        elif c == 1:
            kld.append(0)
            ent_array_1.append(0)
            ent_array_2.append(tmp_ent_array_2[j2])
#            j1 += 1
            j2 += 1
        elif c == 2:
            kld.append(0)
            ent_array_1.append(tmp_ent_array_1[j1])
            ent_array_2.append(0)
            j1 += 1
 #           j2 += 1
    print len(kld)
    j = 1
    hotspots = []
    for d in kld:
        if d >= 6.0:
            hotspots.append(j)
#            if max(consensus_1[j]) < 0.95 or max(consensus_2[j]) < 0.95:
#            print consensus_1[j], consensus_2[j]
        j += 1
    print >> sys.stderr, 'hotspots:\n',  hotspots
    
    plot_entropy(ent_array_1[6:], ent_array_2[6:], kld[6:], hotspots, 'sample')
    import scipy
    import random
    import numpy as np
    ns1 = [ random.random() for i in range(100) ]
    ns2 = [ random.random() for i in range(100) ]
    
    t = scipy.stats.ansari(ns1, ns2)
    print t
    ns1 = np.array(ent_array_1[6:600])
    ns2 = np.array(ent_array_2[6:600])
    w = scipy.stats.wilcoxon(ns1, ns2)
    print w

if __name__ == '__main__':
    main()
