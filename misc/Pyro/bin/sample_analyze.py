#!/usr/bin/env python

import os
import sys
import logging
import logging.handlers
import socket
from Bio import SeqIO

start, stop = 50, 240
win_size, win_shifts = stop-start, 1

iterations, alpha, m = 20000, 0.05, 10
our_hosts = ['bs-mbp08', 'bs-dsvr07', 'bs-dsvr24']
hostname = socket.gethostname().split('.')[0]


LOG_FILENAME = './sample_analyze.log'
# Make a global logging object.
log_inst = logging.getLogger('log_inst')
log_inst.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler('./log', 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter('%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s')
h.setFormatter(f)
log_inst.addHandler(h)

def n_dist(str1, str2):
    '''
    '''
    a = zip(str1, str2)
    s = [ 1 for i in a if i[0] != i[1] and 'N' not in i and '-' not in i]
    return sum(s)

def haplo_from_reads(cor_reads):
    '''
    '''
    
    from Bio import SeqIO
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from operator import itemgetter
    import cPickle

    suffix = cor_reads.split('.')[-1]
    pck_file = cor_reads.rstrip(suffix) + 'pck'
    
    try:
        hap_tmp = cPickle.load(open(pck_file))
        log_inst.info('found pickle file %s, reading from it' % pck_file)
        log_inst.info('different haplotype=%d, considered reads=%d' % (len(hap_tmp), sum(hap_tmp.values())))
    except:
        thresh = 0.999
        hap_tmp = {}
        discarded = []
        already = []
        i = 0
        count = 0
        for r in SeqIO.parse(open(cor_reads), 'fasta'):
            print >> sys.stderr,  '\x1B[1A\x1B[2K', '%d' % i, '%d' % len(already)
            i += 1
            seq_str =  r.seq.tostring()
            if float(r.description.split('|')[1]) > thresh and '-' not in seq_str and seq_str.count('N') < 5:
                count += 1
                already = hap_tmp.keys()
                added = False
                if already == []:
                    hap_tmp[seq_str] = 1
                for ts in already:
                    if n_dist(seq_str, ts) == 0:
                        hap_tmp[ts] += 1
                        added = True
                        break
                if not added:
                    hap_tmp[seq_str] = 1
            else:
                discarded.append(seq_str)

        out = open(pck_file, 'w')
        cPickle.dump(hap_tmp, out)
        out.close()
        sum1 = sum(hap_tmp.values())
        assert  sum1 == count, 'Lost in conversion'
        log_inst.info('pickle file not found, computed')
        log_inst.info('different haplotype=%d, considered reads=%d, discarded=%d' % (len(hap_tmp), sum1, len(discarded)))
    
    hap_tmp = sorted(hap_tmp.iteritems(), key=itemgetter(1), reverse=True)
    
    # check all haplotypes differ at most because of missing data (base N)
    keys = [ h[0] for h in hap_tmp ]
    lk = len(keys)
    for i in range(lk):
        for j in range(i):
            a, b = keys[i], keys[j]
            assert n_dist(a, b) != 0, '%s and %s are equal' % (a, b)
    

    p = [i[1] for i in hap_tmp]

    l_counts = [50, 20, 10]
    colors = ['r', 'g', 'b']
    sum1 = sum(p)
    freq = [1.*i/sum1 for i in p]
    ax = plt.figure()
    freq_2 = freq[:30]
    minmin = 1.1E-5
    plt.vlines(range(len(freq_2)), minmin, freq_2, color='cyan', lw=8, label='haplotypes')
    for rc  in  zip(l_counts, colors):
        plt.axhline(1.*rc[0]/sum1, label='%d reads' % rc[0], color = rc[1])
    ax2 = plt.gca()
    ax2.set_yscale('log')
    ax2.set_xlim(xmin = -0.75)
    plt.legend()
    imtype = 'pdf'
    plt.savefig('pop.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)
    return hap_tmp

def check_original_reads(hap, orig_reads):
    '''
    '''
    from operator import itemgetter
    support = {}
    i = 1
    l = len(hap)
    orig_list = list(SeqIO.parse(open(orig_reads), 'fasta'))
    
    for h in hap:
        this_hap = h[0]
        print >> sys.stderr, i, 'out of', l
        i += 1
        if h[1] > 0:
            for s in orig_list:
                if n_dist(s.seq.tostring(), this_hap) < 1:
                    support[this_hap] = support.get(this_hap, 0) + 1


    support = sorted(support.iteritems(), key=itemgetter(1), reverse=True)
    
    i = 0
    for k in support:
        print '>hap_%d|%d' % (i, k[1])
        print k[0]
        i += 1

def main():
    '''
    '''
    
    args = sys.argv
    log_inst.info(' '.join(args))
    cor_reads = args[1]
    orig_reads = args[2]
    hap = haplo_from_reads(cor_reads)
    check_original_reads(hap, orig_reads)

if __name__ == '__main__':
    main()
