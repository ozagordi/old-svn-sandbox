#!/usr/bin/env python

import os
import sys
import glob

from Bio import SeqIO

import logging
import logging.handlers

from operator import itemgetter
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


LOG_FILENAME = './compare_samples.log'
# Make a global logging object.
x = logging.getLogger('log_inst')
x.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler('./log', 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter('%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s')
h.setFormatter(f)
x.addHandler(h)


def hapseq(str1, str2):
    '''
    '''
    dist = [1  for a in zip(str1, str2) if a[0] != a[1] and 'N' not in a]
    return sum(dist)

    
def local_hap_dict(hap_file, read_file, history):
    '''
    '''
    
    rd = SeqIO.to_dict(SeqIO.parse(open(read_file), 'fasta'))
    support = history / 100
    n = len(rd)
    threshold = 1.*support/n
    
    d = {}
    for s in SeqIO.parse(open(hap_file), 'fasta'):
        fields = s.id.split('|')
        k, freq = fields[0], float(fields[1])
        if '-' not in s.seq and freq >= threshold:
            d[k] = [freq, s.seq]
    
    to_del = []
    for h in d:
        x.debug('analyzing %s' % h)
        match = 0
        for r in rd:
            match += hapseq(rd[r].seq.tostring(), d[h][1])
        if match < support:
            to_del.append(h)
            x.debug('%s has no sufficient support' % h)


    return d, n



def plot_freq_hap(d1, d2, dir1, dir2, n1, n2):
    '''
    '''
    min_reads = 50

    sd1 = sorted(d1.iteritems(), key=itemgetter(1), reverse=True)
    sd2 = sorted(d2.iteritems(), key=itemgetter(1), reverse=True)
    names1 = [ i[0] for i in sd1 ]
    freq1 = [ i[1][0] for i in sd1 ]

    names2 = [ i[0] for i in sd2 ]
    freq2 = [ i[1][0] for i in sd2 ]
    
    ax = plt.figure()
    
    minmin = 10E-5 #min(freq1 + freq2)/10
    plt.vlines(range(len(freq1)), minmin, freq1, color='cyan', lw=8, label=dir1)
    x2 = [ i+.25 for i in range(len(freq2))]
    plt.vlines(x2, minmin, freq2, color='magenta', lw=8, label=dir2)

    ax2 = plt.gca()
    ax2.set_yscale('log')
    ax2.set_xlim(xmin = -0.75)

    plt.xticks(range(len(freq1)) + x2, names1 + names2)
    plt.setp( ax2.get_xticklabels(), rotation=45, horizontalalignment='right', size='x-small' )
    plt.xlabel('haplotypes')
    plt.ylabel('frequency')

    plt.axhline(1.0*min_reads/n1, color='cyan', label='%d reads on %s' % (min_reads, dir1))
    plt.axhline(1.0*min_reads/n2, color='magenta', label='%d reads on %s' % (min_reads, dir2))
    
    plt.legend()
    
    # plt.title('$\Gamma$, KS test history prehistory gives ' + '%3.2E' % KS[1])
    # ph = plt.axvspan((1.-history)*J, J, facecolor='g', alpha=0.3)
    
    # ann_y = 0.9*(max_g - min_g) + min_g    
    # plt.annotate('history', xy=((1-0.7*history)*J, ann_y))

    # pph = plt.axvspan((1.-2*history)*J, (1.-history)*J, facecolor='r', alpha=0.3)
    fstem = dir1+dir2
    
    imtype = 'pdf'
    plt.savefig('%s-hap_dist.%s' % (fstem, imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)    




def analyze_run(ds_file, history=0.20):
    '''
    '''
    import numpy as np
    from scipy.stats import mstats
    
    h = open(ds_file)
    iterations = (l for l in h if l.startswith('iteration'))
    
    theta = []
    gamma = []
    K = []
    unt = []
    for i in iterations:
        isp = i.split()
        iter, tK, tunt = map(int, isp[1:4])
        tth, tg = map(float, isp[4:6])
        K.append(tK)
        unt.append(tunt)
        theta.append(tth)
        gamma.append(tg)
    K_a = np.array(K)
    unt_a = np.array(unt)
    theta_a = np.array(theta)
    gamma_a = np.array(gamma)
    
    J = len(K_a)
    gamma_history = gamma_a[-history*J:]
    gamma_prehistory = gamma_a[-2*history*J:-history*J]
    
#    hm, hstd = np.mean(gamma_history), np.std(gamma_history)
#    pm, pstd = np.mean(gamma_prehistory), np.std(gamma_prehistory)
#    print 'In history gamma is', hm, '+/-', hstd
#    print 'In prehistory gamma is', pm, '+/-', pstd
    KS = mstats.ks_twosamp(gamma_history, gamma_prehistory)
    
#    try:
#        import matplotlib
#        matplotlib.use('pdf')
#        import matplotlib.pyplot as plt
#        import numpy as np
#    except:
#        print 'exit, could not import matplotlib'
#        sys.exit()
    
    ax = plt.figure()
    plt.plot(gamma_a)
    plt.xlabel('iterations')
    plt.ylabel('$\Gamma$')
    plt.title('$\Gamma$, KS test history prehistory gives ' + '%3.2E' % KS[1])
    ph = plt.axvspan((1.-history)*J, J, facecolor='g', alpha=0.3)
    
    min_g, max_g = min(gamma_a), max(gamma_a)
    ann_y = 0.9*(max_g - min_g) + min_g    
    plt.annotate('history', xy=((1-0.7*history)*J, ann_y))

    # pph = plt.axvspan((1.-2*history)*J, (1.-history)*J, facecolor='r', alpha=0.3)
    fstem = ds_file.split('/')[-2]
    
    imtype = 'pdf'
    plt.savefig('%s-gamma_conv.%s' % (fstem, imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)
    
#    plt.show()

    return J*history


def main():
    '''
    '''
    args = sys.argv
    
    
    x.info( ' '.join(args) )
    
    try:
        dir_0 = args[1]
    except:
        x.exception('error in passing command line arguments')
        sys.exit('usage: compare_samples directory_1')
        
    dir_1, dir_2 = os.listdir(dir_0)
    
    h1 = glob.glob('%s%s/w*haplo.out' % (dir_0, dir_1))
    hap_file_1 = h1[0]
    h2 = glob.glob('%s%s/w*haplo.out' % (dir_0, dir_2))
    hap_file_2 = h2[0]
    
    x.info('hap_file_1 is %s' % hap_file_1)
    x.info('hap_file_2 is %s' % hap_file_2)
    
    assert len(h1) == 1 and len(h2) == 1


    # identify ds-out files
    ds_1 = glob.glob('%s%s/w*ds-out' % (dir_0, dir_1))[0]
    ds_2 = glob.glob('%s%s/w*ds-out' % (dir_0, dir_2))[0]
    
    x.info('ds-out 1 is %s' % ds_1)
    x.info('ds-out 2 is %s' % ds_2)

    # plot sampling procedures, return iterations in history
    history_1 = analyze_run(ds_1)
    history_2 = analyze_run(ds_2)

    # identify read files
    rf_1 = glob.glob('%s%s/w*reads.fas' % (dir_0, dir_1))[0]
    rf_2 = glob.glob('%s%s/w*reads.fas' % (dir_0, dir_2))[0]
    
    x.info('read_file 1 is %s' % rf_1)
    x.info('read_file 2 is %s' % rf_2)
    
    # 
    d1, n1 = local_hap_dict(hap_file_1, rf_1, history_1)
    d2, n2 = local_hap_dict(hap_file_2, rf_2, history_2)
    
    plot_freq_hap(d1, d2, dir_1, dir_2, n1, n2)


if __name__ == '__main__':
    main()
