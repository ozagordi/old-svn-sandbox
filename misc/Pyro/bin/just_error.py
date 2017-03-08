#!/usr/bin/env python

import sys
import os
import logging
import logging.handlers
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
#import pythonlib

LOG_FILENAME = './simple_error.log'
# Make a global logging object.
x = logging.getLogger("logfun")
x.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler(LOG_FILENAME, 'w')
h = logging.handlers.RotatingFileHandler(
    LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
ff = logging.Formatter(
    "%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(ff)
x.addHandler(h)
logfun = logging.getLogger("logfun")

# geometric mixture
geo_freq = [0.5]
for i in range(1,10):
    geo_freq.append(geo_freq[i-1] * 0.5)

def plot_dist(mm, ig):
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    max_count = 50
    
    data_1 = [float(mm.count(i))/len(mm) for i in range(max_count)]
    data_2 = [float(ig.count(i))/len(ig) for i in range(max_count)]
    
    print 'Fraction of reads with 0 mismatch', data_1[0]
    print 'Fraction of reads with 1 mismatches', data_1[1]
    print 'Fraction of reads with 2 mismatches', data_1[2]
    print ''
    print 'Fraction of reads with 0 indels', data_2[0]
    print 'Fraction of reads with 1 indels', data_2[1]
    print 'Fraction of reads with 2 indels', data_2[2]
    print ''
    
    blue = '#2171b5'
    dim_w = 0.35
    fig = plt.figure()
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, \
                        wspace=None, hspace=0.3)
    
    max_ind = 10    
    # the first
    ax1 = fig.add_subplot(111)
    
    to_plot_1 = [data_1[i] for i in range(max_ind)]
    to_plot_1.append(sum(data_1[max_ind:]))
    
    to_plot_2 = [data_2[i] for i in range(max_ind)]
    to_plot_2.append(sum(data_2[max_ind:]))
    
    plt.plot(to_plot_1, 'o-',color=blue, label='mismatches')
    plt.plot(to_plot_2, 's-',color='r', label='internal gaps')
    
    plt.ylabel('fraction of reads')
    ax1.set_xlim(-0.5, max_ind+0.5)
    plt.xticks([0,2,4,6,8,10], ('0', '2', '4', '6', '8', r'$\geq$ 10'))

    plt.ylabel('fraction of reads')
    ax1.legend()
    filename = 'error'
    imtype = 'pdf'
    plt.savefig('%s.%s' % (filename, imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)


def find_closest_here(reads_file):
    '''
    The diff_thresh has been set to 0.025 because even when aligning error-free reads
    to the original haplotypes, the distribution of differences of the best 2
    identities goes from 0.023 to 0.091 (~9%)
    '''
    from pythonlib import Alignment
    import tempfile
    import subprocess
    import heapq
    import operator
    
    # diff_thresh = 0.025
    # abs_thresh = 0.85
    
    ref_file = 'ref.fasta'
    out = tempfile.NamedTemporaryFile()
    outname = out.name
    
    cmline = 'needle -asequence %s -bsequence %s \
              -gapopen 6.0 -gapextend 3.0 -auto -adesshow3 -out %s -aformat3 markx10' \
        % (ref_file, reads_file, outname)
    subprocess.call(cmline, shell=True)
    dd = Alignment.alignfile2dict([outname], 'n', 6.0, 3.0, Verbose = False)
    kh = dd.keys()[0]
    d = Alignment.alignfile2dict([outname], 'n', 6.0, 3.0, Verbose = False)[kh]
    out.close()
    this = []
    mm = []
    ident_2 = []
    ig = []
    for k, v in d.items():
        v.summary()
        ig.append(v.int_gaps)
        this.append(float(v.ident)/(v.stop - v.start + 1))
        mm.append(v.mismatch) #v.stop - v.start + 1 - v.ident
        ident_2.append( float(v.ident)/(v.stop - v.start + 1 - v.int_gaps) )
    return ig, this, mm, ident_2


def reads_error(reads_file):
    '''
    '''
    from Bio import SeqIO
    import pickle
    import gzip
        
    if reads_file.endswith('.gz'):
        rh = gzip.open(reads_file)
    else:
        rh = open(reads_file)
        
    reads_list = [s.seq.tostring().strip('N') for s in SeqIO.parse(rh, 'fasta')]
    tot_bases = sum([len(s) for s in reads_list])
    tot_reads = len(reads_list)
    rh.close()
    print 'On ', reads_file, ' ', tot_reads, ' reads found and ', tot_bases, ' bases'

    if not os.path.isdir('pickles'):
        os.mkdir('pickles')
    pck_file = 'pickles/' + reads_file.strip('./').replace('/', '-') + '.pck'
    
    if os.path.exists(pck_file):
        data = pickle.load(open(pck_file))
        logfun.info('estimates fetched from file %s' % pck_file)
        return data
    
    print >> sys.stderr, 'Reads found', tot_reads
    logfun.info('Reads found %d' % tot_reads)
    int_gaps, identity, mismatches, identity_2 = find_closest_here(reads_file)

    data = int_gaps, identity, mismatches, identity_2, tot_bases
    
    oh = open(pck_file, 'w')
    pickle.dump(data, oh)
    oh.close()
    
    return data


def main():    
    args = sys.argv
    # definition of files, directories and all that
    try:
        data_in = args[1]
    except IndexError:
        sys.exit('usage: %s data' % args[0])
        
    if data_in.endswith('.fa') or data_in.endswith('fasta'):
        int_gaps, identity, mismatches, identity_2, tot_bases = reads_error(data_in)
        print 'tot_bases=', tot_bases
        
        # mm = [r[2] for r in data if r[1] > 0.9]
        # ig = [r[0] for r in data if r[1] > 0.9]
        mm = [a[1] for a in zip(identity, mismatches) if a[0] > 0.9]
        ig = [a[1] for a in zip(identity, int_gaps) if a[0] > 0.9]
        
        print >> sys.stderr, 'Error rate per base (percentage)', 100*float(sum(mm))/tot_bases
        print >> sys.stderr, 'In-dels per base (percentage)', 100*float(sum(ig))/tot_bases
        plot_dist(mm, ig)


if __name__ == '__main__':
    main()
