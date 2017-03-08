#!/usr/bin/env python

import sys
import os
import re
from Bio import SeqIO

homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib

def plot_variation(count):
    """
    """
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
        import scipy.stats as stats
        import math
    except:
        print 'exit, could not import matplotlib'
        sys.exit()
    protease = 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'
    blue = '#2171b5'
    fig = plt.figure(figsize=(12, 3))
    ax = plt.subplot(111)
    plt.subplots_adjust(left=0.06, bottom=0.08, right=0.97, top=0.97,
                        wspace=None, hspace=None)
    ax.set_yscale('log')
    
    xx = (1, 100)
    yy = (1, 1)
    plt.scatter(xx, yy, alpha=.0, label='')
    # sanger limit
    thresh = 0.20
    plt.axhline(thresh, color=blue, alpha=0.9)
    note = '%d' % int(100*thresh) + '%'
    # plt.yticks([thresh], note)
    plt.annotate(note, (0.96*xx[1], thresh*1.1), color=blue)
    
    font = {
        # 'family' : 'Andale Mono',
        # 'family' : 'monospace',
        'family': 'monospace',
        'weight' : 'heavy',
        'size'   : 11
        }
    plt.rc('font', **font)  # pass in the font dict as kwargs
    
    for x, c in enumerate(count):
        if c == {}:
            continue
        col = 'black'
        for k, v in c.items():
            xy = (x, v)
            if v > 0.0002:
                if k != protease[x-1]:
                    col = 'r'
                else:
                    col = 'black'
                plt.text(x, v, k, verticalalignment='center', horizontalalignment='center', color=col, label='', family='monospace')
                
    xt = [1] + range(10, 91, 10) + [99]
    for x in xt:
        plt.axvline(x, color=blue, alpha=0.2)
        
    plt.xlim(-0.5, 100.5)
    plt.ylim(0.001, 1.5)
    
    plt.yticks(weight='bold', size=13, family='sans-serif')
    plt.ylabel('frequency', weight='bold', size=13, family='sans-serif')
    plt.xticks(xt, size=13, weight='bold', family='sans-serif')
    plt.annotate('patient', xy=(0.984, 0.52), xycoords='figure fraction', rotation=90, weight='bold',\
                     family='sans-serif', size=13, horizontalalignment='center', verticalalignment='center')
    
    filename = 'var'
    imtype = 'pdf'
    plt.savefig('%s.%s' % (filename, imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()


def align_codons(h):
    '''
    '''
    hap, freq = h
    import subprocess
    # from pythonlib.Alignment import needle_align
    from pythonlib.MarkxIO import Markx10Iterator
    out_file = 'ppp'
    prot = '~/References/HIV-pol-PR.fasta'
    # -sprotein1 -sprotein2
    comline = 'needle %s asis:%s -outfile %s -gapopen 6.0 -gapextend 3.0 -aformat markx10 -auto' % (prot, hap, out_file)
    retcode = subprocess.call(comline, shell=True)
    assert retcode >= 0, 'Child needle was terminated by signal %d' % -retcode
    alignment = Markx10Iterator(open(out_file)).next()
    os.remove(out_file)
    
    al_prot = alignment.get_seq_by_num(0).tostring().upper()
    al_hap = alignment.get_seq_by_num(1).tostring().upper()
    assert not al_prot.startswith('-'), 'Should not start before protease'
    stop = min(len(al_prot.rstrip('-')), len(al_hap.rstrip('-')))
    
#    if '-' in al_prot.strip('-') or '-' in al_hap.strip('-'):
#        return None, None, None
    
    for i, p in enumerate(zip(al_prot, al_hap)):
        if p[1] != '-':
            start = i+1
            break
    
    to_pass = []        
    for i, p in enumerate(zip(al_prot, al_hap)[start-1:stop]):
        if '-' not in p:# p[1] != '-' and p[0] != '-':
            to_pass.append(p[1])
        elif p[0] == '-':
            pass
 #           print >> sys.stderr, 'Warning!'
        else:
 #           print >> sys.stderr, 'Warning!'
            to_pass.append(p[0])
            
    fh = ''.join(to_pass).rstrip('-')
    
    return start, fh, freq


## def consensus(ac_res):
    
##     from Bio.Seq import translate
    
##     count = [ {} for i in range(102) ]
##     for ar in ac_res:
##         start, residues = ar[:2]
##         if start == None and residues == None:
##             continue
##         if start%3 == 0:
##             read = residues
##         elif start%3 == 1:
##             read = residues[2:]
##         elif start%3 == 2:
##             read = residues[1:]
##         try:
##             aa = translate(read)
##         except:
##             print 'error: read', read
##             continue
        
##         if start%3 == 0:
##             start_a = start/3 + 1
##         if start%3:
##             start_a = start/3 + 2
##         for i, amino in enumerate(aa):
##             count[start_a+i][amino] = count[start_a+i].get(amino, 0) + 1


##     to_pass = []
##     for c in count[2:]:
##         try:
##             best = sorted(c.items(), key=lambda (k,v): (v,k))[0][0]
##         except IndexError:
##             continue
##         to_pass.append(best)

##     return ''.join(to_pass)


def count_codons(haps):

    import pickle
    from Bio.Seq import translate
    from operator import itemgetter
    from pythonlib import Alignment
    from pythonlib import mystats
    
    latex = False # print latex table
    count = [ {} for i in range(102) ]
    oh = open('all.dat', 'w')
    hap_freq = {}
    degeneracy = {}
    mask_mupos = []#[10, 11, 22, 25, 32, 46, 58, 62, 67, 74, 89]
    mupos = []
    # These sequences are HXB2 proteases
    wt_protease = 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF'
    wt_protease_nt= 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTA\
TTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTA\
TAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'
    ac_res = map(align_codons, haps)
    
    protease = wt_protease
    for ar in ac_res:
        start, residues, freq = ar # start here is human (from 1)
        start -= 1 # start here is pythonic (from 0)
        if start == None and residues == None: continue
        
        oh.write('%d %s\n' %( round(freq), wt_protease_nt[:start] + residues + wt_protease_nt[len(residues)+start:]))
        
        if start%3 == 0:
            read = residues
        elif start%3 == 1:
            read = residues[2:]
        elif start%3 == 2:
            read = residues[1:]
        try:
            aa = translate(read) # Biopython
        except:
            print 'error: read', read
            continue
        
        if start%3 == 0:
            start_a = start/3 + 1
        if start%3:
            start_a = start/3 + 2
            
        stop_a = len(aa) + start_a + 1
        
        this_hap = str(protease[:start_a-1] + aa + protease[stop_a-2:])
        
        print this_hap.ljust(100), str(freq).ljust(8) # this is used for resistance prediction, whole haplotype and reads
        for i, c in enumerate(this_hap):
            count[i+1][c] = count[i+1].get(c, 0) + freq
        Alignment.needle_align('asis:%s' % wt_protease, 'asis:%s ' % this_hap, 'tmp', 10.0, 0.5)
        d = Alignment.alignfile2dict(['tmp'], 'n', 10.0, 0.5, Verbose = False)['asis']['asis']
        os.remove('tmp')
        
        mutations = []
        
        for i, c in enumerate(zip(d.seq_a, d.seq_b)):
            pos = i + 1
            if '-' in c:
                continue
            if c[0] != c[1]:
                mutations.append(c[0] + str(pos) + c[1])
                if pos not in mask_mupos: mupos.append(pos)
        signature = ', '.join(mutations)
        hap_freq[signature] = hap_freq.get(signature, 0.0) + freq
        degeneracy[signature] = degeneracy.get(signature, 0) + 1
    print ''
    for k, v in hap_freq.items():
        print str(v).ljust(15), ' ', k
    mupos = sorted(mupos)
    spos = {}
    for i, j in enumerate(mupos):
        spos[j]=i
        
    hf_sorted = sorted(hap_freq.items(), key=itemgetter(1), reverse=True)
    tot_reads = sum([h[1] for h in haps])
    tot_hap = sum(hap_freq.values())
    
    print 'Tot reads after', tot_reads
    print 'Tot', tot_hap
    print 'Simpson\'s index on amino acid sequences = %f +/- %f' % mystats.Simpson(hap_freq.values())
    oh = open('degeneracy.pck', 'w')
    pickle.dump(degeneracy, oh)
    oh.close()
    
    for c in count:
        ts = sum(c.values())
        for k in c.keys():
            c[k] /= ts
    plot_variation(count)
    if not latex:
        return hf_sorted
    print ''
    print '|c'*(1+len(spos))
    for i in mupos:
        print '%s%d & '% ( wt_protease[i-1], i),
    print ''

    return hf_sorted

def test_align_codons():
    from nose.tools import eq_
    
    wt_protease_nt= 'CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT'
    data = wt_protease_nt, 0.5
    res = align_codons(data)
    eq_(res[0], 1)
    eq_(res[1], wt_protease_nt)
    eq_(res[2], 0.5)
    
    data = 'CTCAGGTCACTCTTTGGCA', 123.4
    res = align_codons(data)
    eq_(res[0], 2)
    eq_(res[1], data[0])
    eq_(res[2], 123.4)

    data = 'CTCAGGTCACTCTTGGCA', 123.4 # one deletion
    res = align_codons(data)
    eq_(res[0], 2) # in this version, internal deletions are allowed
    eq_(res[1], 'CTCAGGTCACTCTTTGGCA')

    data = 'CTCAGGTCACTCTTTTGGCA', 123.4 # one insertion
    res = align_codons(data)
    eq_(res[0], 2)
    eq_(res[1], 'CTCAGGTCACTCTTTGGCA')
    eq_(res[2], 123.4)


def test_count_codons():
    from nose.tools import eq_

    haps = [('CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT', 100.0),
            ('CATCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT', 200.0),
            ('GCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT', 301.0)]
    
    res = count_codons(haps)
    eq_(res[0][0], 'P1A')
    eq_(res[0][1], 301.0)
    
    eq_(res[1][0], 'P1H')
    eq_(res[1][1], 200.0)

    eq_(res[2][0], '')
    eq_(res[2][1], 100.0)

    eq_(1, 0, 'test with gaps')
    
def main():
    from pythonlib import mystats
    args = sys.argv
    try:
        sup_file = args[1]
    except KeyError:
        sys.exit('usage: %s sup_file' % args[0])

    if sup_file == 'test':
        test_align_codons()
        test_count_codons()
        sys.exit()

    rule_reads = re.compile('ave_reads=(.*)')
    rule_post = re.compile('posterior=(.*) ave')
    all_sup = SeqIO.parse(open(sup_file), 'fasta')
    post_thresh = 0.95
    haps = []
    tot_raw = 0
    tot_good = []
    for s in all_sup:
        ave_reads = float(re.search(rule_reads, s.description).group(1))
        posterior = float(re.search(rule_post, s.description).group(1))
        tot_raw += ave_reads
        if posterior >= post_thresh:
            haps.append([s.seq.tostring().replace('-', ''), ave_reads])
            tot_good.append(ave_reads)
    print 'Tot reads before found in the support file:', tot_raw
    tg = sum(tot_good)
    print 'Tot reads higher than the %f posterior threshold: %f' % (post_thresh, tg)
    tg_int = map(round, tot_good)
    print 'Simpsons index on nucleotide sequences: %f +/- %f\n' % mystats.Simpson(tg_int)
    
    count_codons(haps)
    
if __name__ == '__main__':
    main()
