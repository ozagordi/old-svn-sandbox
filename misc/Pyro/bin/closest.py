#!/usr/bin/env python

import sys
import os
import re
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib

Verbose = False

ref_genomes = ['07-56681',
               '07-56951',
               '08-59712',
               '07-54825',
               '08-04134',
               '08-01315',
               '08-55163',
               '08-57881',
               '08-02659',
               '08-04512']


# geometric mixture
geo_freq = [0.5]
for i in range(1,10):
    geo_freq.append(geo_freq[i-1] * 0.5)


def plot_freq(freq_dict, delta):
    """
    
    """
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
    except:
        print 'exit, could not import matplotlib'
        sys.exit()

    import math
    legend = True
    cb = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow']
    ca = ['#00BFFF', # alternative colour set
          '#BC8F8F',
          '#ADFF2F',
          '#B8860B',
          '#00FF7F',
          '#DC143C',
          '#A0522D',
          '#8A2BE2',
          '#D8BFD8',
          '#87CEEB',
          '#6495ED']
    cn = len(ca)
    
    total = sum(freq_dict.values())
    plt.figure(1, figsize=(12,12))
    plt.subplot(211, axisbg='white')
    
    b = 0.0
    l = 20
    wd = 100
    ci = 0
    
    for c_key in ref_genomes:
        try:
            height = float(freq_dict[c_key])/total
        except KeyError:
            height = 0.0
            """            try:
            h = math.log(height, 0.5)
            except ValueError:
            h = 0.0
            """
        h = height
        tlabel = '%s    %4.2f %%' % (ref_genomes[ci], geo_freq[ci]*100)
        plt.bar(l, h, width=wd, bottom=b, color=ca[ci%cn], alpha=0.3, label=tlabel)
        b += h
        ci += 1


    l = 150

    b = 0
    ci = 0
    for height in geo_freq:
        h = height #math.log(height/minf, 2)
        plt.bar(l, h, width=wd, bottom=b, color=ca[ci%cn], alpha=0.3, label='')#ref_genomes[ci])
        b += h
        ci += 1
    plt.annotate('geometric\ndistribution', xy=(160, 1.02), size=12)

    plt.xlabel('')
    plt.ylabel('frequency')
    plt.axis([0, 600, 0, 1])
    plt.xticks('')    


    plt.legend()

    plt.subplot(212, axisbg='white')
    
    nraw, bins, patches = plt.hist(delta, 20, normed=0, facecolor='green', alpha=0.5, label='raw reads')
    plt.xlabel('best two scores relative delta')
    
    params = {'legend.fontsize': 14, 'font.size': 20}
    plt.rcParams.update(params)
    

    imtype = 'pdf'
    plt.savefig('fig_closest.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()
    
    
    return




def align_info():
    """
    
    """
    from Bio import SeqIO
    from Bio import SeqUtils
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from MarkxIO import Markx10Iterator
    
    f_forward = open('tmp_align_f.needle')
    f_reverse = open('tmp_align_r.needle')
    
    forwardaligniter = Markx10Iterator(f_forward)
    reversealigniter = Markx10Iterator(f_reverse)
    score = {}

    print >> sys.stderr, 'parsing the alignments'
    
    while True:
        try:
            f_align = forwardaligniter.next()
            r_align = reversealigniter.next()
        except:
            break
        if f_align is None or r_align is None:
            break
    
        assert f_align.get_all_seqs()[1].id == r_align.get_all_seqs()[1].id, 'same seq back and forward'
        t_scores = float(f_align._annotations['sw_score']), float(r_align._annotations['sw_score'])
        score[f_align.get_all_seqs()[1].id] = max(t_scores)
        
    return score




def reads2clones_align():
    """
    
    """
    from Bio.Application import generic_run
    from cmline import NeedleCommandline
    go = 6.0
    ge = 3.0
    
    ref_genome = 'tmp.fas'

    needle_run = NeedleCommandline()
    needle_run.set_parameter('-asequence', ref_genome)
    needle_run.set_parameter('-gapopen', go)
    needle_run.set_parameter('-gapextend', ge)
    needle_run.set_parameter('-aformat', 'markx10')
    
    # forward
    needle_run.set_parameter('-bsequence', 'tmp_reads_f.fas')
    outfile = 'tmp_align_f.needle'
    needle_run.set_parameter('-outfile', outfile)

    result_1, messages_1, errors_1 = generic_run(needle_run)
    
    # reverse
    needle_run.set_parameter('-bsequence', 'tmp_reads_r.fas')
    outfile = 'tmp_align_r.needle'
    needle_run.set_parameter('-outfile', outfile)

    result_2, messages_2, errors_2 = generic_run(needle_run)
    
    
    for ar in result_1.available_results():
        print ar, result_1.get_result(ar)
            
    if Verbose:
        for m in messages_1.readlines():
            print >>sys.stderr, m

        for e in errors_1.readlines():
            print >>sys.stderr, e

    for ar in result_2.available_results():
        print ar, result_2.get_result(ar)
            
    if Verbose:
        for m in messages_2.readlines():
            print >>sys.stderr, m

        for e in errors_2.readlines():
            print >>sys.stderr, e
    return

def main():
    
    from Bio import SeqIO
    import pickle
    import operator
    import heapq
    args = sys.argv
    try:
        reads_file, clones_file = args[1].rstrip('/'), args[2]
    except:
        sys.exit('usage: closest.py reads_file clones_file')
    
    reads_dict = {}
    
    f_fasta = open(reads_file)
    seqlist = list(SeqIO.parse(f_fasta, 'fasta'))
    countreads = len(seqlist)

    # forward...
    f_fasta_forward_filename = 'tmp_reads_f.fas'
    f_fasta_forward = open(f_fasta_forward_filename, 'w')
    SeqIO.write(seqlist, f_fasta_forward, 'fasta')
    f_fasta_forward.close()

    # ...and reverse
    for seq in seqlist:
        seq.seq = seq.seq.reverse_complement()
    f_fasta.close()
    f_fasta_reverse_filename = 'tmp_reads_r.fas'
    f_fasta_reverse = open(f_fasta_reverse_filename,'w')
    SeqIO.write(seqlist,f_fasta_reverse,'fasta')
    f_fasta_reverse.close()
    
    # clones
    hc = open(clones_file)
    clones_list = list((SeqIO.parse(hc, 'fasta')))
    
    try:
        wh = open('%s-closest.dict' % reads_file.replace('/', '-'))
        reads_dict = pickle.load(wh)
    except:
        for c in clones_list:
            print >> sys.stderr, c.id, c.seq[:10]
            th = open('tmp.fas', 'w')
            th.write('>%s\n' % c.id)
            th.write('%s' % c.seq.tostring())
            th.close()
            
        # align to single clone
            reads2clones_align()
            t_d = align_info()
            for k in t_d:
                try:
                    reads_dict[k][c.id] = t_d[k]
                except KeyError:
                    reads_dict[k] = {}
                    reads_dict[k][c.id] = t_d[k]

        wh = open('%s-closest.dict' % reads_file.replace('/', '-'), 'w')
        pickle.dump(reads_dict, wh)
        wh.close()
    
    r_keys = reads_dict.keys()
    print >> sys.stderr, 'Lost', len(seqlist) - len(r_keys), 'reads'
    
    clones_freq = {}
    delta = []
    discarded = 0
    for k in reads_dict:
        s = reads_dict[k]
        best2 = heapq.nlargest(2, s.iteritems(), operator.itemgetter(1))
        s1, s2 = best2[0][1], best2[1][1]
        td = (s1-s2)/s1
        delta.append(td)
        
        if td > 0.01:
            best_clone = best2[0][0]
            clones_freq[best_clone] = clones_freq.get(best_clone, 0) + 1
        else:
            discarded += 1
    print clones_freq
    print >> sys.stderr, discarded, 'reads had two matches'
    plot_freq(clones_freq, delta)


if __name__ == '__main__':
    main()
