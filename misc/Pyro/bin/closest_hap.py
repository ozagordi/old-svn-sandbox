#!/usr/bin/env python

import sys
import os
import re
import logging
import logging.handlers
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib

LOG_FILENAME = './closest_hap.log'
# Make a global logging object.
x = logging.getLogger("logfun")
x.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler(LOG_FILENAME, 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
x.addHandler(h)
logfun = logging.getLogger("logfun")

# geometric mixture
geo_freq = [0.5]
for i in range(1,10):
    geo_freq.append(geo_freq[i-1] * 0.5)

def plot_freq(haps_freq, best, haps_file, id_threshold=2, log_base=2):
    """
    
    """
    import operator
    import pickle
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
    except:
        sys.exit('exit, could not import matplotlib')
    
    import math
    
    merged_freq = {}
    not_merged_freq = {}
    for k, v in best.iteritems():
        if len(v[1]) <= id_threshold:
            if v[1] == []:
                merged_freq[v[0]] = merged_freq.get(v[0], []) + [k]
            else:
                not_merged_freq[k] = [v[0]]

    to_plot = {}
    for k, v in merged_freq.iteritems():
        to_plot[k] = sum([haps_freq[i] for i in v])
        logfun.info('plotting %s with freq %f' % (k, to_plot[k]))
        
    to_plot_mut = {}
    for k, v in not_merged_freq.iteritems():
        to_plot_mut[k] = haps_freq[k]
        logfun.info('plotting mutated %s with freq %f' % (k, to_plot_mut[k]))
    

    d1 = sorted(to_plot.iteritems(), key=operator.itemgetter(1), reverse=True)
    lab1 = [ i[0] for i in d1 ]
    v1 = [math.log(i[1], log_base) for i in d1]
    x1 = np.arange(1, len(v1)+1)

    d2 = sorted(to_plot_mut.iteritems(), key=operator.itemgetter(1), reverse=True)
    lab2 = [best[i[0]][0]+ ' ' + '-'.join(best[i[0]][1]) for i in d2]
    v2 = [math.log(i[1], log_base) for i in d2]
    x2 = np.arange(len(x1)+1, len(x1) + len(v2)+1)
    
    # plot 
    plt.figure(1, figsize=(12, 12))
    
    plt.plot(x1, v1, 'o', alpha=0.6, label='perfectly matching', markersize=12, color='c')
    plt.plot(x2, v2, '^', alpha=0.6, label='mismathcing', markersize=12, color='m')
    #plt.plot(x1, est_to_plot, 'o', label='estimated', markersize=14, alpha=0.7, color='m')
    
    # setup the plot
    plt.title(haps_file)
    plt.axis(xmin=0.4, ymax=0.5, xmax=max(x2)+0.5)
    xt = np.concatenate([x1, x2])
    plt.xticks(xt, lab1+lab2)
    
    plt.setp( plt.gca().get_xticklabels(), rotation=30, horizontalalignment='right', size='x-small')
    plt.grid(alpha = 0.3)
    plt.ylabel(r'frequency [$\log_%d$]' % log_base)
    plt.xlabel('reference')
#    plt.legend(loc='upper right', numpoints=1)
    
    params = {'legend.fontsize': 8, 'font.size': 12}
    plt.rcParams.update(params)
    
    imtype = 'pdf'
    fstem = haps_file.lstrip('../').rstrip('.out').replace('/', '-')
    plt.savefig('%s.%s' % (fstem,imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)
    
    plt.show()
    
    return


def mutations(hap, ref):
    '''
    '''
    from Bio.Application import generic_run
    from pythonlib.cmline import NeedleCommandline
    from pythonlib.MarkxIO import Markx10Iterator
    
    needle_run = NeedleCommandline()
    needle_run.set_parameter('-asequence', 'asis:%s' % hap)
    needle_run.set_parameter('-bsequence', 'asis:%s' % ref)
    needle_run.set_parameter('-aformat', 'markx10')
    needle_run.set_parameter('-gapopen', 10)
    needle_run.set_parameter('-gapextend', 1)
#    needle_run.set_parameter('-usa', True)
#    needle_run.set_parameter('-des', True)
    needle_run.set_parameter('-auto', True)
    needle_run.set_parameter('-outfile', 'tmpout')
    
    result_1, messages_1, errors_1 = generic_run(needle_run)
    for ar in result_1.available_results():
        logfun.info(ar + result_1.get_result(ar))

    for m in messages_1.readlines():
        logfun.debug(m)
    
    for e in errors_1.readlines():
        logfun.debug(e)


    oh = open('tmpout')
    aligniter = Markx10Iterator(oh)
    mt = []
    while True:
        try:
            al = aligniter.next()
        except:
            break
        if al is None:
            break
        
        seq_pair = zip(al.get_all_seqs()[0].seq, al.get_all_seqs()[1].seq)
        i = 0
        for sp in seq_pair:
            i += 1
            if '-' not in sp and 'N' not in sp and sp[0] != sp[1]:
               mt.append('%s%d%s' % (sp[0], i, sp[1])) 
    
    return mt


def align_info(al_file):
    """
    
    """
    from Bio import SeqIO
    from Bio import SeqUtils
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from pythonlib.MarkxIO import Markx10Iterator
    
    f_forward = open(al_file)
    
    forwardaligniter = Markx10Iterator(f_forward)

    score = {}
    
    while True:
        try:
            f_align = forwardaligniter.next()
        except:
            break
        if f_align is None:
            break
        
        f_id = f_align.get_all_seqs()[1].id.split(':')[-1]
        score[f_id] = float(f_align._annotations['sw_score'])

    return score



def haps2clones_align(est_hap):
    """
    
    """
    from Bio.Application import generic_run
    from pythonlib.cmline import NeedleCommandline
    go = 6.0
    ge = 3.0
    Verbose = False
    
    needle_run = NeedleCommandline()
    needle_run.set_parameter('-asequence', est_hap)
    needle_run.set_parameter('-gapopen', go)
    needle_run.set_parameter('-gapextend', ge)
    needle_run.set_parameter('-aformat', 'markx10')
#    needle_run.set_parameter('-usa', True)
#    needle_run.set_parameter('-des', True)
    needle_run.set_parameter('-auto', True)
    needle_run.set_parameter('-bsequence', 'tmp_ref.fas')
    outfile_1 = est_hap + '_hap.needle'
    needle_run.set_parameter('-outfile', outfile_1)
    
    result_1, messages_1, errors_1 = generic_run(needle_run)
    
    for ar in result_1.available_results():
        logfun.info(ar + result_1.get_result(ar))

    for m in messages_1.readlines():
        logfun.debug(m)
    
    for e in errors_1.readlines():
        logfun.debug(e)
            
    return outfile_1


def main():
    '''
    This program takes haplotypes (reconstructed with shorah e.g.) in fasta files where
    the frequency is specified in the fasta header: >hap_name|freq. It aligns them to clones
    in a reference file and plot the relative abundance
    '''
    from Bio import SeqIO
    import pickle
    import operator
    import heapq
    import tempfile
    from multiprocessing import Pool
    from shutil import copy
    
    id_threshold = 2
    args = sys.argv
    logfun.info(' '.join(args))
    
    try:
        haps_file, clones_file = args[1].rstrip('/'), args[2]
    except:
        sys.exit('usage: closest_hap.py haps_file clones_file')
    
        
    f_fasta = open(haps_file)
    haps_dict = SeqIO.to_dict(SeqIO.parse(f_fasta, 'fasta'))
    f_fasta.close()
    counthaps = len(haps_dict)
    
    haps_freq = {}
    for k, v in haps_dict.iteritems():
        hi = k.split('|')[0]
        haps_freq[hi] = float(v.description.split('|')[1])
        
    # normalize
    total = sum(haps_freq.values())
    for k in haps_freq:
        haps_freq[k] /= float(total)
        
    copy(clones_file, 'tmp_ref.fas')
    # clones
    hc = open(clones_file)
    clones_dict = SeqIO.to_dict(SeqIO.parse(hc, 'fasta'))
    hc.close()

    exp_freq = {}
    exp_freq_inv = {}
    for sk in clones_dict:
        s = clones_dict[sk]
        exp_freq[s.id]= float(s.description.split('|')[-1])
        exp_freq_inv[s.description.split('|')[-1]] = s.id
    logfun.info('found expected frequencies in the reference file')
    
    try:
        wh = open('%s-closest.dict' % reads_file.lstrip('../').replace('/', '-'))
        haps_closest_dict = pickle.load(wh)
        logfun.info('dictionary of closest references found, loaded')
    except:
        logfun.info('dictionary of closest references NOT found, aligning all to all')
        refs = []
        c_names = []
        for ck in haps_dict:
            c = haps_dict[ck]
            th = tempfile.NamedTemporaryFile(delete=False)
            th.write('>%s\n' % c.id)
            th.write('%s' % c.seq.tostring())
            th.close()
            refs.append(th.name)
            c_names.append(c.id)
            
        # align to all clones
        logfun.info('now to all references')
        pool = Pool()
        al_results = pool.map(haps2clones_align, refs)
        pool2 = Pool()
        logfun.info('parsing the alignment results')
        t_d = pool2.map(align_info, al_results)
        
        haps_closest_dict = {}
        for cn, rd in zip(c_names, t_d):
            for k in rd:
                kk = exp_freq_inv[k.split('#')[0]]
                try:
                    haps_closest_dict[cn][kk] = rd[k]
                except KeyError:
                    haps_closest_dict[cn] = {}
                    haps_closest_dict[cn][kk] = rd[k]
                    
        wh = open('%s-closest.dict' % haps_file.lstrip('../').replace('/', '-'), 'w')
        pickle.dump(haps_dict, wh)
        wh.close()
        
        # remove references tmp files
        for r in refs:
            os.remove(r)
        # remove alignment tmp files
        for r in al_results:
            os.remove(r)
        try:
            os.remove('tmp_hap.fas')
        except:
            pass
        
    h_keys = haps_closest_dict.keys()
    logfun.info('Lost %d reads' % (len(haps_dict) - len(h_keys)))
    
    best = {}
    to_plot = []
    for k, v in haps_closest_dict.iteritems():
        best2 = heapq.nlargest(2, v.iteritems(), operator.itemgetter(1))
        s1, s2 = best2[0][1], best2[1][1]
        best_clone = best2[0][0]
        
        ksp = k.split('|')
        this_hap = haps_dict[k].seq.tostring().replace('-','').upper()
        this_ref = clones_dict[best_clone].seq.tostring().replace('-','').upper()
        best[ksp[0]] = best_clone, mutations(this_hap, this_ref)
            
    logfun.info('%d haplotypes below the threshold' % len(to_plot))
    plot_freq(haps_freq, best, haps_file)

if __name__ == '__main__':
    main()
