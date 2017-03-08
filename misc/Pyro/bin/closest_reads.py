#!/usr/bin/env python

import sys
import os
import re
import logging
import logging.handlers
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib

LOG_FILENAME = './closest_reads.log'
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

'''
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
'''

# geometric mixture
geo_freq = [0.5]
for i in range(1,10):
    geo_freq.append(geo_freq[i-1] * 0.5)

def plot_freq(est_freq, delta, ref_genomes, exp_freq=None, log_base=2):
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

    if exp_freq: #compute expected and estimated
        lab = est_freq.keys()
        
        total = sum(est_freq.values())
        for i in est_freq:
            est_freq[i] /= float(total)
        
        total = sum(exp_freq.values())
        for i in exp_freq:
            exp_freq[i] /= float(total)
            
        d = sorted(exp_freq.iteritems(), key=operator.itemgetter(1), reverse=True)
        lab = [i[0] for i in d]
        exp_to_plot = [math.log(exp_freq[l], log_base) for l in lab]
        est_to_plot = [math.log(est_freq[l], log_base) for l in lab]
        
    else: #compute estimated only
        total = sum(est_freq.values())
        for i in est_freq:
            est_freq[i] /= float(total)
        
        d = sorted(est_freq.iteritems(), key=operator.itemgetter(1), reverse=True)
        lab = [i[0] for i in d]
        est_to_plot = [math.log(est_freq[l], log_base) for l in lab]

    po = open('est_freq.pck', 'w')
    pickle.dump(est_freq, po)
    po.close()
    x1 = np.arange(1, len(est_to_plot)+1)
    

    # plot 
    plt.figure(1, figsize=(12, 12))
    plt.subplot(211, axisbg='white')
    
    plt.figure(figsize=(8, 8))
    ax = plt.axes([0.11, 0.11, 0.85, 0.85])
    if exp_freq:
        plt.plot(x1, exp_to_plot, '^-', label='expected', markersize=14, alpha=0.7, color='c')
    
    # plt.plot(x2, geo_values, '-', alpha=0.3, label='geometric series', color='r')
    # plt.plot(x2, log_act, '^', label='present', markersize=14, alpha=0.7, color='c')
    plt.plot(x1, est_to_plot, 'o', label='estimated', markersize=14, alpha=0.7, color='m')


    # setup the plot
    plt.axis(xmin=0.4, ymax=0.5, xmax=max(x1)+0.5)    
    plt.xticks(range(1, 1+len(lab)), lab)
    plt.setp( plt.gca().get_xticklabels(), rotation=30, horizontalalignment='right', size='x-small')
#    plt.yticks(np.arange(0, -11, -1))
    plt.grid(alpha = 0.3)
    plt.ylabel(r'frequency [$\log_%d$]' % log_base)
    plt.xlabel('reference')
    plt.legend(loc='upper right', numpoints=1)
    
    params = {'legend.fontsize': 8, 'font.size': 12}
    plt.rcParams.update(params)
    
    imtype = 'pdf'
    plt.savefig('fig_closest.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)
    
    plt.show()
    
    return




def align_info(al_pair):
    """
    
    """
    from Bio import SeqIO
    from Bio import SeqUtils
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from pythonlib.MarkxIO import Markx10Iterator
    
    f_forward = open(al_pair[0])
    f_reverse = open(al_pair[1])
    
    forwardaligniter = Markx10Iterator(f_forward)
    reversealigniter = Markx10Iterator(f_reverse)
    score = {}
    
    while True:
        try:
            f_align = forwardaligniter.next()
            r_align = reversealigniter.next()
        except:
            break
        if f_align is None or r_align is None:
            break
        
        f_id = f_align.get_all_seqs()[1].id.split(':')[-1]
        r_id = r_align.get_all_seqs()[1].id.split(':')[-1]
        assert f_id == r_id, 'same seq back and forward %s %s' % (f_id, r_id)
        t_scores = float(f_align._annotations['sw_ident']), float(r_align._annotations['sw_ident'])
        score[f_id] = max(t_scores)
    
    return score



def reads2clones_align(ref_genome):
    """
    
    """
    from Bio.Application import generic_run
    from pythonlib.cmline import NeedleCommandline
    go = 6.0
    ge = 3.0
    Verbose = False
    
    needle_run = NeedleCommandline()
    needle_run.set_parameter('-asequence', ref_genome)
    needle_run.set_parameter('-gapopen', go)
    needle_run.set_parameter('-gapextend', ge)
    needle_run.set_parameter('-aformat', 'markx10')
    needle_run.set_parameter('-usa', True)
    needle_run.set_parameter('-auto', True)

    
    # forward
    needle_run.set_parameter('-bsequence', 'tmp_reads_f.fas')
    outfile_1 = ref_genome + '_f.needle'
    needle_run.set_parameter('-outfile', outfile_1)
    
    result_1, messages_1, errors_1 = generic_run(needle_run)
    
    # reverse
    needle_run.set_parameter('-bsequence', 'tmp_reads_r.fas')
    outfile_2 = ref_genome + '_r.needle'
    needle_run.set_parameter('-outfile', outfile_2)

    result_2, messages_2, errors_2 = generic_run(needle_run)
    
    for ar in result_1.available_results():
        logfun.info(ar + result_1.get_result(ar))

    for m in messages_1.readlines():
        logfun.debug(m)
    
    for e in errors_1.readlines():
        logfun.debug(e)
    
    for ar in result_2.available_results():
        logfun.info(ar + result_2.get_result(ar))
    
    for m in messages_2.readlines():
        logfun.debug(m)
    
    for e in errors_2.readlines():
        logfun.debug(e)
        
    return outfile_1, outfile_2


def main():
    '''
    '''
    from Bio import SeqIO
    import pickle
    import operator
    import heapq
    import tempfile
    from multiprocessing import Pool
    
    args = sys.argv
    logfun.info(' '.join(args))
    
    try:
        reads_file, clones_file = args[1].rstrip('/'), args[2]
    except:
        sys.exit('usage: closest_reads.py reads_file clones_file')
    
    reads_dict = {}
    
    f_fasta = open(reads_file)
    seqlist = list(SeqIO.parse(f_fasta, 'fasta'))
    countreads = len(seqlist)
    
    # forward...
    logfun.info('write reads forward')
    f_fasta_forward_filename = 'tmp_reads_f.fas'
    f_fasta_forward = open(f_fasta_forward_filename, 'w')
    SeqIO.write(seqlist, f_fasta_forward, 'fasta')
    f_fasta_forward.close()

    # ...and reverse
    logfun.info('write reads backward')
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
    ref_genomes = [ s.id for s in clones_list ]
    
    try:
        exp_freq = {}
        for s in clones_list:
            exp_freq[s.id]= float(s.description.split('|')[-1])
        logfun.info('found expected frequencies in the reference file')
    except:
        exp_freq = None
        logfun.info('expected frequencies in the reference file not found')
        pass
    
    try:
        wh = open('%s-closest.dict' % reads_file.lstrip('../').replace('/', '-'))
        reads_dict = pickle.load(wh)
        logfun.info('dictionary of closest references found, loaded')
    except:
        logfun.info('dictionary of closest references NOT found, aligning all to all')
        refs = []
        c_names = []
        for c in clones_list:
            th = tempfile.NamedTemporaryFile(delete=False)
            th.write('>%s\n' % c.id)
            th.write('%s' % c.seq.tostring())
            th.close()
            refs.append(th.name)
            c_names.append(c.id)
        
        # align to all clones
        logfun.info('now to all clones')
        pool = Pool()
        al_results = pool.map(reads2clones_align, refs)
        pool2 = Pool()
        logfun.info('parsing the alignment results')
        t_d = pool2.map(align_info, al_results)
        
        for cn, rd in zip(c_names, t_d):
            for k in rd:
                try:
                    reads_dict[k][cn] = rd[k]
                except KeyError:
                    reads_dict[k] = {}
                    reads_dict[k][cn] = rd[k]
        
        wh = open('%s-closest.dict' % reads_file.lstrip('../').replace('/', '-'), 'w')
        pickle.dump(reads_dict, wh)
        wh.close()
        
        # remove references tmp files
        for r in refs:
            os.remove(r)
        # remove alignment tmp files
        for r in al_results:
            os.remove(r[0])
            os.remove(r[1])
            
    r_keys = reads_dict.keys()
    logfun.info('Lost %d reads' % (len(seqlist) - len(r_keys)))
    
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
    logfun.info(clones_freq)
    logfun.info(str(discarded) + ' reads had two matches')
    plot_freq(clones_freq, delta, ref_genomes, exp_freq)

if __name__ == '__main__':
    main()
