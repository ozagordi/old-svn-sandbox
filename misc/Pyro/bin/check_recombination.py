#!/usr/bin/env python

import sys
import os
import re
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib

Verbose = True

def unzip(a):
    return tuple(map(list,zip(*a)))

def calc_ident(a, b):
    
    a = a.rstrip('-')
    b = b.rstrip('-')
    l = max( len(a)-len(a.lstrip('-')), len(b)-len(b.lstrip('-')) )
    a = a[l:]
    b = b[l:]
    
    al2num = ( lambda c: 1 if c[0] == c[1] else 0 )
    ones = [ al2num(c) for c in zip(a, b) ]
    
    return float(sum(ones))/len(ones)


def invert_keys(al_set):
    
    c_ids = al_set.keys()
    r_ids = al_set[c_ids[0]].keys()
    
    ial_set = {}
    for c in c_ids:
        for r in r_ids:
            try:
                ial_set[r][c] = al_set[c][r]
            except KeyError:
                ial_set[r] = {}
                ial_set[r][c] = al_set[c][r]
    
    return ial_set


def plot_amb_skewness(skewness, skewness_amb):
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
    legend = False
    plot_log= False
    nraw, bins, patches = plt.hist(skewness, 100, normed=0, facecolor='m', alpha=0.5, label='unambiguous', log=plot_log)
    nraw, bins, patches = plt.hist(skewness_amb, 100, normed=0, facecolor='c', alpha=0.5, label='ambiguous', log=plot_log)
    plt.xlabel('best two score delta')
    
    params = {'legend.fontsize': 14, 'font.size': 20}
    plt.rcParams.update(params)
    
    plt.legend() #loc = 'lower right', numpoints = 1)
    imtype = 'pdf'
    plt.savefig('fig_recomb_skew.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()
    
    


def plot_amb_scatter(tot_score, sum_score, tot_score_amb, sum_score_amb):
    """
    
    """
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        import numpy as np
    except:
        print 'exit, could not import matplotlib'
        sys.exit()

    import math
    legend = True
    lower =  int(min(tot_score + sum_score + tot_score_amb + sum_score_amb))
    upper =  int(max(tot_score + sum_score + tot_score_amb + sum_score_amb) + 0.5)
    print lower, upper
    
    line = np.arange(lower, upper+1)
#    line2 = [2*i for i in line]
    plt.scatter(tot_score_amb, sum_score_amb, label='ambiguous', color='c')
    plt.scatter(tot_score, sum_score, label='unambiguous', color='m')
#    plt.hexbin(tot_score, sum_score, cmap=cm.jet)

    plt.plot(line, line, color='b')
    plt.xlabel('score of the entire read')
    plt.ylabel('sum of the scores')
    
    params = {'legend.fontsize': 12, 'font.size': 16}
    plt.rcParams.update(params)
    plt.legend(loc = 'lower right', numpoints = 1)

    imtype = 'pdf'
    plt.savefig('fig_recomb_scatter.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()
    
    
    return


def plot_amb_hist(total_delta, amb_delta):
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
    plot_log = True
    legend = True
    
    nraw, bins, patches = plt.hist(total_delta, 100, normed=0, facecolor='m', alpha=0.5, label='unambiguous', log=plot_log)
    nraw, bins, patches = plt.hist(amb_delta, 100, normed=0, facecolor='c', alpha=0.5, label='ambiguous', log=plot_log)
    plt.xlabel('best two score delta')
    plt.xlim(xmin=-0.4, xmax=0.6)
    plt.xticks(np.arange(-0.4, 0.6, 0.05), rotation=90)

    params = {'legend.fontsize': 14, 'font.size': 20}
    plt.rcParams.update(params)
    if legend:
        plt.legend()

    imtype = 'pdf'
    plt.savefig('fig_recomb_hist.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()
    
    
    return



def main():
    
    from Bio import SeqIO
    import cPickle
    from multiprocessing import cpu_count
    from pythonlib import Alignment
    from pythonlib import pprocess
    import operator
    import heapq
    import time
    
    try:
        n_proc = cpu_count()
    except NotImplementedError:
        n_proc = 4
    
    HPP = cPickle.HIGHEST_PROTOCOL
    min_len = 200
    args = sys.argv
    
    try:
        reads_file, clones_file = args[1].rstrip('/'), args[2]
    except:
        sys.exit('usage: check_recombination.py reads_file clones_file')
    
    reads_dict = {}
    reads_dict_1 = {}
    reads_dict_2 = {}
    gaps_dict = {}
    gaps_dict_1 = {}
    gaps_dict_2 = {}
    
    f_fasta = open(reads_file)
    tmp_seqlist = list(SeqIO.parse(f_fasta, 'fasta'))
    f_fasta.close()
    countreads = len(tmp_seqlist)
    print >> sys.stderr, ' %d reads in the original file '.center(60, '-') % countreads
    
    seqlist = [ s for s in tmp_seqlist if len(s) > min_len ]
    
    print >> sys.stderr, ' %d reads are longer than %d '.center(60, '-') % (len(seqlist), min_len)
    
    try:
        t = time.time()
        print >> sys.stderr, ' loading file '.center(60, '-')
        wh = open('%s-check-reads_total.pck' % reads_file.replace('.', 'U').replace('/', '-'))
        al_set_total = cPickle.load(wh)
        wh.close()
        
        print >> sys.stderr, ' loading file '.center(60, '-')
        wh = open('%s-check-reads_1.pck' % reads_file.replace('.', 'U').replace('/', '-'))
        al_set_1 = cPickle.load(wh)
        wh.close()
        
        print >> sys.stderr, ' loading file '.center(60, '-')
        wh = open('%s-check-reads_2.pck' % reads_file.replace('.', 'U').replace('/', '-'))
        al_set_2 = cPickle.load(wh)
        wh.close()
        print >> sys.stderr, ' pickle objects loaded in %d seconds '.center(60, '-') % (time.time() - t)
        
    except:
        print >> sys.stderr, ' pickle objects not found, aligning '.center(60, '-')
        # reads are considered already aligned
        f_fasta_forward_filename = 'tmp_reads.fas'
        f_fasta_forward = open(f_fasta_forward_filename, 'w')
        SeqIO.write(seqlist, f_fasta_forward, 'fasta')
        f_fasta_forward.close()

        # split in 2, first segment
        f_fasta = open(f_fasta_forward_filename)
        tmp = list(SeqIO.parse(f_fasta, 'fasta'))
        f_fasta.close()
        for seq in tmp:
            l = len(seq)
            middle = int(float(l)/2)
            seq.seq = seq.seq[:middle]
        out_file = 'tmp_reads_1.fas'
        f_fasta_forward = open(out_file, 'w')
        SeqIO.write(tmp, f_fasta_forward, 'fasta')
        f_fasta_forward.close()
        del tmp
    
        # split in 2, second segment
        f_fasta = open(f_fasta_forward_filename)
        tmp = list(SeqIO.parse(f_fasta, 'fasta'))
        f_fasta.close()
        for seq in tmp:
            l = len(seq)
            middle = int(float(l)/2)
            seq.seq = seq.seq[middle:]
        out_file = 'tmp_reads_2.fas'
        f_fasta_forward = open(out_file, 'w')
        SeqIO.write(tmp, f_fasta_forward, 'fasta')
        f_fasta_forward.close()
        del tmp
    
        # clones
        hc = open(clones_file)
        clones_list = list((SeqIO.parse(hc, 'fasta')))
        
        i = 0
        tmp_files = []
        for c in clones_list:
            print >> sys.stderr, ' clone %s '.center(60, '-') % c.id
            tfn = 'tmp%d.fas' % i
            tmp_files.append(tfn)
            th = open(tfn, 'w')
            th.write('>%s\n' % c.id)
            th.write('%s' % c.seq.tostring())
            th.close()
            i += 1
        # parallelism
        queue = pprocess.Queue(limit=n_proc)
            
        ral_par = queue.manage(pprocess.MakeParallel(Alignment.needle_align))
        for tf in tmp_files:
            # align total
            ral_par(tf, 'tmp_reads.fas', tf.split('.')[0] + '-total.needle')

        for tf in tmp_files:
            # align total
            ral_par(tf, 'tmp_reads_1.fas', tf.split('.')[0] + '-1.needle')

        for tf in tmp_files:
            # align total
            ral_par(tf, 'tmp_reads_2.fas', tf.split('.')[0] + '-2.needle')

        for res in queue:
            if Verbose:
                print >> sys.stderr, res[0], res[1]

        # alignment with whole reads
        files = [f for f in os.listdir('./') if f.startswith('tmp') and f.endswith('-total.needle')]        
        al_set_total = Alignment.alignfile2set(files, 'total_read', 6.0, 3.0)
        
        wh = open('%s-check-reads_total.pck' % reads_file.replace('.', 'U').replace('/', '-'), 'w')
        cPickle.dump(al_set_total, wh, HPP)
        wh.close()
        
        # alignment with first half
        files = [f for f in os.listdir('./') if f.startswith('tmp') and f.endswith('-1.needle')]        
        al_set_1 = Alignment.alignfile2set(files, 'total_read', 6.0, 3.0)
        
        wh = open('%s-check-reads_1.pck' % reads_file.replace('.', 'U').replace('/', '-'), 'w')
        cPickle.dump(al_set_1, wh, HPP)
        wh.close()
        
        # alignment with second half
        files = [f for f in os.listdir('./') if f.startswith('tmp') and f.endswith('-2.needle')]        
        al_set_2 = Alignment.alignfile2set(files, 'total_read', 6.0, 3.0)
        
        wh = open('%s-check-reads_2.pck' % reads_file.replace('.', 'U').replace('/', '-'), 'w')
        cPickle.dump(al_set_2, wh, HPP)
        wh.close()
        
        # except ends
    

    count = 0
    for i in al_set_total:
        for j in al_set_total[i]:
            count += 1
            

    ial_set_total = invert_keys(al_set_total)
    ial_set_1 = invert_keys(al_set_1)
    ial_set_2 = invert_keys(al_set_2)
    
    del al_set_total
    del al_set_1
    del al_set_2
    
    r_keys = ial_set_total.keys()
    lost =  len(seqlist) - len(r_keys)
    assert lost == 0, 'lost' + str(lost) + 'reads'
    
    delta = []
    ambiguous = 0
    total_delta = []
    amb_delta = []
    tot_score = []
    sum_score = []
    tot_score_amb = []
    sum_score_amb = []
    outliers = []
    best_out = {}
    thresh_inc = 0.05
    print >> sys.stderr, 'Total reads', len(r_keys)
    
    skewness = []
    skewness_amb = []
    
    for k in ial_set_total:
        total = ial_set_total[k]
        s1 = ial_set_1[k]
        s2 = ial_set_2[k]
        
        l_tot = [(s[0], s[1].score) for s in total.iteritems()]        
        best2_total = heapq.nlargest(2, iter(l_tot), operator.itemgetter(1))
        
        l_1 = [(s[0], s[1].score) for s in s1.iteritems()]
        best2_s1 = heapq.nlargest(2, iter(l_1), operator.itemgetter(1))
        
        l_2 = [(s[0], s[1].score) for s in s2.iteritems()]
        best2_s2 = heapq.nlargest(2, iter(l_2), operator.itemgetter(1))
        
        clone_t = best2_total[0][0]
        clone_1 = best2_s1[0][0]
        clone_2 = best2_s2[0][0]

        ial_set_total[k][clone_t].summary()
        ial_set_1[k][clone_1].summary()
        ial_set_2[k][clone_2].summary()
        
        len_t = ial_set_total[k][clone_t].stop - ial_set_total[k][clone_t].start + 1
        len_1 = ial_set_1[k][clone_1].stop - ial_set_1[k][clone_1].start + 1
        len_2 = ial_set_2[k][clone_2].stop - ial_set_2[k][clone_2].start + 1
                
        bt = best2_total[0][1]/len_t
        b1 = best2_s1[0][1]/len_1
        b2 = best2_s2[0][1]/len_2
        
        relative_gain = (b1 + b2 - 2*bt)/(b1 + b2)
        
        # if 0.4 < relative_gain and relative_gain < 0.8:
        if abs(len_t - len_1 - len_2) > 5000:
            print best2_total
            print best2_s1
            print best2_s2
            
            print len_t
            print len_1
            print len_2
            print ial_set_total[k][clone_t].seq_a
            print ial_set_total[k][clone_t].seq_b
            
            print ial_set_1[k][clone_1].seq_a
            print ial_set_1[k][clone_1].seq_b
            
            print ial_set_2[k][clone_2].seq_a
            print ial_set_2[k][clone_2].seq_b
            sys.exit()

        if best2_total[0][0] != best2_s1[0][0] or best2_total[0][0] != best2_s2[0][0]:
            amb_delta.append(relative_gain)
            skewness_amb.append( abs(b1-b2)/(bt) )
            tot_score_amb.append(bt)
            sum_score_amb.append((b1 + b2)/2)
            ambiguous += 1
            if relative_gain > 0.05 and b1 > 4.5 and b2 > 4.5:
                tk = k#.split('#')[0]
                outliers.append(tk)
                best_out[tk] = [ best2_total[0][0], best2_s1[0][0], best2_s2[0][0] ]
        else:
            total_delta.append(relative_gain)
            skewness.append( abs(b1-b2)/(bt) )
            tot_score.append(bt)
            sum_score.append((b1 + b2)/2)

            
    #    print >> sys.stderr, discarded, 'reads had two matches'
    print >> sys.stderr, ambiguous, 'potentially ambiguous'
    print >> sys.stderr, len(outliers), 'outliers'
    # write the outliers
    handle = open(reads_file)
    tmp_dict = SeqIO.to_dict(SeqIO.parse(handle, 'fasta'))
    reads_dict = {}
    for k in tmp_dict.keys():
        k1 = k.split('#')[0]
        reads_dict[k1] = tmp_dict[k]
    out_list = [ reads_dict[r] for r in outliers ]
    out_handle = open('outliers.fas', 'w')
    SeqIO.write(out_list, out_handle, 'fasta')
    out_handle.close()
#    for k in best_out:
#        print >> sys.stderr, k, best_out[k]
    plot_amb_hist(total_delta, amb_delta)
    #plot_amb_scatter(tot_score, sum_score, tot_score_amb, sum_score_amb)
    # plot_amb_skewness(skewness, skewness_amb)
if __name__ == '__main__':
    main()
