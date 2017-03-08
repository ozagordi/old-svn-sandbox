#!/usr/bin/env python
''' Now developed with Aptana'''

import sys
import os
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib
from Bio import SeqIO

def shannon_entropy(d):
    '''Shannon entropy of a dictionary'''
    from math import log
    if d == {}: return None
    sv = sum(d.values())
    l = [float(v)/sv for v in d.values()]
    return sum([-p*log(p, 2) for p in l if p > 0])

def far_counts(far_handle):
    '''Takes a far file and returns the vector
    of Shannon entropy
    '''
    seqs = list(SeqIO.parse(far_handle, 'fasta'))
    L = len(seqs[0])
    counts = [{} for i in range(L)]
    assert len(counts) == L
    for s in seqs:
        seq_start = False
        for i, b in enumerate(s.seq.tostring().rstrip('-')):
            if b != '-' and s[i-1] == '-':
                seq_start = True
            if seq_start:
                counts[i][b] = counts[i].get(b, 0) + 1
                
    return counts

def check_alignment():
    '''This counts the number of internal long indels 
    '''
    from pythonlib.MarkxIO import Markx10Iterator
    
    f_forward = open('tmp_align_f.needle')
    f_reverse = open('tmp_align_r.needle')
    
    forwardaligniter = Markx10Iterator(f_forward)
    reversealigniter = Markx10Iterator(f_reverse)
    count_forward = 0
    count_reverse = 0
    count = {}
    read_count = 0
    # iterates through the alignments
    while True:
        try:
            f_align = forwardaligniter.next()
            r_align = reversealigniter.next()
        except:
            break
        if f_align is None or r_align is None:
            break
        
        assert f_align.get_all_seqs()[1].id == r_align.get_all_seqs()[1].id, 'same seq back and forward'
        descr = f_align.get_all_seqs()[1].id
                
        if float(f_align._annotations['sw_score']) > float(r_align._annotations['sw_score']):
            tmp = f_align.get_seq_by_num(1).tostring().upper()
            refseq = f_align.get_seq_by_num(0).tostring().upper()
            count_forward += 1
        else:
            tmp = r_align.get_seq_by_num(1).tostring().upper()
            refseq = r_align.get_seq_by_num(0).tostring().upper()
            count_reverse += 1
            
        if tmp.strip('-').count('-----') > 0:
            read_count += 1
        
        for i in range(2, 31):
            word = str('-'*i)
            c = tmp.strip('-').count(word)
            count[i] = count.get(i, 0) + c

    print read_count
    return count

def plot_entropy(counts, thresh=0.5):
    ''' This takes the counts as input, computes
    the entropy and plots it. At the same time it
    saves a csv file with sequence information
    '''
    import matplotlib
    matplotlib.use('pdf')   
    import matplotlib.pyplot as plt
    import csv
    import numpy as np
    bases = ['A', 'C', 'G', 'T']
    entropy = np.array([shannon_entropy(d) for d in counts])
    pos = np.arange(1, len(entropy)+1)
    plt.figure(figsize=(12, 2))
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.87)
                #wspace=None, hspace=None)
    plt.title('Manhattan plot of Shannon entropy')
    plt.xlim(xmin=0, xmax =len(entropy)+1)
    plt.ylim(ymin=-0.05, ymax=2.1)
    
    
    plt.scatter(pos, entropy)
    writer = csv.writer(open('mutations.csv', 'wb'), dialect='excel')
    writer.writerow(['pos\\bases'] + bases)
    for i, e in enumerate(entropy):
        p = pos[i]
        if e > thresh:
            counts_here = [float(counts[i].get(b, 0))/sum(counts[i].values()) for b in bases]
            counts_rounded = [round(x, 2) for x in counts_here]
            writer.writerow([p] + counts_rounded)
            xy = p, e
            #xytext = p, e + 0.02
            plt.annotate(str(p), xy, xycoords='data',
                        xytext=(-8, 20), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->"))
            

    #plt.show()
    imtype='pdf'
    plt.savefig('shann_entr.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                orientation='portrait', papertype=None, format=imtype,\
                transparent=False)


def main():
    
    import numpy as np
    thresh = 0.5
    args = sys.argv
    far_file = args[1]
    
    counts = far_counts(open(far_file))
    plot_entropy(counts)
    
if __name__ == '__main__':
    from StringIO import StringIO
    data = ">read1\n--ACGT--\n>read2\n-CAGGT\n"
    far_counts(StringIO(data))
    main()
