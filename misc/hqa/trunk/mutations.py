#!/usr/bin/env python

import sys
import os
# showing moritz
# comment 2
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib
from Bio import SeqIO


import time
import optparse
import random, math

#ref_genome = os.path.join(homedir, 'References/HIV-HXB2.fasta')
ref_genome = os.path.join(homedir, 'References/105REF_capital.fasta')

def mutation_error(aligned_reads):
    """Compute the error rate from stop codons analysis. See:
    Cuevas et al. The effect of ribavirin on the mutation rate
    and spectrum of Hepatitis C virus in vivo.
    Journal of Virology (2009) pp. 
    """
    nsmt = ['AAA', 'AAG', 'AGA',
            'CAA', 'CAG', 'CGA',
            'GAA', 'GAG', 'GGA',
            'TAC', 'TAT', 'TCA', 'TCG', 'TGC', 'TGG', 'TGT', 'TTA', 'TTG']
    stops = ['TAG', 'TGA', 'TAA']
    
    correction = {
        'AAA': 3, 'AAG': 3, 'AGA': 3,
        'CAA': 3, 'CAG': 3, 'CGA': 3,
        'GAA': 3, 'GAG': 3, 'GGA': 3,
        'TAC': 1.5, 'TAT': 1.5, 'TCA': 1.5, 'TCG': 3, 'TGC': 3, 'TGG': 3, 'TGT': 3, 'TTA': 1.5, 'TTG': 3
        }
    
    ref = list(SeqIO.parse(open(ref_genome), 'fasta'))[0]
    tot = 0
    n = 0
    corr = 0.
    positions = [0]*(len(ref.seq)/3)
    exclude = [95, 96, 195, 285, 316, 317, 330, 336, 339, 163, 164, 171, 309, 323, 329, 340] #288, 949, 950, 951]
    
    for r in aligned_reads:
        start = aligned_reads[r][1]
        if start%3 == 0:
            read = aligned_reads[r][0]
        elif start%3 == 1:
            read = aligned_reads[r][0][2:]
        elif start%3 == 2:
            read = aligned_reads[r][0][1:]
            
        end = len(read) - len(read)%3 - 1
        
        for i in range(0, end, 3):
            if read[i:i+3] in stops:
                if (start+i)/3 not in exclude:
                    positions[(start+i)/3] += 1
                    parent = ref.seq[start+i:start+i+3].tostring()
                    sc = read[i:i+3]
                    tot += 1
                    try:
                        corr += correction[parent]
                    except:
                        pass
            if read[i:i+3] in nsmt:
                n +=1
    to_p = [ i[0] for i in enumerate(positions) if i[1] >= 1 and i[1] <= 5]
    print to_p
    print len(to_p)
    print tot, corr, n
    
    print >> sys.stderr, 'mu is', '%6.2E' % (corr/n)
    
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import numpy as np
    import scipy.stats as stats
    
    pos2 = [ i for i in positions if i < 20]
    
    ax = plt.subplot(211)
    plt.plot(positions)
    ax = plt.subplot(212)
    ax.set_yscale('log')
    plt.hist(pos2, 50)
    filename = 'histo_pos'
    imtype = 'pdf'
    plt.savefig('%s.%s' % (filename, imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)
    
    sys.exit()

def align_reads(filename):
    """reads the file with reads, align them with the reference,
    returns a dictionary with reads (in-dels are discarded)
    and starting position with respect to the reference
    """
    from pythonlib import EmbossStandalone
    from pythonlib.MarkxIO import Markx10Iterator
    
    needle_exe = 'needle'

    aligned_reads = {}

    f_fasta = open(filename)
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

    print >> sys.stderr, 'Found', countreads, 'reads'
    
    if not os.path.isfile('tmp_align_f.needle'):
        print >> sys.stderr, 'needle forward'
        EmbossStandalone.needle(needle_exe, ref_genome, f_fasta_forward_filename,
                                out='tmp_align_f.needle', gapopen=6.0, gapext=3.0,
                                aglobal3='False', adesshow3 = 'True')
    
    if not os.path.isfile('tmp_align_r.needle'):
        print >> sys.stderr, 'needle backward'
        EmbossStandalone.needle(needle_exe, ref_genome, f_fasta_reverse_filename,
                                out='tmp_align_r.needle', gapopen=6.0, gapext=3.0,
                                aglobal3='False', adesshow3 = 'True')
    
    f_forward = open('tmp_align_f.needle')
    f_reverse = open('tmp_align_r.needle')
    
    forwardaligniter = Markx10Iterator(f_forward)
    reversealigniter = Markx10Iterator(f_reverse)
    count_forward = 0
    count_reverse = 0
    
    while True:
        
        # pos += 1
        # print >> sys.stderr,  '\x1B[1A\x1B[2K', pos    
        try:
            f_align = forwardaligniter.next()
            r_align = reversealigniter.next()
        except:
            break
        
        if f_align is None or r_align is None:
            break
    
        assert f_align.get_all_seqs()[1].id == r_align.get_all_seqs()[1].id, 'same seq back and forward'
        
        this_id = f_align.get_all_seqs()[1].id
        
        if float(f_align._annotations['sw_score']) > float(r_align._annotations['sw_score']):
            tmp = f_align.get_seq_by_num(1).tostring().upper()
            refseq = f_align.get_seq_by_num(0).tostring().upper()
            count_forward += 1
        else:
            tmp = r_align.get_seq_by_num(1).tostring().upper()
            refseq = r_align.get_seq_by_num(0).tostring().upper()
            count_reverse += 1
        
        q_align_start = len(tmp) - len(tmp.lstrip('-'))
        q_align_end = len(tmp.rstrip('-'))
        
        m_align_start = len(refseq) - len(refseq.lstrip('-'))
        m_align_end = len(refseq.rstrip('-'))
        
        align_start = max(m_align_start, q_align_start)
        align_end = min(m_align_end, q_align_end)
        
        this_read = []
        for c in zip(refseq[align_start:align_end+1], tmp[align_start:align_end+1]):
            if c[0] != '-' and c[1] != '-':
                this_read.append(c[1])
            elif c[1] == '-':
                this_read.append(c[0])
            elif c[0] == '-':
                pass
        aligned_reads[this_id] = [''.join(this_read), align_start]

    
    return aligned_reads
        
def print_amino_al(aligned_reads, file):
    """ 
    """
    from Bio.Seq import translate
    
    aa_reads = {}
    start_pos = {}
    end = 0
    for r in aligned_reads:
        start = aligned_reads[r][1]
        
        if start%3 == 0:
            read = aligned_reads[r][0]
        elif start%3 == 1:
            read = aligned_reads[r][0][2:]
        elif start%3 == 2:
            read = aligned_reads[r][0][1:]
        
        aa_reads[r] = translate(read)
        
        start_pos[r] = start/3
        
        if start%3:
            start_pos[r] += 1
        
        if len(aa_reads[r]) + start_pos[r] > end:
            end = len(aa_reads[r]) + start_pos[r]

    print >> sys.stderr, 'sorting'
    items = [ (v, k) for k, v in start_pos.items() ]
    items.sort()
    
    out_file = file.split('.')[0] + '.paa'
    oh = open(out_file, 'w')
    
    for i in items:
        id = i[1]
        this_read = aa_reads[id]
        tl = len(this_read)
        to_write = '-'*start_pos[id] + this_read + '-'*(end - tl - start_pos[id])
        oh.write('>%s\n' % id)
        oh.write(to_write)
        oh.write('\n')
    oh.close()
    
    return


def codons(s, frame=0):
    """ check this function
    """
    end = len(s[frame:]) - len(s[frame:])%3 - 1
    codons = [ s[i:i+3] for i in range(frame, end, 3) ]
    
    return codons


def compute_codon_usage(aligned_reads):
    """
    """
    codon_usage = {}
    
    for r in aligned_reads:
        start = aligned_reads[r][1]
        
        if start%3 == 0:
            read = aligned_reads[r][0]
        elif start%3 == 1:
            read = aligned_reads[r][0][2:]
        elif start%3 == 2:
            read = aligned_reads[r][0][1:]
        
        this_codons = codons(read)
        for c in this_codons:
            try:
                codon_usage[c] += 1
            except KeyError:
                codon_usage[c] = 1.
    
    return codon_usage


def count_co(sites, file):
    """
    """
    from Bio import SeqIO
    h = open(file.split('.')[0] + '.paa')
    al = list(SeqIO.parse(h, 'fasta'))
    
    print >> sys.stderr, 'Found', len(al), 'protein translated reads'
    
    occ = []
    for s in al:
        occ.append( [s.seq[c-1] for c in sites] )

    co_freq = {}
    freq = [{} for k in sites ]
    
    for co in occ:
        word = ''.join(co)
        
        if word.endswith('-') or word.startswith('-'):
            continue
        
        try:
            co_freq[word] += 1
        except KeyError:
            co_freq[word] = 1.

        for c in enumerate(co):
            pos, let = c
            try:
                freq[pos][let] += 1
            except KeyError:
                freq[pos][let] = 1.

    
    co_this_sum = sum(co_freq.values())
    for k in co_freq:
        co_freq[k] /= co_this_sum
    
    this_sum = [ sum(k.values()) for k in freq ]
    assert this_sum[0] == this_sum[-1], 'sums must be equal'

    for pos in freq:
        for k in pos:
            pos[k] /= this_sum[0]

    print >> sys.stderr, 'Sum of amminoacids is: ', int(this_sum[0])
    print >> sys.stderr, 'Min frequency: ', 1/this_sum[0]
    print >> sys.stderr, 'Strange counts (if any):'
    odds = {}
    for k in co_freq:
        factor = 1.
        for c in enumerate(k):
            pos, let = c
            factor *= freq[pos][let]
        odds[k] = float(co_freq[k])/factor
        if abs(odds[k]) > 1.0:
            print >> sys.stderr, '  ', k, co_freq[k] * co_this_sum
    print >> sys.stderr,''
    
    if len(sites) > 1:
        return odds
    else:
        return freq[0]




def main():
    """run when not interactive
    """
    
    import pickle
    import math
#    from Bio import Translate
#    from Bio.Seq import Seq, translate
#    from Bio.Alphabet import IUPAC
#    standard_translator = Translate.unambiguous_dna_by_id[1]
    
    args = sys.argv
    
    try:
        file = args[1]
        pick_file = file.split('.')[0] + '.pck'
    except:
        print >> sys.stderr, 'usage: mutations.py reads_file'
    
    try:
        h = open(pick_file)
        aligned_reads = pickle.load(h)
        h.close()
    except:
        aligned_reads = align_reads(file)
        h = open(pick_file, 'w')
        pickle.dump(aligned_reads, h)
        h.close()
    
    if args[2] == 'stop':
        mu = mutation_error(aligned_reads)
        print mu
        sys.exit()

    out_file = file.split('.')[0] + '.paa'    
    if not os.path.exists(out_file):
        print_amino_al(aligned_reads, file)
    
    try:
        sites = map(int, args[2:])
    except:
        sites = None
    
    if sites:
        odds = count_co(sites, file)
    else:
        print >> sys.stderr, 'done'
    
    print >> sys.stderr, 'Odds are:'
    for d in odds:
        if len(d) > 1:
            lo = math.log(odds[d])
        else:
            lo = odds[d]
        print >> sys.stderr, '\t', d, '\t', '%9.6f' % lo, '*'*int(lo+0.5)
    
if __name__ == '__main__':
    main()
