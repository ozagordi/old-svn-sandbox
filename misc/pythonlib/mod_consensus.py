#!/usr/bin/env python
import sys

d2i = {
    'A': 0, 'a': 0,
    'C': 1, 'c': 1,
    'G': 2, 'g': 2,
    'T': 3, 't': 3,
    '-': 4,
    'N': 5, 'n': 5
    }
i2d = ['A', 'C', 'G', 'T', '-', 'N']
B = len(i2d)
f_code = {}
f_code['-'] = ['-']
f_code['A'] = ['A']
f_code['C'] = ['C']
f_code['G'] = ['G']
f_code['T'] = ['T']
f_code['R'] = ['G', 'A']
f_code['Y'] = ['T', 'C']
f_code['K'] = ['G', 'T']
f_code['M'] = ['A', 'C']
f_code['S'] = ['G', 'C']
f_code['W'] = ['A', 'T']
f_code['B'] = ['G', 'T', 'C']
f_code['D'] = ['G', 'A', 'T']
f_code['H'] = ['A', 'C', 'T']
f_code['V'] = ['G', 'C', 'A']
f_code['N'] = ['A', 'G', 'C', 'T']

thresh = 0.4 #1./(B-1)

def choose(good):
    '''
    '''
    
    for k, v in f_code.items():
        if set(good) == set(v):
            return k


def far_consensus_old(handle):
    """
    Parse reads from a file with aligned reads in fasta format
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna
    import random
    import operator
    # handle = open(in_file)
    format = 'fasta'
    # print >> sys.stderr, 'Parsing aligned reads from file', in_file
    mod_reads = {}
    alignment = SeqIO.parse(handle, format)
    al_len = 0
    mstart_arr = []
    for s in alignment:
        id = s.id
        al_len += 1
        
        #  print >> sys.stderr, '\x1B[1A\x1B[2K parsing sequence', id, al_len
        
        try: # check all same length
            old_length = gen_length
        except: # this is done only for the first sequence
            gen_length = len(s.seq)
            old_length = gen_length
            consensus = [None] * (2*gen_length + 1)
            for i, c in enumerate(consensus): #range(5000): #gen_length + 1):
                consensus[i] = {}
#                c = [0] * B
        
        gen_length = len(s.seq)
        
        if gen_length != old_length:
            print 'All reads must have the same length'
            print id, 'is', gen_length,' bps long, different from the previous one'
            sys.exit()
            
        if id in mod_reads.keys():
            sys.exit('reads should all have different names')
            
        mod_reads[id] = []
        # sequences go from 1 to gen_length
        mst = s.seq.tostring()
        
        mls = list( mst.rstrip('-') )
        mstop = len (mls)
        flanking = True
        mod_reads[id].append('X')
        
        for j, c in enumerate(mls):
            if c != '-' and flanking:
                mstart = j
                flanking = False

            if flanking == True:
                mod_reads[id].append('X')
            else:
                mod_reads[id].append(c)
                base = d2i[c]
                consensus[j][c] = consensus[j].get(c, 0) + 1

        mstart_arr.append(mstart)
        # all sequences parsed

    cseq = []
    
    for pos in consensus[1:]:
        ord = sorted(pos.items(), key=operator.itemgetter(1), reverse=True)
        if ord != []:
            cseq.append(ord[0][0])

    return ''.join(cseq)

'''
cons = SeqRecord(Seq(str(cons_str), generic_dna), id=in_file,
description='')
cons_rec = [cons]
ohandle = open("consensus.faa", "w")
SeqIO.write(cons_rec, ohandle, "fasta")
handle.close()

mstart_histo = int_histogram(mstart_arr, 0, gen_length)
cov_h = open('match_start.dat', 'w')

for i in mstart_histo:
cov_h.write('%d\t%d\n' % (i[0], i[1]))
return cons_str
'''
def main():
    in_file = iopen(sys.argv[1])
    far_consensus_old(in_file)

if __name__ == '__main__':
    main()
