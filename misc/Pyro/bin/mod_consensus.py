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

thresh = 1./(B-1)


def int_histogram(data, start, stop):
    """
    Returns the histogram vector for data, where data only contains integer
    """
    for i in data:
        if type(i) != type(1):
            sys.exit('data must be integer')
    
    histo = [0]*(stop + 1 - start)
    step  = 1
    
    for d in data:
        if start <= d and d <= stop:
            histo[d-start] += 1
    
    histo_data = [ [i+start, histo[i]] for i in range(stop + 1 - start) ]
    
    return histo_data


def far_consensus(in_file):
    """
    Parse reads from a file with aligned reads in fasta format
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna
    import random
    
    handle = open(in_file)
    format = 'fasta'
    print >> sys.stderr, 'Parsing aligned reads from file', in_file
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
            consensus = [None] * (5000) #2*gen_length + 1)
            for i in range(5000): #gen_length + 1):
                consensus[i] = [0] * B
        
        gen_length = len(s.seq)
        
        #if gen_length != old_length:
        #    print 'All reads must have the same length'
        #    print id, 'is', gen_length,' bps long, different from the previous one'
        #    sys.exit()
            
        if id in mod_reads.keys():
            sys.exit('reads should all have different names')
            
        mod_reads[id] = []
        # sequences go from 1 to gen_length
        mst = s.seq.tostring()
        
        mls = list( mst.rstrip('-') )
        mstop = len (mls)
        flanking = True
        mod_reads[id].append('X')

        j = 0
        for c in mls:
            j += 1
            
            
            if c != '-' and flanking:
                mstart = j
                flanking = False

            if flanking == True:
                mod_reads[id].append('X')
            else:
                mod_reads[id].append(c)
                base = d2i[c]
                consensus[j][base] += 1

        #        print mstart, mstop
        mstart_arr.append(mstart)
        # all sequences parsed

    cseq = []
    
    coverage = [0]*len(consensus)
    j = 1
    for pos in consensus[1:]:
        tmax = max(pos)
        tot_bases = sum(pos)
        coverage[j] = tot_bases
        j += 1
        if tmax > thresh * tot_bases:
            b = []
            for c in range(len(pos)):
                if pos[c] == tmax:
                    b.append(c)
            cseq.append(i2d[random.choice(b)])
        else:
            pass
#            cseq.append('N')

    final_cons = []
    for cb in cseq:
        if cb != '-':
            final_cons.append(cb)
    print 'Length of the consensus', len(final_cons)
    
    cons_str = ''.join(final_cons)
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
    return

def main():
    in_file = sys.argv[1]
    far_consensus(in_file)

if __name__ == '__main__':
    main()
