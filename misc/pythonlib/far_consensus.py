def far_consensus(in_file):
    """
    Parse reads from a file with aligned reads in fasta format
    """
    from Bio import SeqIO
    import sys
    import random
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
            consensus = [ [0] * B for i in range(2*gen_length) ]
        
        gen_length = len(s.seq)
                    
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

    final_cons = cseq
    print >> sys.stderr, 'Length of the consensus', len(final_cons)
    
    cons_str = ''.join(final_cons)
    
    return cons_str, coverage, mstart_arr

def write_consensus(cons_str):
    
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna
    
    cons = SeqRecord(Seq(str(cons_str), generic_dna), id=in_file,
                     description='')
    cons_rec = [cons]
    ohandle = open("consensus.faa", "w")
    SeqIO.write(cons_rec, ohandle, "fasta")
    handle.close()
    
