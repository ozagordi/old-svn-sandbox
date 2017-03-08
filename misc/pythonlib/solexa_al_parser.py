import sys

from Bio import SeqIO
from Bio.Seq import Seq


verbose = False
thresh = 0
def_gl = 10000
def_sl = 36

def parse(al_file, ref_file=None, match_start=0):
    """
    Solexa align file
    """
    
    file_desc = {}
    
    try:
        rh = open(ref_file)
        rs = SeqIO.parse(rh, 'fasta').next().seq
        ref_sequence = rs.tostring().upper()
        genome_length = len(ref_sequence)+1
        print >> sys.stderr, 'reference length = ', genome_length - 1
    except:
        print >> sys.stderr, 'no reference assigned, length = ', def_gl
        ref_sequence = None
        genome_length = def_gl + 1
    
    handle = open(al_file)
    file_stem = al_file.rstrip('.txt')
    try:
        line, file_desc = parse_file_header(handle)
        seq_length = int(file_desc['SEQ_LENGTH'])
    except:
        seq_length = def_sl
        line = handle.readline()
    print >> sys.stderr, 'Number of cycles', seq_length
    if verbose:
        print file_desc
    assert not line.startswith('#'), (line, 'should not be here')

    [subst_read, subst_seq, read_start, coverage] = parse_alignment_info(handle, line, seq_length, file_stem, match_start, ref_sequence)
    return [subst_read, subst_seq, read_start, coverage]


def parse_file_header(handle) :
    """
    pfh
    """
    file_desc = {}
    line = handle.readline().rstrip('\n')
    assert line.startswith('#'), 'The file header should finish here'
    while line.startswith('#') :
        line = line.rstrip('\n')
        desc, value = line.split(' ')[0].lstrip('#'), '-'.join(line.split(' ')[1:])#.rstrip('\n')
        file_desc[desc] = value
        line = handle.readline().rstrip('\n')
    return line, file_desc


def parse_alignment_info(handle, line, seq_length, file_stem, match_start, ref_sequence):

    fgh = open('%s-pot_gaps.txt' % file_stem, 'w')
    assert not line.startswith('#'), line
    try:
        len_ref = len(ref_sequence)
    except:
        len_ref = def_gl
    pot_gaps = False
    subst_read = {}
    subst_seq  = {}
    read_start = [0]*len_ref
    coverage = [0]*(len_ref+seq_length+1)
    lnum = 0
    print 'len_ref', len_ref
    print ''
    while line:
        tl = line.rstrip('\n').split()
        print '\x1B[1A\x1B[2k', lnum
        lnum += 1
        try:
            if tl[2] != '1' or int(tl[3]) < match_start or int(tl[3]) > match_start + len_ref or float(tl[1]) < thresh: # only if...
                take = False
                sys.exit()
            else:
                take = True
        except:
            take = False
        
        if take:
            if len(tl[0]) != len(tl[5]):
                print '\t\t', len(tl[0]), len(tl[5]), '\n'
                print line
                sys.exit()
            if tl[4] == 'R': # match the reverse strand
                read = Seq(tl[0]).reverse_complement().tostring()
            else:
                read = tl[0]
            
            mstart = int(tl[3]) - match_start
            try:
                read_start[int(tl[3])] += 1
            except:
                pass
            try:
                seq = ref_sequence[mstart : mstart + seq_length]
            except:
                seq = tl[5].upper()
            if '.' not in read and 'N' not in read:
                for p in range(mstart, mstart + seq_length):
                    coverage[p] += 1
            else:
                for p in range(seq_length):
                    try:
                        coverage[mstart + p] += (read[p] != '.' and read[p] != 'N')
                    except:
                        pass
            pos = 1
            mpos = []
            for c in zip(read, seq):
                
                if tl[4] == 'F':
                    a_pos = pos
                elif tl[4] == 'R':
                    a_pos = -pos
                
                if c[0].upper() != c[1].upper():
                    mpos.append(a_pos)
                try:
                    subst_read[c[0]+c[1]][a_pos] += 1
                except:
                    subst_read[c[0]+c[1]] = [0]*(seq_length+1)
                    subst_read[c[0]+c[1]][a_pos] += 1

                try:
                    subst_seq[c[0]+c[1]][a_pos + mstart] += 1
                except:
                    subst_seq[c[0]+c[1]] = [0]*(len_ref+1)
                    subst_seq[c[0]+c[1]][a_pos + mstart] += 1
                    
                pos += 1
                
            for p in range(1,len(mpos)-1):
                if mpos[p-1] == mpos[p] - 1 and mpos[p] == mpos[p+1] - 1 : # two consecutive mismatches
                    fgh.write('%d\t%s\t%s\n' % (mpos[p-1], read, seq))
                    pot_gaps = True
                    break
        else:
            pass
        line = handle.readline()
    if pot_gaps:
        print >> sys.stderr, 'potential gaps found, reads and seq written to file'
    fgh.close()
    return [subst_read, subst_seq, read_start, coverage]
