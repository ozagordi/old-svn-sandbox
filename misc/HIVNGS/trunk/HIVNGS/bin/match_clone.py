#!/usr/bin/env python
'''This builds on my Alignment library to find the closest
matching haplotype in alignment file obtained with needleall
of clones to support file, e.g.
needleall viral_mix_5_seqs.fasta PR_9-support.fas -aformat=markx10 -ausashow3 -adesshow3
'''
__author__ = "Osvaldo Zagordi"
import sys
import os.path
homedir = os.path.expanduser('~/')
sys.path.append(os.path.join(homedir, 'pythonlib'))


def align_to_clones(s, clones_file):
    '''
    '''
    import subprocess
    import Alignment
    
    outfile = 'tmp.needle'
    needle_cline = 'needle asis:%s %s -aformat=markx10 -auto -out %s' % (s.seq, clones_file, outfile)  
    try:
        retcode = subprocess.call(str(needle_cline), shell=True)
        if retcode < 0:
            sys.exit('Child needle was terminated by signal %d' % -retcode)
#               else:
#                   print >> sys.stderr, 'Child needle returned %i' % retcode
    except OSError:
        sys.exit('Execution of needle failed: %s' % ee)
        pass
    
    alf = Alignment.alignfile2dict([outfile], 'support', 10.0, 0.5, Verbose = False)
    os.remove(outfile)
    mm = {}
    assert len(alf) == 1, len(alf)
    for seq, al in alf.items():
        assert len(al) == 5, len(al)
        for clone, res in al.items():
            res.summary()
            mm[clone] = res.mismatches
    
    return mm

def main():
    import re
    import pickle
    from Bio import SeqIO
    threshold = 0.9
    args = sys.argv
    # definition of files, directories and all that
    try:
        sup_file = args[1]
        clones_file = args[2]
    except IndexError:
        sys.exit('usage: %s sup_file clones_file' % args[0])
        
    freqs = {}
    for i, s in enumerate(SeqIO.parse(open(sup_file), 'fasta')):
        m_obj = re.search('posterior=(.*)\s*ave_reads=(.*)', s.description)
        post, ave_reads = map(float, (m_obj.group(1), m_obj.group(2)))
        print i, 'sequence'
        if post < threshold or ave_reads < 1.: continue
        if post > 1.0: print 'WARNING: posterior=', post
        
        align_res = align_to_clones(s, clones_file)
        assert align_res.values().count(0) < 2
        sorted_al_res = sorted(align_res, key=align_res.get)
        best_clone = sorted_al_res[0]
	best_mm = align_res[best_clone]
        if best_mm <= 0:
            freqs[best_clone] = freqs.get(best_clone, 0) + ave_reads
        else:
            freqs['unmatched'] = freqs.get('unmatched', 0) + ave_reads
        
    for k, v in freqs.items():
        print '| %8.4f | %s |' % (100*v/sum(freqs.values()), k)
if __name__ == '__main__':
    main()
    
