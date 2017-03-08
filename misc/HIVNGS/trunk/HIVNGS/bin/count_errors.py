#!/usr/bin/env python
'''This builds on my Alignment library to count errors
in reads coming from an experiment in which a single
clones was sequenced. Revision of just_error.py
'''
__author__ = "Osvaldo Zagordi"
import sys
import os.path
homedir = os.path.expanduser('~/')
sys.path.append(os.path.join(homedir, 'pythonlib'))


def parse_alignment(al_file_forw, al_file_rev):
    '''
    '''
    import Alignment
    
    alf = Alignment.alignfile2dict([al_file_forw], 'forward', 10.0, 0.5, Verbose = False)
    alr = Alignment.alignfile2dict([al_file_rev], 'reverse', 10.0, 0.5, Verbose = False)
    assert len(alf) == 1
    assert len(alr) == 1
    alif = alf.values()[0]
    alir = alr.values()[0]
    insertions = {}
    deletions = {}
    mismatches = {}
    bases = 0
    reads = 0
    for k in alif.keys():
        if alif[k].score > alir[k].score:
            v = alif[k]
        else:
            v = alir[k]  
        v.summary()
        if v.insertions + v.deletions > 10:
            continue
        
        insertions[v.insertions] = insertions.get(v.insertions, 0) + 1
        deletions[v.deletions] = deletions.get(v.deletions, 0) + 1
        mismatches[v.mismatches] = mismatches.get(v.mismatches, 0) + 1
        bases += v.stop - v.start + 1
        reads += 1
    data = insertions, deletions, mismatches, bases, reads
    
    return data

def main():
    import pickle
    
    args = sys.argv
    # definition of files, directories and all that
    try:
        forw_file = args[1]
        rev_file = args[2]
    except IndexError:
        sys.exit('usage: %s alignment_forw_file alignment_reverse_file' % args[0])
    file_only = os.path.split(forw_file)[1]
    
    pck_file = ''.join(file_only.split('.')[:-1]).split('_forw')[0] + '.pck'
    print pck_file
    if os.path.exists(pck_file):
        data = pickle.load(open(pck_file))
    else:
        oh = open(pck_file, 'w')
        data = parse_alignment(forw_file, rev_file)
        pickle.dump(data, oh)
        oh.close()
    
    insertions, deletions, mismatches, bases, reads = data
    mm = sum([k*v for k, v in mismatches.items()])
    ins = sum([k*v for k, v in insertions.items()])
    dels = sum([k*v for k, v in deletions.items()])
    print 'tot_reads\t', reads
    print 'tot_bases\t', bases
    print 'insertions\t', ins
    print 'deletions\t', dels
    print 'mismatch\t', mm
    print ''   
    print >> sys.stderr, 'Error rate per base (percentage)', 100*float(mm)/bases
    print >> sys.stderr, 'Ins per base (percentage)', 100*float(ins)/bases
    print >> sys.stderr, 'Dels per base (percentage)', 100*float(dels)/bases
    #plot_dist(mm, ig)


if __name__ == '__main__':
    main()
    