#!/usr/bin/env python
__author__ = "Osvaldo Zagordi"
__version__ = "$Revision: 0.1 $"
__date__ = "$Date: 2009/07/129$"
__copyright__ = ""
__license__ = ""

import sys
import os.path
homedir = os.path.expanduser('~/')
sys.path.append(os.path.join(homedir, 'pythonlib'))

dna_code = {'A': set(['A']),
            'C': set(['C']),
            'G': set(['G']),
            'T': set(['T']),
            
            'R': set(['G', 'A']),
            'Y': set(['T', 'C']),
            'M': set(['A', 'C']),
            'K': set(['G', 'T']),
            'S': set(['G', 'C']),
            'W': set(['A', 'T']),
            
            'H': set(['A', 'C', 'T']),
            'B': set(['C', 'G', 'T']),
            'V': set(['A', 'C', 'G']),
            'D': set(['A', 'G', 'T']),
            'N': set(['A', 'C', 'G', 'T']),
            '-': set(['A', 'C', 'G', 'T'])
}

class AlignDict(dict):
    '''
    A set of alignment instances. Subclass dict
    in order to have a 2D dictionary.
    '''
    
    def __init__(self, set_name, gap_open, gap_extend, default = dict):
        '''
        Initialization needs a name, gap open, gap extension
        '''
        self.name = set_name
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.default = default
        
    def __getitem__(self, key):
        '''
        subclassing dictionary getitem
        '''
        if key not in self:
            self[key] = {} #self.default()
        return dict.__getitem__(self, key)


class AlignInstance:
    '''
    Represents a single pairwise alignment
    '''
    def __init__(self, id_a, id_b, seq_a, seq_b, score, descr_a=None, descr_b=None):
        '''
        Initialization needs id of both sequences, zipped sequence
        pair, and score
        ''' 
        self.id_a = id_a
        self.id_b = id_b
        
        self.score = score
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.descr_a = descr_a
        self.descr_b = descr_b
        
    def summary(self):
        '''
        summary must be called to calculate start, stop, internal gaps, and identity
        '''
        start = None
        stop = None
        i = 0
        import sys
        import itertools

        it_pair = itertools.izip(self.seq_a, self.seq_b)
        while True:
            i += 1
            try:
                p = it_pair.next()
            except:
                break
            if p is None:
                break
            if start == None:
                if '-' not in p:
                    start = i
            else:
                if '-' not in p:
                    stop = i
        
        self.start = start
        self.stop = stop
        
        self.insertions = 0 #gap in the first sequence
        self.deletions = 0 #gap in the second sequence 
        self.ident = 0
        self.mismatches = 0
        
        it_pair = itertools.izip(self.seq_a[start-1:stop], self.seq_b[start-1:stop])
        while True:
            i += 1
            try:
                p = it_pair.next()
            except:
                break
            if p is None:
                break
            if p == ('-', '-'):
                print >> sys.stderr, ' double gap in %s %s' % (self.id_a, self.id_b)
            if p [0] == '-':
                self.insertions += 1
            elif p[1] == '-':
                self.deletions += 1
            if p[0].upper() == p[1].upper():
                self.ident += 1
            if dna_code[p[0].upper()] & dna_code[p[1].upper()] == set([]):
                self.mismatches += 1
                
        return


    def print_info(self):
        '''
        Print all information to stderr
        '''
        import sys
        print >> sys.stderr, 'id_a:', self.id_a
        print >> sys.stderr, 'id_b:', self.id_b
        print >> sys.stderr, 'start:', self.start
        print >> sys.stderr, 'stop:', self.stop
        print >> sys.stderr, 'score:', self.score
        print >> sys.stderr, 'mismatches:', self.mismatches
        print >> sys.stderr, 'ident:', self.ident

def alstart(seq_a, seq_b):
    import tempfile
    import subprocess
    import os
    
    f = tempfile.NamedTemporaryFile()#delete=False)
    fn = f.name
    f.close()
    needle_align(seq_a, seq_b, fn)
    af = alignfile2dict([fn], 'noname', 6.0, 3.0)
    ak1 = af.keys()[0]
    ak2= af[ak1].keys()[0]
    a = af[ak1][ak2]
    a.summary()
    os.remove(fn)
    return a.start


def alignfile2dict(al_files, name, gap_open, gap_extend, Verbose = False): #, r_file):
    """
    Takes a list of files and returns an alignment 2D dictionary
    """
    import sys
    import MarkxIO
    
    al_set = AlignDict(name, gap_open, gap_extend)
    for f_file in al_files:
        f_forward = open(f_file)
        forwardaligniter = MarkxIO.Markx10Iterator(f_forward)
        if Verbose:
            print >> sys.stderr, 'parsing the alignments', f_file
        while True:
            try:
                f_align = forwardaligniter.next()
            except:
                break
            if f_align is None:
                break
            score = float(f_align._annotations['sw_score'])
            # id_a, id_b = f_align.get_all_seqs()[0].id.split('#')[0],\
                # f_align.get_all_seqs()[1].id.split('#')[0]
            id_a, id_b = list(f_align)[0].id.split('#')[0],\
                list(f_align)[1].id.split('#')[0]
            descr_a = list(f_align)[0].description
            descr_b = list(f_align)[1].description
            a, b = f_align.get_seq_by_num(0).tostring().upper(),\
                f_align.get_seq_by_num(1).tostring().upper()
            al_set[id_a][id_b] = AlignInstance(id_a, id_b, a, b, score, descr_a, descr_b)
            
        f_forward.close()
    return al_set


def needle_align(a_seq, b_seq, out_file, go = 6.0, ge = 3.0, reverse1=False, reverse2=False):
    """
    """
#    from Bio import Application
    from pythonlib import cmline # import NeedleCommandline
    import subprocess
    import os
    import sys
    Verbose = False
    needle_run = cmline.NeedleCommandline()
    needle_run.set_parameter('-asequence', a_seq)
    needle_run.set_parameter('-bsequence', b_seq)
    needle_run.set_parameter('-outfile', out_file)
    
    needle_run.set_parameter('-gapopen', go)
    needle_run.set_parameter('-gapextend', ge)
    needle_run.set_parameter('-aformat', 'markx10')
    needle_run.set_parameter('-auto')
    if reverse1:
            needle_run.set_parameter('-sreverse1')
    if reverse2:
            needle_run.set_parameter('-sreverse2')

    line_1 = str(needle_run)
    
    #print >> sys.stderr, line_1, 'running'
    if os.path.exists(out_file):
        return out_file, 'exists'
    try:
        retcode = subprocess.call(line_1, shell=True)
        if retcode < 0:
            if Verbose:
                print >>sys.stderr, "'%s'" % line_1
                print >>sys.stderr, "Child diri_sampler was terminated by signal", -retcode
        else:
            if Verbose:
                print >>sys.stdout, "Child diri_sampler returned %i" % retcode
    except OSError, ee:
        if Verbose:
            print >>sys.stderr, "Execution of diri_sampler failed:", ee

    return a_seq, b_seq, 'aligned'

#def novo_align(seq_file, index_file, out_file, Verbose=False):
def novo_align(arg1):
    '''
    '''
    import subprocess
    Verbose = False
    seq_file, index_file, out_file = arg1
    
    line_1 = 'novoalign -f %s -d %s > %s' % (seq_file, index_file, out_file)
    try:
        retcode = subprocess.call(line_1, shell=True)
        if retcode < 0 and Verbose:
            print >>sys.stderr, "'%s'" % line_1
            print >>sys.stderr, "Child diri_sampler was terminated by signal", -retcode
        else:
            if Verbose:
                print >>sys.stdout, "Child diri_sampler returned %i" % retcode
    except OSError, ee:
        print >>sys.stderr, "Execution of diri_sampler failed:", ee
    
    return

def main():
    '''
    Quick test
    '''
    print ' running a quick test '.center(60, '#')
    print ' some information is printed '.center(60, '#')
    A = AlignSet('First_set', 6.0, 3.0)
    a = 'AAAAACACACGCAT-------'
    b = '---AACAC-C-CATTTTTTTT'
    print a
    print b
    ai = AlignInstance('a', 'b', a, b, 120.0)
    ai.summary()
    ai.print_info()
    A['a']['b'] = ai
    print A['a']['b'].ident
    print A['a'].values()[0].score

if __name__ == '__main__':
    main()
