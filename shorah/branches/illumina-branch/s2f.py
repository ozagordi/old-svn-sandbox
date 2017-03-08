#!/usr/bin/env python

# Copyright 2007, 2008
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
# ETH Zurich

# This file is part of ShoRAH.
# ShoRAH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ShoRAH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os

#homedir = os.path.expanduser('~/')
#sys.path.append(homedir)
import pythonlib
from pythonlib import EmbossStandalone

from Bio import SeqIO

from pythonlib.MarkxIO import Markx10Iterator

import time
import optparse
import random, math

needle_exe = 'needle'
acclength = 2.
amb_thresh = 2
min_length = 1#00
dna_code = ['A','C','T','G']
alphabet = ['A', 'C', 'G', 'T', 'N']

f_code = {}
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


optparser = optparse.OptionParser()
optparser.add_option("-f","--readfile",type="string",default="",dest="input")
optparser.add_option("-r","--ref",type="string",default="",dest="ref")
optparser.add_option("-o","--output",help="output suffix must be '.far' or none for stdout",type="string",dest="o")
optparser.add_option("-t","--threshold",help="threshold of similarity for throwing reads away... <default=0.7>",type="float",dest="threshold",default=0.7)
optparser.add_option("-k","--keep_files",help="dont delete temporary files <default=off>",action="store_true",dest="keep")
(options,args) = optparser.parse_args()

if options.o == None:
    f_output = sys.stdout
else:
    if(options.o.split('.')[-1] != 'far'):
        sys.exit('The suffix of output file must be .far')
    f_output = open(options.o, 'w')

try:
    keep_files = options.keep
except:
    keep_files = False


# check the reference genome
fas_reference = options.ref
rh = open(fas_reference, 'rU')
ref = list(SeqIO.parse(rh, 'fasta') )
assert len(ref) == 1, 'One and only one sequence must be in the reference file'
gen_length = len(ref[0].seq)
if gen_length > 20:
    gen_ref_start = ref[0].seq.tostring()[:20]
else:
    gen_ref_start = ref[0].seq.tostring()

assert 'N' not in ref[0].seq, "Found an ambiguous position in the reference '%s'" % ref[0].id
print >> sys.stderr, 'The reference genome length is', gen_length

if options.input.split('.')[-1] in ['fas', 'fasta']:
    if os.path.isfile(options.input):
        f_fasta_filename = options.input
    else:
        print 'fasta file "%s" not found...' % options.input
        sys.exit()
else:
    print 'format of input-file not supported...'
    sys.exit()


f_fasta = open(f_fasta_filename)
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

length = 0.
length2 = 0.
n = 0.

for read in seqlist:
    if read.seq.count('N') < amb_thresh:
        len_seq = len(read.seq)
        length += float(len_seq)
        length2 += float(len_seq*len_seq)
#        readdict[read.] = [seq,len_seq]
        n += 1.
            
meanlr = length/n
stdlr = math.sqrt((n*length2 - length*length)/(n*n - n))
allowed_length = [meanlr - acclength * stdlr, meanlr + (1+acclength) * stdlr]
print >> sys.stderr, 'Allowed interval for length is', allowed_length


if not os.path.isfile('tmp_align_f.needle'):
    print >> sys.stderr, 'needle forward'
    EmbossStandalone.needle(needle_exe, options.ref, f_fasta_forward_filename,
                            out='tmp_align_f.needle', gapopen=6.0, gapext=3.0, aglobal3='False')
"""
else:
    print >>sys.stderr, 'The alignment file tmp_align_f.needle is already present'
    statinfo = os.stat('tmp_align_f.needle')
    age_sec = time.time() - statinfo.st_mtime
    if age_sec > 3600:
        print >>sys.stderr, 'Warning: it was modified more than an hour ago'
    age = time.gmtime(age_sec)
    
    print >>sys.stderr, 'If you want to run the alignment again, remove it'
    print >>sys.stderr, "using existing 'tmp_align_f.needle'..."
"""
if not os.path.isfile('tmp_align_r.needle'):
    print >> sys.stderr, 'needle backward'
    EmbossStandalone.needle(needle_exe, options.ref, f_fasta_reverse_filename,
                            out='tmp_align_r.needle', gapopen=6.0, gapext=3.0, aglobal3='False')
"""
else:
    print >>sys.stderr, 'The alignment file tmp_align_r.needle is already present'
    statinfo = os.stat('tmp_align_f.needle')
    age_sec = time.time() - statinfo.st_mtime
    if age_sec > 3600:
        print >>sys.stderr, 'Warning: it was modified more than an hour ago'
    age = time.gmtime(age_sec)
    print >>sys.stderr, 'If you want to run the alignment again, remove it'
"""

#print dir(normalalign.next())

counted_ident = []

countreads_afterreverse = 0
count_forward = 0
count_reverse = 0

reads = {}
ref = {}
global_ins = {}
insertions = {}
ID = 0
maxlen = [0, 0]
disc_seq = []

discarded_fn = 'discarded.fas'
disc_h = open(discarded_fn, 'w')

f_forward = open('tmp_align_f.needle')
f_reverse = open('tmp_align_r.needle')

forwardaligniter = Markx10Iterator(f_forward)
reversealigniter = Markx10Iterator(f_reverse)
shorts = 0
ref_start = 100000
ref_end = 0
pos = 0

print >> sys.stderr, 'parsing the alignments'
# iterates through the alignments
while True:
    
#    pos += 1
#    print >> sys.stderr,  '\x1B[1A\x1B[2K', pos    
    try:
        f_align = forwardaligniter.next()
        r_align = reversealigniter.next()
    except:
        break
    if f_align is None or r_align is None:
        break
    
    assert f_align.get_all_seqs()[1].id == r_align.get_all_seqs()[1].id, 'same seq back and forward'
    
    basecount = 0.
    idencount = 0.
    accept_read = False
    
    if float(f_align._annotations['sw_score']) > float(r_align._annotations['sw_score']):
        tmp = f_align.get_seq_by_num(1).tostring().upper()
        refseq = f_align.get_seq_by_num(0).tostring().upper()
        count_forward += 1
    else:
        tmp = r_align.get_seq_by_num(1).tostring().upper()
        refseq = r_align.get_seq_by_num(0).tostring().upper()
        count_reverse += 1
    
    # get start and stop values for the read
    mins = []
    for b in alphabet:
        if b in tmp:
            mins.append(tmp.find(b))
    q_align_start = min(mins)
    maxvalues = [ tmp.rfind(c) for c in alphabet ]
    q_align_end = max(maxvalues)
    
    # get start and stop values for the matched reference
    mins = []
    for b in alphabet:
        if b in refseq:
            mins.append(refseq.find(b))
    m_align_start = min(mins)
    maxvalues = [ refseq.rfind(c) for c in alphabet ]
    m_align_end = max(maxvalues)
    
    align_start = max(m_align_start, q_align_start)
    align_end = min(m_align_end, q_align_end)
    
    if align_start < ref_start:
        ref_start = align_start
    if align_end > ref_end:
        ref_end = align_end
    
    tmpl = list(tmp)
    
    # take only the overlap with the reference
    # if the read begins before the ref_genome...
    if m_align_start > 0:
        for i in range(m_align_start):
            tmpl[i] = '-'
            
    # ...or if it continues also after it
    if q_align_end > m_align_end:
        tmpl = tmpl[:m_align_end+1]
        refseq = refseq[:m_align_end+1]
    
    tmp = ''.join(tmpl)
    
    assert len(tmp) == len(refseq),  'Read and match lengths are %d %d' % (len(tmp), len(refseq))
    # alignment must always start from the beginning of the reference
    assert refseq.replace('-','').startswith(gen_ref_start), 'not considering the beginning'
    # last character of the matched reference cannot be a gap
    assert not refseq.endswith('-'), 'last character of the matched reference cannot be a gap\nrefseq:\n%s' % refseq
    
    #detect positions of insertions
    pos_ins = {}
    i = 0
    count = 0
    
    while i <= align_end:
        if refseq[i] == '-':
            # there's an insertion, how long?
            ext = len(refseq[i:]) - len(refseq[i:].lstrip('-'))
            pos_ins[count] = ext
            i += ext-1
        elif refseq[i] != '-':
            count += 1
        i += 1
        
    
    for i in range(align_start, align_end+1):
        if tmp[i] == refseq[i]:
            idencount += 1.
        if tmp[i] in dna_code:
            basecount += 1.
    try:
        counted_ident.append(idencount/basecount)
    except:
        print tmp
        print refseq
        sys.exit('Read %s is empty' %  f_align.get_all_seqs()[1].id)
    
    if options.threshold < idencount/basecount:
        reads[ID] = tmp
        ref[ID] = refseq
        insertions[ID] = pos_ins
        
        for pi in pos_ins:
            try:
                if global_ins[pi] < pos_ins[pi]:
                    global_ins[pi] = pos_ins[pi]
            except:
                global_ins[pi] = pos_ins[pi]
        
        if len(refseq) > maxlen[0]:
            maxlen[0] = len(refseq)
            maxlen[1] = ID
        ID += 1
        countreads_afterreverse += 1
    else:
        pass
    """
    ts = r_align.get_seq_by_num(1).tostring()
    tid =  f_align.get_all_seqs()[1].id
    disc_seq.append(SeqRecord(Seq(ts, generic_nucleotide), id=tid, description=''))
    tseq = f_align.get_all_seqs()[1]
    disc_h.write('>%s\n' % tseq.id)
    disc_h.write(tseq.seq.tostring().replace('-', '') + '\n')
    """
disc_h.close()

print >>sys.stderr,'%d reads were above the threshold (discarded), %d reads left' % (countreads - countreads_afterreverse,countreads_afterreverse)
try:
    print >>sys.stderr,'dropped %d reads with wron length' % sff_droppedreads_length
except:
    pass
print >>sys.stderr,'forward: %d, reverse: %d' % (count_forward,count_reverse)

#print >>sys.stderr, ref_end

assert len(reads) == len(ref)

################################
# propagating the gaps         #
# from the pairwise alignments #
################################
print >>sys.stderr, 'propagating gaps'
ids = ref.keys()

new_reads = {}
begin = {}
pad = '-'
fasta_length = 80

for ID in ids:
    to_insert = {}
    
    for ti in global_ins:
        try:
            already = insertions[ID][ti]
        except:
            already = 0
        to_insert[ti] = global_ins[ti] - already
    
    i = 0
    this_read = []
    start = False
    padded = 0
    
    for c in zip(reads[ID], ref[ID]):
        
        if i in to_insert and c[1] != '-':
            n = to_insert[i]
            this_read.append(pad*n)
            padded += n
        this_read.append(c[0])
        
        if c[1] != '-':
            i += 1
        
        if start == False and c[0] != '-':
            start = True
            begin[ID] = i
        
    assert padded == sum(to_insert.values()), 'ID = %d\tpadded = %d\t sum = %d' % (ID, padded, sum(to_insert.values()))
    
    new_reads[ID] = ''.join(this_read)
    

print >> sys.stderr, 'sorting'
items = [ (v, k) for k, v in begin.items() ]
items.sort()

for i in items:
    ID = i[1]
    f_output.write('>read_%d\n' % ID)
    cc = 0
    for c in new_reads[ID]:
        f_output.write(c)
        cc += 1
        if cc % fasta_length == 0:
            f_output.write('\n')
    if cc % fasta_length != 0:
            f_output.write('\n')

if not keep_files:
    try:
        os.remove('tmp_reads_f.fas')
    except:
        pass
    try:
        os.remove('tmp_reads_r.fas')
    except:
        pass
    try:
        os.remove('tmp_align_f.needle')
    except:
        pass
    try:
        os.remove('tmp_align_r.needle')
    except:
        pass
