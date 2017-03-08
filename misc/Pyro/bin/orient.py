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
import sys,os
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib
#from pythonlib import EmbossStandalone
from pythonlib.MarkxIO import Markx10Iterator
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
import subprocess
import time
import optparse
import random, math

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
#optparser.add_option("-k","--keep_files",help="dont delete temporary files <default=off>",action="store_true",dest="keep")
(options,args) = optparser.parse_args()

if options.o == None:
    f_output = sys.stdout
else:
    f_output = open(options.o, 'w')

# check the reference genome
fas_reference = options.ref
rh = open(fas_reference, 'rU')
ref = list(SeqIO.parse(rh, 'fasta') )
assert len(ref) == 1, 'One and only one sequence must be in the reference file'

gen_length = len(ref[0].seq)
# assert 'N' not in ref[0].seq, "Found an ambiguous position in the reference '%s'" % ref[0].id
print >> sys.stderr, 'The reference genome length is', gen_length

if options.input.split('.')[-1] in ['fas', 'fasta']:
    if os.path.isfile(options.input):
        f_fasta_filename = options.input
    else:
        print 'fasta file "%s" not found...'%options.input
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
f_fasta_reverse = open(f_fasta_reverse_filename, 'w')
SeqIO.write(seqlist,f_fasta_reverse,'fasta')
f_fasta_reverse.close()

print >> sys.stderr, 'Found', countreads, 'reads'

length = 0.
length2 = 0.
n =0.

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
    print >>sys.stderr, 'Aligning back...'
    cmline_forw = NeedleCommandline(asequence=options.ref, bsequence = f_fasta_forward_filename,
                            outfile='tmp_align_f.needle', gapopen=6.0, gapextend=3.0, aformat='markx10')
    child_process_forw = subprocess.call(str(cmline_forw), shell=True)

if not os.path.isfile('tmp_align_r.needle'):
    print >>sys.stderr, '...and forth'
    cmline_rev = NeedleCommandline(asequence=options.ref, bsequence = f_fasta_reverse_filename,
                            outfile='tmp_align_r.needle', gapopen=6.0, gapextend=3.0, aformat='markx10')
    child_process_rev = subprocess.call(str(cmline_rev), shell=True)

diff_ident = []
diff_score = []

f_ident = []
r_ident = []
f_score = []
r_score = []

counted_ident = []
good_score = []
good_ident = []

countreads_afterreverse = 0
count_forward = 0
count_reverse = 0

reads = {}
ref = {}
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
    
    f_ident.append(float(f_align._annotations['sw_ident']))
    r_ident.append(float(r_align._annotations['sw_ident']))
    f_score.append(float(f_align._annotations['sw_score']))
    r_score.append(float(r_align._annotations['sw_score']))
    
    basecount = 0.
    idencount = 0.
    accept_read = False
    
    if float(f_align._annotations['sw_score']) > float(r_align._annotations['sw_score']):
        tmp = f_align.get_seq_by_num(1).tostring().upper()
        refseq = f_align.get_seq_by_num(0).tostring().upper()
        count_forward += 1
        good_score.append(float(f_align._annotations['sw_score']))
        good_ident.append(float(f_align._annotations['sw_ident']))
    else:
        tmp = r_align.get_seq_by_num(1).tostring().upper()
        refseq = r_align.get_seq_by_num(0).tostring().upper()
        count_reverse += 1
        good_score.append(float(r_align._annotations['sw_score']))
        good_ident.append(float(r_align._annotations['sw_ident']))
    
    f_output.write('>%s\n' % f_align.get_all_seqs()[1].id)
    f_output.write('%s\n' % tmp.replace('-', ''))
