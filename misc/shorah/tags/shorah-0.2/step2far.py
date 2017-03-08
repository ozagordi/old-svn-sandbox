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


from pythonlib import SFFParser, EmbossStandalone
from Bio import SeqIO
import MyAlignIO
import sys,os
import time
import optparse
import random,math

needle_exe = './needle'
acclength = 3.
amb_thresh = 2
dna_code = ['A','C','T','G']

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
optparser.add_option("-t","--threshold",help="threshold of errors for throwing reads away... <default=0.4>",type="float",dest="threshold",default=0.4)
optparser.add_option("-k","--keep_files",help="dont delete temporary files <default=off>",action="store_true",dest="keep")
(options,args) = optparser.parse_args()

if options.o == None:
    f_output = sys.stdout
else:
    if(options.o.split('.')[-1] != 'far'):
        sys.exit('The suffix of output file must be .far')
    f_output = open(options.o,'w')

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
assert 'N' not in ref[0].seq, "Found an ambiguous position in the reference '%s'" % ref[0].id
print 'The reference genome length is', gen_length


sff_droppedreads_quality = 0

if options.input.split('.')[-1] == 'sff':
    sff_reader = SFFParser.SFFReader(options.input)
    reads = sff_reader.reads

    f_tmp = open('tmp_reads_f.fas','w')
    print 'converting sff-file to fasta format...'

    length = 0.
    length2 = 0.
    n =0.

    readdict = {}

    for read in reads:
        seq = ''.join(read.bases).strip('N')
        if seq.count('N') < amb_thresh:
            len_seq = len(seq)
            length += float(len_seq)
            length2 += float(len_seq*len_seq)
            readdict[read.name] = [seq,len_seq]
            n += 1.
            
    meanlr = length/n
    stdlr = math.sqrt((n*length2 - length*length)/(n*n - n))
    allowed_length = [meanlr - acclength * stdlr, meanlr + acclength * stdlr]


    for ID in readdict.iterkeys():

        """
        print 'bases:\t',read.bases
        print 'clip_adapter:\t',read.clip_adapter_left,read.clip_adapter_right
        print 'clip_qual:\t',read.clip_qual_left,read.clip_qual_right
        print 'eight_byte_padding:\t',read.eight_byte_padding
        print 'flow_index_per_base:\t',read.flow_index_per_base
        print 'flowgram_values:\t',read.flowgram_values
        print 'mydataoffset:\t',read.mydataoffset
        print 'mytype:\t',read.mytype
        print 'name:\t',read.name
        print 'name_length:\t',read.name_length
        print 'number_of_bases:\t',read.number_of_bases
        print 'quality_scores:\t',read.quality_scores
        print 'read_header_length:\t',read.read_header_length
        sys.exit(1)
        """
        if readdict[ID][1] > allowed_length[0] and readdict[ID][1] < allowed_length[1]:
            s = ''
            print >>f_tmp,'>%s_%d'%(ID,readdict[ID][1])
            for letter in readdict[ID][0]:
                if len(s) == 80:
                    print >>f_tmp,s
                    s = ''
                if letter in dna_code:
                    s += letter
                elif f_code.has_key(letter):
                    s += random.choice(f_code[letter])
            if len(s) != 0:
                print >>f_tmp,s
            else:
                sff_droppedreads_length += 1

    f_tmp.close()
    f_fasta_filename = 'tmp_reads_f.fas'
elif options.input.split('.')[-1] in ['fas','fasta']:
    if os.path.isfile(options.input):
        f_fasta_filename = options.input
    else:
        print 'fasta file "%s" not found...'%options.input
        sys.exit()
else:
    print 'format of input-file not supported...'
    sys.exit()


f_fasta = open(f_fasta_filename)
seqlist = list(SeqIO.parse(f_fasta,'fasta'))
for seq in seqlist:
    seq.seq = seq.seq.reverse_complement()
f_fasta.close()
countreads = len(seqlist)

f_fasta_reverse_filename = 'tmp_reads_r.fas'
f_fasta_reverse = open(f_fasta_reverse_filename,'w')
SeqIO.write(seqlist,f_fasta_reverse,'fasta')
f_fasta_reverse.close()

if not os.path.isfile('tmp_align_f.needle'):
    EmbossStandalone.needle(needle_exe,options.ref,f_fasta_filename,out='tmp_align_f.needle',gapopen=6.0,gapext=3.0,aglobal3='False')
else:
    print >>sys.stderr, 'The alignment file tmp_align_f.needle is already present'
    statinfo = os.stat('tmp_align_f.needle')
    age_sec = time.time() - statinfo.st_mtime
    if age_sec > 3600:
        print >>sys.stderr, 'Warning: it was modified more than an hour ago'
    age = time.gmtime(age_sec)
    
    print >>sys.stderr, 'If you want to run the alignment again, remove it'
    print >>sys.stderr, "using existing 'tmp_align_f.needle'..."

if not os.path.isfile('tmp_align_r.needle'):
    EmbossStandalone.needle(needle_exe,options.ref,f_fasta_reverse_filename,out='tmp_align_r.needle',gapopen=6.0,gapext=3.0,aglobal3='False')
else:
    print >>sys.stderr, 'The alignment file tmp_align_r.needle is already present'
    statinfo = os.stat('tmp_align_f.needle')
    age_sec = time.time() - statinfo.st_mtime
    if age_sec > 3600:
        print >>sys.stderr, 'Warning: it was modified more than an hour ago'
    age = time.gmtime(age_sec)
    print >>sys.stderr, 'If you want to run the alignment again, remove it'


f_normal = open('tmp_align_f.needle')
normalaligniter = iter(MyAlignIO.parse(f_normal,'markx10'))

#print dir(normalalign.next())

diff_ident = []
diff_score = []

countreads_afterreverse = 0
count_forward = 0
count_reverse = 0

reads = {}
ref = {}
ID = 0
maxlen = [0,0]
f_reverse = open('tmp_align_r.needle')
for reversealign in MyAlignIO.parse(f_reverse,'markx10'):
    normalalign = normalaligniter.next()

    assert normalalign.get_all_seqs()[1].id == reversealign.get_all_seqs()[1].id, 'same seq back and forward'

    #print normalalign.get_all_seqs()[1].id
    #print reversealign.get_all_seqs()[1].id
    #print normalalign._annotations['sw_score']
    #print reversealign._annotations['sw_score']
    #print normalalign._annotations['sw_ident']
    #print reversealign._annotations['sw_ident']
    #print normalalign.get_seq_by_num(1)
    #print reversealign.get_seq_by_num(1)
    # diff_ident.append(float(normalalign._annotations['sw_ident']) - float(reversealign._annotations['sw_ident']))
    # diff_score.append(float(normalalign._annotations['sw_score']) - float(reversealign._annotations['sw_score']))
    basecount = 0.
    errorcount = 0.
    accept_read = False
    if float(normalalign._annotations['sw_ident']) > float(reversealign._annotations['sw_ident']):
        tmp = normalalign.get_seq_by_num(1).tostring().upper()
        refseq = normalalign.get_seq_by_num(0).tostring().upper()
        count_forward += 1
    else:
        tmp = reversealign.get_seq_by_num(1).tostring().upper()
        refseq = reversealign.get_seq_by_num(0).tostring().upper()
        count_reverse += 1
    align_start = min([tmp.find('A'),tmp.find('C'),tmp.find('T'),tmp.find('G')])
    align_end = max([tmp.rfind('A'),tmp.rfind('C'),tmp.rfind('T'),tmp.rfind('G')])
    for i in range(align_start,align_end):
        if tmp[i] != refseq[i]:
            errorcount += 1.
        if tmp[i] in dna_code:
            basecount += 1.
    if options.threshold > errorcount/basecount:
        if len(refseq) > maxlen[0]:
            maxlen[0] = len(refseq)
            maxlen[1] = ID
        reads[ID] = tmp
        ref[ID] = refseq
        ID += 1
        
        #tmp = tmp.replace('-','').strip('N')
        #print >>f_output,">%s"%normalalign.get_all_seqs()[1].id.split('#')[0]
        #print >>f_output,tmp
        countreads_afterreverse += 1


print >>sys.stderr,'%d reads were above the threshold, %d reads left'%(countreads - countreads_afterreverse,countreads_afterreverse)
try:
    print >>sys.stderr,'dropped %d reads with wron length'%sff_droppedreads_length
except:
    pass
print >>sys.stderr,'forward: %d, reverse: %d'%(count_forward,count_reverse)


#assert len(reads) == len(ref)


################################
# propagating the gaps         #
# from the pairwise alignments #
################################

pos = 0
while pos < maxlen[0]:
    # list of ref-sequences with no gaps in it at current pos
    gaplist = ref.keys()
    removed_from_gaplist = False
    # list of ref-sequences which already terminated
    shortend = []
    for id in ref.iterkeys():
        if len(ref[id]) > maxlen[0]:
            maxlen[0] = len(ref[id])
            maxlen[1] = id
        try:
            if ref[id][pos] == '-':
                gaplist.remove(id)
                removed_from_gaplist = True
        except IndexError:
            shortend.append(id)
            gaplist.remove(id)

        # debugging output:
        #print "\t%d\t%d\t%d\t%d\t%d\t%d"%(pos,gl_size_start,len(gaplist),gl_size_start-len(gaplist),len(shortend),max[0])


    # only when there are gaps
    if removed_from_gaplist:
        if maxlen[1] in gaplist:
            maxlen[0] += 1
        # go through the list and insert gaps in ref and read
        for id in gaplist:
            ref[id] = ref[id][:pos] + '-' + ref[id][pos:]
            reads[id] = reads[id][:pos] + '-' + reads[id][pos:]
    # insert gaps in ref and read, when already terminated
    for id in shortend:
        ref[id] = ref[id] + '-'
        reads[id] = reads[id] + '-'

    pos += 1


# sort according to startpos

startpos = {}
for ID in reads.iterkeys():
    start_align = min([reads[ID].find('A'),reads[ID].find('C'),reads[ID].find('T'),reads[ID].find('G')])
    if not startpos.has_key(start_align):
        startpos[start_align] = []
    startpos[start_align].append(ID)


for i in range(maxlen[0]):
    if startpos.has_key(i):
        for ID in startpos[i]:
            print >> f_output, ">read#%s"%ID
            s = ''
            for letter in reads[ID]:
                if len(s) == 80:
                    print >> f_output,s
                    s = ''
                s += letter
            if len(s) > 0:
                print >> f_output,s




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
