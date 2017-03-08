#!/usr/bin/env python

import os
import sys
import textwrap
import logging
import logging.handlers

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

import socket
hostname = socket.gethostname().split('.')[0]
if hostname != 'bs-mbp08':
    sys.path.pop(3)

import Bio
#assert Bio.__version__ == '1.54b'
from Bio import SeqIO


def scoreseq(ref, seq):
    from subprocess import Popen, PIPE
    
    com = "needle asis:%s asis:%s -gapopen 6.0 -gapextend 3.0 -stdout -auto | grep Score" % (ref, seq)
    pipe = Popen(com, shell='/bin/bash', bufsize=1024, stdout=PIPE, close_fds=True).communicate()[0]
    scoreplus = float(pipe.strip().split()[-1])
    
    com = "needle asis:%s asis:%s -sreverse2 -gapopen 6.0 -gapextend 3.0 -stdout -auto | grep Score" % (ref, seq)
    pipe = Popen(com, shell='/bin/bash', bufsize=1024, stdout=PIPE, close_fds=True).communicate()[0]
    scoreminus = float(pipe.strip().split()[-1])

    return scoreplus, scoreminus
    


LOG_FILENAME = './sff_qual_histo.log'
# Make a global logging object
x = logging.getLogger("logfun")
x.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler('./log', 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s\t%(asctime)s\n\t%(message)s\n")
h.setFormatter(f)
x.addHandler(h)
logfun = logging.getLogger("logfun")

# command line and save it to log
args = sys.argv
logfun.info(' '.join(args))

file_in = args[1]
assert file_in.endswith('.sff')
file_stem = file_in.split('.')[-2]
file_ref = args[2]
try:
    phred_thresh = int(args[3])
except IndexError:
    phred_thresh = 10
#logfun.info('Converting file %s into %s' % (file_in, file_out))
#logfun.info('Discarding reads with one phred score < %d' % phred_thresh)
refseq = list(SeqIO.parse(open(file_ref), 'fasta'))[0].seq.tostring()
qname_in = os.path.basename(file_in)
dirname_in = os.path.dirname(file_in)

tot, good = 0, 0
read_l = []
min_pos_plus = []
min_pos_minus = []
for record in SeqIO.parse(file_in, 'sff'):
    tot += 1
    left_clip =  record.annotations["clip_qual_left"]
    right_clip = record.annotations["clip_qual_right"]
    scores = record.letter_annotations['phred_quality'][left_clip:right_clip]
    recstr = record.seq[left_clip:right_clip].tostring()
    min_score = min(scores)
    
    if min_score < phred_thresh:
        scoreplus, scoreminus = scoreseq(recstr, refseq)
        if scoreplus > scoreminus:
            min_pos_plus.append(scores.index(min_score))
        else:
            min_pos_minus.append(scores.index(min_score))
        
    good += 1
    read_l.append(right_clip - left_clip + 1)
    if tot % 1000 == 0:
        logfun.info('Now checked %d reads' % tot)    
min_pos_plus = np.array(min_pos_plus)
min_pos_minus = np.array(min_pos_minus)
print min_pos_plus
logfun.info('%d out of %d reads saved [%f %%]' % (good, tot, 100*float(good)/tot) )
logfun.info('read length is %d +/- %d' % (int(0.5+np.mean(read_l)), int(0.5+np.std(read_l))))

plt.title('Position of the first base of quality below %d' % phred_thresh)
plt.hist(min_pos_plus, bins=50, color='orange', alpha=0.7, label='forward')
plt.hist(min_pos_minus, bins=50, color='blue', alpha=0.7, label='reverse')
plt.xlabel('position on the read')
plt.ylabel('reads')
plt.legend()
filename = 'first_low_qual_%s' % file_stem
imtype = 'pdf'
plt.savefig('%s.%s' % (filename, imtype), dpi=None, facecolor='w', edgecolor='w',\
            orientation='landscape', papertype=None, format=imtype,\
            transparent=False)
#plt.show()
