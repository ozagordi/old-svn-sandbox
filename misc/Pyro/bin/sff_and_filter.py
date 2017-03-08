#!/usr/bin/env python

import os
import sys
import textwrap
import logging
import logging.handlers

import numpy as np
#import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt

import socket
hostname = socket.gethostname().split('.')[0]
if hostname != 'bs-mbp08':
    sys.path.pop(3)

import Bio
#assert Bio.__version__ == '1.54b'
from Bio import SeqIO

LOG_FILENAME = './sff_conversion.log'
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
file_out = args[2]
h_out = open(file_out, 'w')
try:
    phred_thresh = int(args[3])
except IndexError:
    phred_thresh = 10
logfun.info('Converting file %s into %s' % (file_in, file_out))
logfun.info('Discarding reads with one phred score < %d' % phred_thresh)

qname_in = os.path.basename(file_in)
dirname_in = os.path.dirname(file_in)

qname_out = os.path.basename(file_out)
dirname_out = os.path.dirname(file_out)

tot, good = 0, 0
read_l = []
min_pos = []
for record in SeqIO.parse(file_in, 'sff'):
    tot += 1
    left_clip =  record.annotations["clip_qual_left"]
    right_clip = record.annotations["clip_qual_right"]
    scores = record.letter_annotations['phred_quality'][left_clip:right_clip]
    min_score = min(scores)
    if min_score < phred_thresh:
        min_pos.append(scores.index(min_score))
    
    good += 1
    read_l.append(right_clip - left_clip + 1)
    h_out.write('>%s\n' % record.id)
    h_out.write(textwrap.fill(record.seq[left_clip:right_clip].tostring(), 80))
    h_out.write('\n')

min_pos = np.array(min_pos)

logfun.info('%d out of %d reads saved [%f %%]' % (good, tot, 100*float(good)/tot) )
logfun.info('read length is %d +/- %d' % (int(0.5+np.mean(read_l)), int(0.5+np.std(read_l))))

plt.title('Position of the first base of quality below %d' % phred_thresh)
plt.hist(min_pos, bins=50, color='orange')
plt.xlabel('position on the read')
plt.ylabel('reads')
plt.show()