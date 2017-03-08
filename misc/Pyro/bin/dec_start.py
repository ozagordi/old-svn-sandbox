#!/usr/bin/env python

import os
import sys
import logging
import logging.handlers
import socket
from Bio import SeqIO


start, stop = 0, 3300
win_size, win_shifts = 400, 1

iterations, alpha = 5000, 0.001
our_hosts = ['bs-mbp08', 'bs-dsvr07', 'bs-dsvr24']
hostname = socket.gethostname().split('.')[0]

args = sys.argv

LOG_FILENAME = './dec_start.log'
# Make a global logging object.
x = logging.getLogger("logfun")
x.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler('./log', 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
x.addHandler(h)
logfun = logging.getLogger("logfun")

logfun.info(' '.join(args))
logfun.info('start=%d, stop=%d' % (start, stop))
logfun.info('iterations=%d, alpha=%f' % (iterations, alpha))

try:
    reads_file = os.path.realpath(args[1])
    ref_file = os.path.realpath(args[2])
    threshold = float(args[3])
except:
    sys.exit('usage: dec_start.py reads_file ref_file threshold')

if reads_file == '-h':
    sys.exit('usage: dec_start.py reads_file ref_file threshold')

assert hostname in our_hosts, 'unknown computer'

if hostname == 'bs-mbp08':
    shorah_dir = '/Users/ozagordi/sandbox/shorah/trunk/'
else:
    shorah_dir = '/nas/ozagordi/bin/shorah/trunk/'

diri_exe = os.path.join(shorah_dir, 'diri_sampler')
s2f_exe = os.path.join(shorah_dir, 's2f.py')
dec_exe = os.path.join(shorah_dir, 'dec.py')
assert  os.path.exists(diri_exe) and os.path.exists(dec_exe) and os.path.exists(s2f_exe), \
    '%s %s %s' % (diri_exe, dec_exe, s2f_exe)

base_dir = os.getcwd()

# wild_type.fas must me there
#assert os.path.exists(os.path.join(data_dir, 'wild_type.fas')), 'wild_type.fas does not exist'

try:
    os.remove('diri_sampler')
    os.link(diri_exe, 'diri_sampler')
except:
    os.link(diri_exe, 'diri_sampler')

logfun.info('\t# about to run\n')
# check if running alignment is needed
if os.path.exists('./reads.far'):
    logfun.info('reads.far exists, using it')
    pass
else:
    logfun.info('\t# running step2far')
    cmline = 'nice %s -r %s -f %s -o reads.far -t %f -k ' % \
        (s2f_exe, ref_file, reads_file, threshold)
    logfun.info(cmline + '\n')
    os.system(cmline)

#ys.exit()
# take only the center window

good_reads = ( record[start:stop] for record in SeqIO.parse(open('reads.far'), 'fasta') ) 
outfile = 'window.far'
out_handle = open(outfile, "w")
count = SeqIO.write(good_reads, out_handle, "fasta")
out_handle.close()

logfun.info('Taking %i reads overlapping with the given window' % count)

# finally run
logfun.info('\t# running dec.py')
cmline = 'nice %s -f window.far -w %d -s %d -j %d -a %f -k &> dec.log ' % (dec_exe, win_size, win_shifts, iterations, alpha)
logfun.info(cmline + '\n')
os.system(cmline)
logfun.info('\t# Exiting')

