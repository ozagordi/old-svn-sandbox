#!/usr/bin/env python

import os
import sys

from Bio import SeqIO

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

args = sys.argv

file = args[1]

try:
    phred_thresh = int(args[2])
except:
    phred_thresh = 20

good_ids = ( record.id for record in \
                 SeqIO.parse(open(file), 'qual') \
                 if min(record.letter_annotations["phred_quality"]) >= phred_thresh  and len(record.seq) > 200)

qname = os.path.basename(file)
dirname = os.path.dirname(file)
#fname = qname.split('.')[0] + '.fna'
fname = 'patient_2/reads.fas'
handle = open(os.path.join(dirname, fname))
print >> sys.stderr, 'original reads are in', handle.name

record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
print >> sys.stderr, "Found %d reads" % len(record_dict)

good_reads = (record_dict[k] for k in good_ids if k in record_dict.keys())
outfile = 'ppp' #os.path.join('filtered_data_%d' % phred_thresh, fname.split('/')[0]+'.fas')
out_handle = open(outfile, "w")
count = SeqIO.write(good_reads, out_handle, "fasta")
out_handle.close()

print >> sys.stderr, "Saved %i reads" % count

print 1.0*count/len(record_dict)
