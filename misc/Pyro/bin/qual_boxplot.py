#!/usr/bin/env python

import sys
from Bio import SeqIO
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np

step = 20
filename = sys.argv[1]

qual = [ record.letter_annotations["phred_quality"] for record in SeqIO.parse(open(filename), 'qual') if len(record.seq) >= 200 ]
max_len = 200 #max([len(rq) for rq in qual])

print 'Found', len(qual), 'reads'
print 'Max len is', max_len

av_q = {}
for pos in np.arange(0, max_len, step):
    av_q[pos+1] = [ np.mean(qv[pos:pos+step]) for qv in qual]

# Set figure scale, font and similia
fig = plt.gca()
fig.hold(True)
positions = []
data = []
for k, v in av_q.iteritems():
    positions.append(k+step-1)
    data.append(v)
plt.boxplot(data, notch=0, sym='o', vert=1, whis=1.5, positions=positions, widths=7)
#plt.setp( plt.gca().get_xticklabels(), rotation=90, horizontalalignment='right', size='x-small')
plt.xlim(step-10, 10+max(positions))
plt.ylim(10, 42)
plt.xlabel('bp')
plt.ylabel('phred score')
font = {'size'   : 15,
        #        'family' : 'monospace',
        #        'weight' : 'bold'
        }
plt.rc('font', **font)  # pass in the font dict as kwargs
#plt.show()

# save the figure

imtype = 'pdf'
plt.savefig('box_freq.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                orientation='portrait', papertype=None, format=imtype,\
                transparent=False)

