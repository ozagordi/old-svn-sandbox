#!/usr/bin/env python

import os
import sys

from Bio import SeqIO

import matplotlib
#matplotlib.use('pdf')
import matplotlib.pyplot as plt
import numpy as np
import gzip

quantity1 = 'untouched'
window = 500
filename = sys.argv[1]
if filename.endswith('.gz'):
    dtype = [('', np.int32)]*3 + [('', np.float32)]*2
    f = gzip.open(filename)
    y = np.loadtxt(f, dtype=dtype, skiprows=1)
    y = y.view(np.dtype([('iter', np.int32), ('K', np.int32), ('untouched', np.int32),
                         ('theta', np.float32), ('gamma', np.float32)]))

else:
    with open(filename) as f:
        r = f.readline().split()
        dtype = [('', np.int32)]*3 + [('', np.float32)]*2
        y = np.loadtxt(f, dtype=dtype, skiprows=0)
        y = y.view(np.dtype([('iter', np.int32), ('K', np.int32), ('untouched', np.int32),
                             ('theta', np.float32), ('gamma', np.float32)]))


#variance
variance = np.array([ y[quantity1][i:i+window].var()/y[quantity1][i:i+window].mean() for i in range(0, len(y)-window, 100 ) ])

# Set figure scale, font and similia
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.semilogy(variance, 'bv-')
ax1.set_xlabel('iterations')
ax1.set_ylabel('variance')
#ax1.set_ylim(1E-3)

plt.show()
