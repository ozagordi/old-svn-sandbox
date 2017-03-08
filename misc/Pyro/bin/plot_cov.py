#!/usr/bin/env python


import sys
import matplotlib.pyplot as plt
import numpy as np
covfile = sys.argv[1]
h = open(covfile)
pos, cov = np.loadtxt(covfile, dtype= ('i4', 'i4'), delimiter=' ', skiprows=1, unpack=True)
assert len(pos) == len(cov)

plt.title('Coverage for %s' % covfile.split('/')[-1])
plt.xlabel('position [bps]')
plt.ylabel('reads')
plt.yscale('log')
plt.plot(pos, cov, lw=2.5)

my_orange = '#E69F00'
plt.xticks()
ymin = plt.ylim()[0]
a = plt.gca()
a.set_xticklabels((), ())#[1000, 2000], ['a', 'b'])
for p in [2981, 3548, 6815, 7367]:
    plt.axvline(x=p, c=my_orange, lw=1.5, alpha=0.7)
    plt.text(p-25, ymin*1.05, str(p), rotation=90, verticalalignment='bottom', horizontalalignment='right')#, color=my_orange)

plt.show()
