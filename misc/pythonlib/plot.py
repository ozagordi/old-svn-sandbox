#!/usr/bin/env python
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

args = sys.argv
try:
    filename = args[1]
except:
    sys.exit('usage: %s filename' % args[0].split('/')[-1])
try:
    cols = [int(abs(int(i)))-1 for i in args[2:]]
except:
    cols = [0, 1]

if len(cols) == 2:
    y2col = None
    x, y = np.loadtxt(filename, unpack=True, usecols=cols)
elif len(cols) == 3:
    y2col = cols[-1]
    x, y, y2 = np.loadtxt(filename, unpack=True, usecols=cols)
else:
    sys.exit('too many columns')

signs = [1, 1, 1]
for i, a in enumerate(args[2:]):
    try:
        signs[i] = int(a)/abs(int(a))
    except:
        pass
print 'Using columns: ', [i+1 for i in cols]
x *= signs[0]
y *= signs[1]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x, y)
for label in ax1.get_xticklabels(): 
    label.set_rotation(30)
    
if y2col:
    y2 *= signs[2]
    # Make the y-axis label and tick labels match the line color.
    color_1 = ax1.get_lines()[0].get_color()
    #ax1.set_ylabel('exp', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color(color_1)
        
    ax2 = ax1.twinx()
    ax2.plot(x, y2, color='r')
    #ax2.set_ylabel('sin', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')

plt.show()