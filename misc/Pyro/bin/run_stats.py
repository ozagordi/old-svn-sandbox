#!/usr/bin/env python

# Copyright 2007, 2008, 2009
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
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

import sys
import gzip

class mygzipfile(gzip.GzipFile):
    ''' Subclassing GzipFile such that supports __enter__ and __exit__
    and it can be used in with_as statement
    '''
    def __enter__(self):
        if self.fileobj is None:
            raise ValueError("I/O operation on closed GzipFile object")
        return self
    
    def __exit__(self, *args):
        self.close()

def plot_stats(filenames, quantity1, quantity2=None):
    ''' Plot quantity for ds-out file
    '''
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import numpy as np
    last = -1000

    run_name = filenames[0].split('/')[-3]
    
    y = np.array([])
    y = y.view(np.dtype([('iter', np.int32), ('K', np.int32), ('untouched', np.int32), \
                               ('theta', np.float32), ('gamma', np.float32)]))
    
    for tf in filenames:
        with mygzipfile(tf) as f:
            r = f.readline().split()
            dtype = [('', np.int32)]*3 + [('',np.float32)]*2
            ty = np.loadtxt(f, dtype=dtype)
            ty = ty.view(np.dtype([('iter', np.int32), ('K', np.int32), ('untouched', np.int32),
                             ('theta', np.float32), ('gamma', np.float32)]))
        y = np.append(y, ty[last:])
        '''
        iter = y['iter']
        K = y['K']
        untouch = y['untouch']
        theta = y['theta']
        gamma = y['gamma']
        '''
    print y.shape
    data_1 = y[quantity1]
    try:
        data_2 = y[quantity2]
    except:
        pass
    print data_1.mean()
    # Set figure scale, font and similia
    fig = plt.figure()
    plt.title('Analysis of the run, file %s' % run_name)
    if quantity2:
        ax1 = fig.add_subplot(211)
    else:
        ax1 = fig.add_subplot(111)
        
    ax1.hist(y[quantity1], bins=100)#, 'bv-')

    # Make the y-axis label and tick labels match the line color.
    q1l = quantity1
    if quantity1 == 'gamma':
        q1l = r'$\Gamma$'
    if quantity1 == 'theta':
        q1l = r'$\theta$'
        
    ax1.set_xlabel(q1l, color='b')
    
    if quantity2:
        ax2 = fig.add_subplot(212)
        ax2.hist(y[quantity2], bins=100)#, 'rx-')
        q2l = quantity2
        if quantity2 == 'gamma':
            q2l = r'$\Gamma$'
        if quantity2 == 'theta':
            q2l = r'$\theta$'
        ax2.set_xlabel(q2l, color='r')

        
    # save the figure

    imtype = 'pdf'
    if quantity2:
        plt.savefig('run_%s_%s_%s.%s' % (run_name, quantity1, quantity2, imtype), dpi=None, facecolor='w', edgecolor='w',\
                        orientation='portrait', papertype=None, format=imtype,\
                        transparent=False)
    else:
        plt.savefig('run_%s_%s.%s' % (run_name, quantity1, imtype), dpi=None, facecolor='w', edgecolor='w',\
                        orientation='portrait', papertype=None, format=imtype,\
                        transparent=False)

def main():
    ''' What does a main do?
    '''
    import sys
    import os
    
    args = sys.argv
    poss_quant = ['theta', 'gamma', 'K', 'untouched']
    try:
        dirname = args[1]
        quantity1 = args[2]
    except:
        sys.exit('usage: plot_sampling.py files quantity1 [optional: quantity2] [%s]' % '|'.join(poss_quant))
        
    try:
        quantity2 = args[3]
        if quantity2 not in poss_quant:
            sys.exit('usage: plot_sampling.py files quantity1 [optional: quantity2] [%s]' % '|'.join(poss_quant))
    except:
        quantity2 = None
        
    if quantity1 not in poss_quant:
        sys.exit('usage: plot_sampling.py files quantity1 [optional: quantity2] [%s]' % '|'.join(poss_quant))
        
    sam_dir = os.path.join(dirname, 'sampling')
    filenames = [os.path.join(sam_dir, f) for f in os.listdir(sam_dir) if f.startswith('w1-540') or f.startswith('w181-720')]
    print filenames
    plot_stats(filenames, quantity1, quantity2)
    
if __name__ == "__main__":
    main()
