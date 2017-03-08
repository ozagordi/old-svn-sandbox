#!/usr/bin/env python
# Copyright 2007, 2008, 2009, 2010
# Niko Beerenwinkel,
# Arnab Bhattacharya,
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
import os
blue = '#2171b5'

def plot_posterior(freq_file, sup_dict, threshold=0.2):
    ''' Plot posterior and frequency boxplot for files freq_file.csv
    and support produced by diri_sampler
    '''
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import numpy as np
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

    if freq_file.endswith('.gz'):
        with mygzipfile(freq_file) as f:
            r = f.readline().split()[2:]
            history = len(r)
            dtype = [('haplotypes', 'S16')] + [('', np.int32)]*(history+1)
            y = np.loadtxt(f, dtype=dtype)
            y = y.view(np.dtype([('haplotypes', 'S16'), ('support', np.int32, 1),
                                 ('reads', np.int32, history)]))
    else:
        with open(freq_file) as f:
            r = f.readline().split()[2:]
            history = len(r)
            dtype = [('haplotypes', 'S16')] + [('', np.int32)]*(history+1)
            y = np.loadtxt(f, dtype=dtype)
            y = y.view(np.dtype([('haplotypes', 'S16'), ('support', np.int32, 1),
                                 ('reads', np.int32, history)]))
            
    haplotypes = y['haplotypes']
    support = y['support']
    reads = y['reads']
    
    assert history == len(reads[0]), 'Some iteration is missing'
    
    cov = np.sum(reads, axis=0)
    coverage = cov[0]
    chk = np.equal(cov/cov[0], np.ones(history))
    if not np.all(chk):
        print 'Some read is missing in %s' % freq_file
        
    n_supported=0
    sup_values = []
    tmp_reads = []
    for i, h in enumerate(haplotypes):
        this_hap = sup_dict['hap_%d' % i]
        if support[i] >= threshold*history:# and reads[i].mean() > 1.0 and a.int_gaps <= 2:
            tmp_reads.append(reads[i,:])
            sup_values.append(float(support[i])/history)
            n_supported += 1
            
    print >> sys.stderr, 'Supported: ',n_supported
    x = np.arange(1, n_supported+1)
    reads_t = np.array(tmp_reads).transpose().astype(float)/coverage
    
    print >> sys.stderr, 'Now the figure'
    
    # Set figure scale, font and similia
    fig = plt.figure(figsize=(12, 9))
    # fig.set_yscale('log', basey=10)
    fig.hold(True)
    
    # first plot
    print >> sys.stderr, 'First plot'
    ax1 = plt.subplot(111)
    # plt.grid('True')
    plt.subplots_adjust(left = 0.09, right = 0.93, top = 0.98, bottom = 0.05,\
                            hspace = 0.01)
    ax1.set_yscale('log', basey=10)
    d = plt.boxplot(reads_t, positions=x, sym='', hold=True)
    plt.setp(d['boxes'], lw=4, color='r')
    plt.setp(d['medians'], lw=0, color='black')
    plt.setp(d['caps'], lw=0) #color='black)
    plt.setp(d['whiskers'], lw=0) #color='black', lw=4)
    ax1.set_xlim(0.5, n_supported+0.5)
    plt.yticks(color='r')
    plt.ylabel('frequency', color='r')
    ax1.xaxis.grid(True, which='major')
    plt.setp( plt.gca().get_xticklabels(), rotation=90)#, horizontalalignment='right', size='x-small')

    # second plot
    print >> sys.stderr, 'Second plot'
    ax2 = plt.twinx()
    # plt.grid('True')
    plt.plot(x, sup_values, 'o-', ms=12, c=blue) #ls='-', c='r', 'bo')
    ax2.set_xlim(0.5, n_supported+0.5)
    mm = min(0.8, 0.9*min(sup_values))
    ax2.set_ylim(ymin=mm, ymax = 1.05)
    ax2.set_xticks(x)
    
    plt.yticks(color=blue)
    plt.ylabel('posterior', color=blue)
    font = {'size'   : 18,
        # 'family' : 'monospace',
        'weight' : 'bold'
        }
    plt.rc('font', **font)  # pass in the font dict as kwargs
    matplotlib.rcParams['xtick.major.size'] = 6
    matplotlib.rcParams['ytick.major.size'] = 6

    # save the figure
    imtype = 'pdf'
    plt.savefig('box_freq_support.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='portrait', papertype=None, format=imtype,\
                    transparent=False)

def main():

    import re
    import gzip
    from Bio import SeqIO
    # import subprocess
    # rule_reads = re.compile('ave_reads=(.*)')
    # rule_post = re.compile('posterior=(.*) ave')
    rule_key = re.compile('w(\d*)-(\d*)')
    
    args = sys.argv
    try:
        run_dir = args[1]
        win_number = int(args[2])
    except:
        sys.exit('usage: %s run_directory window_number' % args[0].split('/')[-1])
    sup_dir = os.path.join(run_dir, 'support')
    freq_dir = os.path.join(run_dir, 'freq')
    
    freq_keys = dict([(int(rule_key.match(f).group(1)), f) for f in os.listdir(freq_dir)] )
    sup_keys =  dict([(int(rule_key.match(f).group(1)), f) for f in os.listdir(sup_dir)] )
    
    for k, v in freq_keys.items():
        assert freq_keys[k][:6] == sup_keys[k][:6], v
    
    s_freq_keys = sorted(freq_keys)
    s_sup_keys =  sorted(sup_keys)
    freq_file = os.path.join(freq_dir, freq_keys[s_freq_keys[win_number-1]])
    sup_file = os.path.join(sup_dir, sup_keys[s_sup_keys[win_number-1]])
    print >> sys.stderr, 'Doing', freq_keys[s_freq_keys[win_number-1]]
    print sup_file
    
    if sup_file.endswith('.gz'):
        h = gzip.open(sup_file)
    else:
        h = open(sup_file)
    sup = SeqIO.parse(h, 'fasta')
    
    sup_dict = dict([(s.id.split('|')[0], s.seq.tostring().replace('-', '').strip('N')) for s in sup])
    
    plot_posterior(freq_file, sup_dict)

if __name__ == '__main__':
    main()
