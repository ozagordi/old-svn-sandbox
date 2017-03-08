#!/usr/bin/env python

import sys
import os

homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib
from Bio import SeqIO


#PI wild type
PIWT = {
    10: 'L',
    11: 'V',
    16: 'G',
    20: 'K',
    24: 'L',
    30: 'D',
    32: 'V',
    33: 'L',
    34: 'E',
    36: 'M',
    46: 'M',
    47: 'I',
    48: 'G',
    50: 'I',
    53: 'F',
    54: 'I',
    60: 'D',
    62: 'I',
    63: 'L',
    64: 'I',
    71: 'A',
    73: 'G',
    76: 'L',
    77: 'V',
    82: 'V',
    84: 'I',
    85: 'I',
    88: 'N',
    89: 'L',
    90: 'L',
    93: 'I'
    }

# PIs drug resistant mutation
"""
PIRM = {
    'ataza': [10, 16, 20, 24, 32, 33, 34, 36, 46, 48,
              50, 53, 54, 60, 62, 64, 71, 73, 82, 84,
              85, 88, 90, 93],
    
    'fosam': [10, 32, 46, 47, 50, 54, 73, 76, 82, 84,
              90],

    'darun': [11, 32, 33, 47, 50, 54, 73, 76, 84, 89],

    'indin': [10, 20, 24, 32, 36, 46, 54, 71, 73, 76,
              77, 82, 84, 90],

    'lopin': [10, 20, 24, 32, 33, 46, 47, 50, 53, 54,
              63, 71, 73, 76, 82, 84, 90],

    'nelfi': [10, 30, 36, 46, 71, 77, 82, 84, 88, 90],

    'saqui': [10, 24, 48, 54, 62, 71, 73, 77, 82, 84,
              90],

    'tipra': [10, 13, 20, 33, 35, 36, 43, 46, 47, 54,
              58, 69, 74, 82, 83, 84, 90]
}
"""
PIRM = {
    'ataza': { 10: 'IFVC',
               16: 'E',
               20: 'RMITV',
               24: 'I',
               32: 'I',
               33: 'IFV',
               34: 'Q',
               36: 'ILV',
               46: 'IL',
               48: 'V',
               50: 'L',
               53: 'LY',
               54: 'LVMTA',
               60: 'E',
               62: 'V',
               64: 'LMV',
               71: 'VITL',
               73: 'CSTA',
               82: 'ATFI',
               84: 'V',
               85: 'V',
               88: 'S',
               90: 'M',
               93: 'LM' },
    
    'fosam': { 10: 'FIRV',
               32: 'I',
               46: 'IL',
               47: 'V',
               50: 'V',
               54: 'LVM',
               73: 'S',
               76: 'V',
               82: 'AFST',
               84: 'V',
               90: 'M' },

    'darun': { 11: 'I',
               32: 'I',
               33: 'F',
               47: 'V',
               50: 'V',
               54: 'ML',
               73: 'S',
               76: 'V',
               84: 'V',
               89: 'V' },

    'indin': { 10: 'IRV',
               20: 'MR',
               24: 'I',
               32: 'I',
               36: 'I',
               46: 'IL',
               54: 'V',
               71: 'VT',
               73: 'SA',
               76: 'V',
               77: 'I',
               82: 'AFT',
               84: 'V',
               90: 'M' },

    'lopin': { 10: 'FIRV',
               20: 'MR',
               24: 'I',
               32: 'I',
               33: 'F',
               46: 'IL',
               47: 'VA',
               50: 'V',
               53: 'L',
               54: 'VLAMTS',
               63: 'P',
               71: 'VT',
               73: 'S',
               76: 'V',
               82: 'AFTS',
               84: 'V',
               90: 'M' },

    'nelfi': { 10: 'FI',
               30: 'N',
               36: 'I',
               46: 'IL',
               71: 'VT',
               77: 'I',
               82: 'AFTS',
               84: 'V',
               88: 'DS',
               90: 'M' },

    'saqui': { 10: 'IRV',
               24: 'I',
               48: 'V',
               54: 'VL',
               62: 'V',
               71: 'VT',
               73: 'S',
               77: 'I',
               82: 'AFTS',
               84: 'V',
               90: 'M' },

    'tipra': { 10: 'V',
               13: 'V',
               20: 'MR',
               33: 'F',
               35: 'G',
               36: 'I',
               43: 'T',
               46: 'L',
               47: 'V',
               54: 'AMV',              
               58: 'E',
               69: 'K',
               74: 'P',
               82: 'LT',
               83: 'D',
               84: 'V',
               90: 'M' }
    }

drug_names = {
    'ataza': 'Atazanavir',
    'fosam': 'Fosamprenavir',
    'darun': 'Darunavir',
    'indin': 'Indinavir',
    'lopin': 'Lopinavir',
    'nelfi': 'Nelfinavir',
    'saqui': 'Saquinavir',
    'tipra': 'Tipranavir'
    }

def counts(self, reverse=False):
    """ return list of keys sorted by value
    """
    aux = [ (self[k], k) for k in self ]
    aux.sort()
    if reverse: aux.reverse()
    return [k for v, k in aux]

def drm(file, drug):
    """
    """
    from Bio import SeqIO
    
    h = open(file)
    plst = list(SeqIO.parse(h, 'fasta'))
    print >> sys.stderr, 'Found', len(plst), 'reads'
    mutations = {}
    for c in PIRM[drug]:
        mutations[c] = {}
    sites = PIRM[drug]
    
    
    for prot in plst:
        letters = [ prot.seq[s-1] for s in sites ]
        for l in zip(sites, letters):
            try:
                mutations[l[0]][l[1]] += 1
            except KeyError:
                mutations[l[0]][l[1]] = 1.
    
    for k in sorted(mutations):
        del(mutations[k]['-'])
#        print k, mutations[k]

    return mutations



def plot_drm(file, drug, plot_log=True):
    """
    """
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
        import scipy.stats as stats
        import math
    except:
        print 'exit, could not import matplotlib'
        sys.exit()
    
    ca = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow']
    cb = ['#00BFFF', # alternative colour set
          '#BC8F8F',
          '#ADFF2F',
          '#B8860B',
          '#00FF7F',
          '#DC143C',
          '#A0522D',
          '#8A2BE2',
          '#D8BFD8',
          '#87CEEB',
          '#6495ED']
    cn = len(ca)
    
    acol = {
        'G': 'green',
        'S': 'green',
        'T': 'green',
        'Y': 'green',
        'C': 'green',
        'N': 'green',
        'A': 'orange',
        'V': 'orange',
        'L': 'orange',
        'I': 'orange',
        'P': 'orange',
        'W': 'orange',
        'F': 'orange',
        'M': 'orange',
        'D': 'red',
        'E': 'red',
        'K': 'blue',
        'R': 'blue',
        'H': 'blue',
        'Q': 'brown',
        '*': 'black'}
    
    ax = plt.subplot(111)
    ax.set_yscale('log')
    mutations = drm(file, drug)
    xx = (8, 92)
    yy = (1, 1)
    plt.scatter(xx, yy, alpha=.0, label='')

    # sanger limit
    thresh = 0.20
    plt.axhline(thresh, color='blue')
    note = '%d' % int(100*thresh) + '%'
    #   plt.yticks([thresh], note)
    plt.annotate(note, (1*xx[1], thresh*1.1), color='blue')
    plt.text(0.98, 0.97, 'wild-type', horizontalalignment='right', verticalalignment='center', family='normal',
             color='black', transform = ax.transAxes)
    plt.text(0.98, 0.935, 'resistant', horizontalalignment='right', verticalalignment='center', family='normal',
             color='red', transform = ax.transAxes)
    plt.text(0.98, 0.9, 'unknown', horizontalalignment='right', verticalalignment='center', family='normal',
             color='green', transform = ax.transAxes)
    
    font = {'family' : 'monospace',
            'weight' : 'bold'
            # 'size'   : 10
            }
    plt.rc('font', **font)  # pass in the font dict as kwargs

    for k in sorted(mutations):
        
        # vertical bar
        plt.axvline(x=int(k), linewidth = 8.0, alpha=0.1, color='grey', label='')
        
        this = mutations[k]
        ss = sum(this.values())
        
        aux = [int(k) for c in this]
        xs = np.array(aux)
        
        aux = [ this[c]/ss for c in this ]
        ys = np.array(aux)
        
        plt.scatter(xs,ys, alpha=0.0)
        
        for c in this:
            xy = (int(k)-1, 0.9*this[c]/ss)
            if c in PIRM[drug][k]:
                col = 'red'
            elif c == PIWT[k]:
                col = 'black'
            else:
                col = 'green'
            """
            print this, k, drug,PIRM_ann[drug][k][0]
            print c, col
            """
            plt.annotate(c, xy, color=col, label='')
            
    font_t = {
        'weight' : 'normal',
        'size'   : 14
        }
    plt.xticks(sorted(mutations), size=10, weight='semibold', family='normal')
    plt.yticks(weight='semibold', family='normal', size=10)
    plt.ylabel('frequency', weight='semibold', family='normal', size=10)
    plt.title('Variation on %s resistance sites' % drug_names[drug], **font_t)
    filename = file.split('.')[0] + '_' + drug
    imtype = 'pdf'
    plt.savefig('%s.%s' % (filename, imtype), dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()



def main():
    """
    """
    
    args = sys.argv
    
    try:
        file = args[1]
        drug = args[2]
    except:
        print >> sys.stderr, 'usage: plot_drm.py file drug'
        print >> sys.stderr, 'drug names are:'
        for d in PIRM:
            print >> sys.stderr, '\t\t', d
        sys.exit()
    
    plot_drm(file, drug)

if __name__ == '__main__':
    main()
