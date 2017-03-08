#!/usr/bin/env python

import sys
import os
import re
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib

Verbose = False
Pedantic = True

'''
ref_genomes = ['07-56681',
               '07-56951',
               '08-59712',
               '07-54825',
               '08-04134',
               '08-01315',
               '08-55163',
               '08-57881',
               '08-02659',
               '08-04512']
'''

# geometric mixture
geo_freq = [0.5]
for i in range(1,10):
    geo_freq.append(geo_freq[i-1] * 0.5)

def plot_pie_hap(freq_dict, haplo_file=None):
    """
    
    """
    import operator
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
    except:
        print 'exit, could not import matplotlib'
        sys.exit()

    import math
    
    cb = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow']
    ca = ['#00BFFF', # alternative colour set
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
    xt_lab = []
    xt_val = []
    plt.figure(figsize=(8, 8))
    ax = plt.axes([0.15, 0.1, 0.83, 0.83])
    
    print >> sys.stderr, 'plotting result'
    assert len(freq_dict) == 1, 'only one window here please'
        
    wind = freq_dict.keys()[0]
    d = sorted(freq_dict[wind].iteritems(), key=operator.itemgetter(1), reverse=True)
    lab = [i[0] for i in d]
    data = [i[1] for i in d]
    values = [math.log(i, 2) for i in data]
    geo_values = [math.log(i, 2) for i in geo_freq]
    x1 = np.arange(1, len(values)+1)
    x2 = np.arange(1, len(geo_values)+1)
    
    ass_file = '/'.join(haplo_file.split('/')[0:-1]) + '/assign.txt'
    actual = np.array([0.0]*len(geo_values))
    for line in open(ass_file):
        r, h = map(int, line.strip().split())
        actual[h] += 1.0
    log_act = [math.log(i/sum(actual), 2) for i in actual]
    
    plt.plot(x2, geo_values, '-', alpha=0.3, label='geometric series', color='r')
    plt.plot(x2, log_act, '^', label='present', markersize=14, alpha=0.7, color='c')
    plt.plot(x1, values, 'o', label='estimated', markersize=14, alpha=0.7, color='m')
    
    plt.axis([0.5, 10.5, -10.5, -0.5])
    plt.xticks(x2)
    plt.yticks(np.arange(-1, -11, -1))
    plt.grid(alpha = 0.3)
    plt.ylabel(r'frequency [$\log_2$]')
    plt.xlabel('haplotype index')
    plt.legend(loc='upper right', numpoints=1)
    
    fontsize = 22
    params = {'font.size': fontsize}
    plt.rcParams.update(params)

    #    plt.xticks('')

    imtype = 'pdf'
    plt.savefig('geo_plot.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='portrait', papertype=None, format=imtype,\
                    transparent=False)
    
    plt.show()
    
    
    return

###########################




def needle_align(a_file, b_file, out_file):
    """

    """
    from Bio.Application import generic_run
    from cmline import NeedleCommandline
    go = 6.0
    ge = 3.0
    
    
    needle_run = NeedleCommandline()
    needle_run.set_parameter('-asequence', a_file)
    needle_run.set_parameter('-bsequence', b_file)
    needle_run.set_parameter('-gapopen', go)
    needle_run.set_parameter('-gapextend', ge)
    needle_run.set_parameter('-outfile', out_file)
    needle_run.set_parameter('-aformat', 'markx10')
    result, messages, errors = generic_run(needle_run)
    
    if Pedantic:
        for ar in result.available_results():
            print ar, result.get_result(ar)
            
        for m in messages.readlines():
            print >>sys.stderr, m

        for e in errors.readlines():
            print >>sys.stderr, e

    return


def reads2clones_align(sample_dir, max_freq=2):
    """
    
    """
    from Bio.Application import generic_run
    from cmline import NeedleCommandline
    go = 6.0
    ge = 3.0
    
    if max_freq >10:
        sys.exit('10 reference genomes maximum')
    
    # align uncorrected reads
    
    for f in range(max_freq):
        ref_genome = './raw/Fastas/clone' + ref_genomes[f] + '.fsta'
        outfile = '%s/reads-%s.needle' % ( sample_dir, ref_genomes[f] )
        needle_run = NeedleCommandline()
        needle_run.set_parameter('-asequence', ref_genome)
        needle_run.set_parameter('-bsequence', '%s/reads.fas' % sample_dir)
        needle_run.set_parameter('-gapopen', go)
        needle_run.set_parameter('-gapextend', ge)
        needle_run.set_parameter('-outfile', outfile)
        needle_run.set_parameter('-aformat', 'markx10')
        result, messages, errors = generic_run(needle_run)
        
        for ar in result.available_results():
            print ar, result.get_result(ar)
            
        if Verbose:
            for m in messages.readlines():
                print >>sys.stderr, m

            for e in errors.readlines():
                print >>sys.stderr, e
 
   # align corrected reads

    for f in range(max_freq):
        ref_genome = './raw/Fastas/clone' + ref_genomes[f] + '.fsta'
        outfile = '%s/reads-cor-%s.needle' % ( sample_dir, ref_genomes[f] ) 
        needle_run = NeedleCommandline()
        needle_run.set_parameter('-asequence', ref_genome)
        needle_run.set_parameter('-bsequence', '%s/reads.cor.fas' % sample_dir)
        needle_run.set_parameter('-gapopen', go)
        needle_run.set_parameter('-gapextend', ge)
        needle_run.set_parameter('-outfile', outfile)
        needle_run.set_parameter('-aformat', 'markx10')
        result, messages, errors = generic_run(needle_run)
        
        for ar in result.available_results():
            print ar, result.get_result(ar)
            
        if Verbose:
            for m in messages.readlines():
                print >>sys.stderr, m

            for e in errors.readlines():
                print >>sys.stderr, e

    return


def align_info(sample_dir, max_freq=2):
    """
    
    """
    from Bio import SeqIO
    from Bio import SeqUtils
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from MarkxIO import Markx10Iterator
    
    verbose = False
    thresh = 0.03
    
    al_files = [ '%s/reads-%s.needle' % ( sample_dir, ref_genomes[f])  for f in range(max_freq) ]
    al_cor_files = [ '%s/reads-cor-%s.needle' % ( sample_dir, ref_genomes[f])  for f in range(max_freq) ]
    
    print al_files, al_cor_files
    # ########################
    # get info for uncorrected
    mm_raw = [[]]
    
    for f in range(max_freq):
        print 'opening', al_files[f]
        handle = open(al_files[f])
        mm_raw[f] = []
        mm_raw.append([])
        
        for al in Markx10Iterator(handle):
            refseq = al.get_seq_by_num(0).tostring().upper()
            tmp = al.get_seq_by_num(1).tostring().upper()
            
            # WARNING: not reliable on short sequences (a base might be absent, giving min=0)
            align_start = min([tmp.find('A'),tmp.find('C'),tmp.find('T'),tmp.find('G')])
            align_end = max([tmp.rfind('A'),tmp.rfind('C'),tmp.rfind('T'),tmp.rfind('G')])
            mm = [1 for c in zip(refseq[align_start:align_end], tmp[align_start:align_end]) if c[0] != c[1] ]
            
            mm_raw[f].append( float(len(mm)) / ( align_end - align_start + 1 ) )


    # ########################
    # get info for corrected
    mm_cor = [[]]
    
    for f in range(max_freq):
        print 'opening', al_cor_files[f]
        handle = open(al_cor_files[f])
        mm_cor[f] = []
        mm_cor.append([])
        
        for al in Markx10Iterator(handle):
            refseq = al.get_seq_by_num(0).tostring().upper()
            tmp = al.get_seq_by_num(1).tostring().upper()
            
            # WARNING: not reliable on short sequences (a base might be absent, giving min=0)
            align_start = min([tmp.find('A'),tmp.find('C'),tmp.find('T'),tmp.find('G')])
            align_end = max([tmp.rfind('A'),tmp.rfind('C'),tmp.rfind('T'),tmp.rfind('G')])
            
            mm = [1 for c in zip(refseq[align_start:align_end], tmp[align_start:align_end]) if c[0] != c[1] ]
            mm_cor[f].append( float(len(mm)) / ( align_end - align_start + 1 ) )
    
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
    except:
        print 'exit'
        sys.exit()
    
    
    plt.figure(1, figsize=(10,12)) #, dpi=100)
    ax = plt.subplot(411, axisbg='lightgrey')
    plt.subplots_adjust(hspace=0.15)
    freq = 0
    
    nraw, bins, patches = plt.hist(mm_raw[freq], 20, normed=0, facecolor='green', alpha=0.5, label='raw reads')
    ncor, bins, patches = plt.hist(mm_cor[freq], 20, normed=0, facecolor='red', alpha=0.5, label='corrected reads')
    plt.title('%d raw, %d corrected reads' % (len(mm_raw[freq]), len(mm_cor[freq])))
    aaa = plt.axis(xmin=0, xmax=0.55)

    ax = plt.subplot(412, axisbg='lightgrey')
    freq = 1
    
    nraw, bins, patches = plt.hist(mm_raw[freq], 20, normed=0, facecolor='green', alpha=0.5, label='raw reads')
    ncor, bins, patches = plt.hist(mm_cor[freq], 20, normed=0, facecolor='red', alpha=0.5, label='corrected reads')
    aaa = plt.axis(xmin=0, xmax=0.55)

    ax = plt.subplot(413, axisbg='lightgrey')
    freq = 2
    
    nraw, bins, patches = plt.hist(mm_raw[freq], 20, normed=0, facecolor='green', alpha=0.5, label='raw reads')
    ncor, bins, patches = plt.hist(mm_cor[freq], 20, normed=0, facecolor='red', alpha=0.5, label='corrected reads')
    aaa = plt.axis(xmin=0, xmax=0.55)

    ax = plt.subplot(414, axisbg='lightgrey')
    freq = 3
    
    nraw, bins, patches = plt.hist(mm_raw[freq], 20, normed=0, facecolor='green', alpha=0.5, label='raw reads')
    ncor, bins, patches = plt.hist(mm_cor[freq], 20, normed=0, facecolor='red', alpha=0.5, label='corrected reads')
    aaa = plt.axis(xmin=0, xmax=0.55)
    
    plt.xlabel('distance')
    plt.legend()

    imtype = 'pdf'
    plt.savefig('corr.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()
    return

def freq_window(haplo_file):
    """
    
    """
    
    import pickle
    import heapq
    import operator
    from Bio import SeqIO
    from Bio import SeqUtils
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from MarkxIO import Markx10Iterator
    
    verbose = False
    thresh = 0.03
    freq_dict = {}
    hap_key = haplo_file.split('/')[-1].rstrip('haplo.out')
    freq_dict[hap_key] = {}
    
    h = open(haplo_file)
    hap_seq = list(SeqIO.parse(h, 'fasta'))
    count = 0
    for hap_item in hap_seq:
        oh = open('tmp_1hap.fas', 'w')
        oh.write('>%s\n' % hap_item.id)
        oh.write(hap_item.seq.tostring().replace('-', ''))
        oh.write('\n')
        oh.close()
        needle_align('tmp_1hap.fas', 'seq.fas', 'tmp_al.needle')
        hal = open('tmp_al.needle')
        min_dist = 1.0
        this_freq = float(hap_item.id.split('|')[1])
        al_records = list(Markx10Iterator(hal))
        assert len(al_records) == 10, 'not all records'
        distances = {}
        for al in al_records:
            assert len(al.get_all_seqs()) == 2, '2 sequences per alignment'
            tmp = al.get_seq_by_num(0).tostring().upper()
            clone = al.get_seq_by_num(1).tostring().upper()
            query, match = list(al)
            # WARNING: not reliable on short sequences (a base might be absent, giving min=0)
            align_start = min([tmp.find('A'),tmp.find('C'),tmp.find('T'),tmp.find('G')])
            align_end = max([tmp.rfind('A'),tmp.rfind('C'),tmp.rfind('T'),tmp.rfind('G')])
            mm = [ 1 for c in zip(clone[align_start:align_end], tmp[align_start:align_end]) if c[0] != c[1] ]
            this_dist = float(len(mm)) / (align_end - align_start + 1)
            assert sum(mm) == len(mm), 'sum equals length'
            assert 0.0 <= this_dist and this_dist <= 1.0, 'not a distance'
            distances[match.id.split('(')[0]] = this_dist
            


        best2 = heapq.nsmallest(2, distances.iteritems(), operator.itemgetter(1))
        s1, s2 = best2[0][1], best2[1][1]
        td = (s2-s1)
        
        if s1 == 0.0 and td > 0.005:
            best_clone = best2[0][0]
            freq_dict[hap_key][best_clone] = freq_dict[hap_key].get(best_clone, 0) + this_freq
            count += 1
            
    wh = open('%s.dict' % haplo_file.replace('/', '-'), 'w')
    pickle.dump(freq_dict, wh)
    wh.close()
    print >> sys.stderr, 'count=', count
    
    try:
        os.remove('tmp_al.needle')
        os.remove('tmp_1hap.fas')
    except:
        pass
    
    return freq_dict



################
def freq_hap(sample_dir):
    """
    
    """
    
    import pickle
    from Bio import SeqIO
    from Bio import SeqUtils
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from MarkxIO import Markx10Iterator
    
    verbose = False
    thresh = 0.03
    freq_dict = {}
    
    tmp_list = os.listdir(sample_dir)
    wfiles = [ fn for fn in tmp_list if fn.startswith('w') and fn.endswith('haplo.out') ]
    totwf = len(wfiles)
    i = 0
    
    for hap in wfiles:
        i += 1
        print >> sys.stderr, 'Haplo file number', i, 'out of', totwf
        freq_dict[hap.rstrip('haplo.out')] = {}
        h = open(sample_dir + '/' + hap)
        hap_seq = list(SeqIO.parse(h, 'fasta'))
        
        for hap_item in hap_seq:
            oh = open('tmp_1hap.fas', 'w')
            oh.write('>%s\n' % hap_item.id)
            oh.write(hap_item.seq.tostring().replace('-', ''))
            oh.write('\n')
            oh.close()
            needle_align('tmp_1hap.fas', 'all_clones.fas', 'tmp_al.needle')
            hal = open('tmp_al.needle')
            min_dist = 1.0
            this_freq = float(hap_item.id.split('|')[1])
            al_records = list(Markx10Iterator(hal))
            assert len(al_records) == 10, 'not all records'
            
            for al in al_records:
                assert len(al.get_all_seqs()) == 2, '2 sequences per alignment'
                tmp = al.get_seq_by_num(0).tostring().upper()
                clone = al.get_seq_by_num(1).tostring().upper()
                query, match = list(al)
                # WARNING: not reliable on short sequences (a base might be absent, giving min=0)
                align_start = min([tmp.find('A'),tmp.find('C'),tmp.find('T'),tmp.find('G')])
                align_end = max([tmp.rfind('A'),tmp.rfind('C'),tmp.rfind('T'),tmp.rfind('G')])
                mm = [ 1 for c in zip(clone[align_start:align_end], tmp[align_start:align_end]) if c[0] != c[1] ]
                this_dist = float(len(mm)) / (align_end - align_start + 1)
                assert sum(mm) == len(mm), 'sum equals length'
                assert 0.0 <= this_dist and this_dist <= 1.0, 'not a distance'
                if this_dist < min_dist:
                    min_dist = this_dist
                    good_id = match.id.split('(')[0]

            if min_dist <= thresh:
                try:
                    freq_dict[hap.rstrip('haplo.out')][good_id] += this_freq
                except KeyError:
                    freq_dict[hap.rstrip('haplo.out')][good_id] = this_freq

    
    wh = open('%s.dict' % sample_dir.replace('/', '-'), 'w')
    pickle.dump(freq_dict, wh)
    wh.close()
    
    try:
        os.remove('tmp_al.needle')
        os.remove('tmp_1hap.fas')
    except:
        pass
    
    return freq_dict


def plot_freq_hap(freq_dict):
    """
    
    """
    try:
        import matplotlib
        matplotlib.use('pdf')
        import matplotlib.pyplot as plt
        import numpy as np
    except:
        print 'exit, could not import matplotlib'
        sys.exit()

    import math
    legend = False
    cb = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow']
    ca = ['#00BFFF', # alternative colour set
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
    xt_lab = []
    xt_val = []
    plt.subplot(111, axisbg='white')
    lmax = 0
    plt.title('clone mix')
    for k in freq_dict.keys():
        b = 0.0
        m = re.search('w(\d*)-(\d*)', k)
        l = 0.2*(int(m.group(2)) + int(m.group(1)))

        if l > lmax:
            lmax = l
        wd = 0.1*(int(m.group(2)) - int(m.group(1)))
        xt_val.append(l + 0.5*wd)
        xt_lab.append(m.group(1))
        #        minf = min(map(freq_dict[k])
        ci = 0
        
        for c_key in ref_genomes:
            try:
                height = freq_dict[k][c_key]
            except KeyError:
                height = 0.0
                """            try:
                h = math.log(height, 0.5)
                except ValueError:
                h = 0.0
                """
            h = height
            plt.bar(l, h, width=wd, bottom=b, color=ca[ci%cn], alpha=0.3)
            b += h
            ci += 1


    l = lmax + 5*wd

    b = 0
    ci = 0
    for height in geo_freq:
        h = height #math.log(height/minf, 2)
        plt.bar(l, h, width=wd, bottom=b, color=ca[ci%cn], alpha=0.3, label=ref_genomes[ci])
        b += h
        ci += 1
    plt.annotate('geometric\ndistribution', xy=(lmax+3*wd, 1.02), size=12)
    plt.xticks(xt_val, xt_lab ,rotation=90, size='x-small')
    plt.xlabel('window start')
    plt.ylabel('frequency')
    plt.axis([0, lmax + 10*wd, 0, 1])
    
    if legend:
        plt.subplot(212, frameon=False)
        ci = 0
        for height in geo_freq:
            h = 0
            tlabel = '%s    %4.2f %%' % (ref_genomes[ci], geo_freq[ci]*100)
            plt.bar(l, h, width=0, bottom=b, color=ca[ci%cn], alpha=0.3, label= tlabel)
            b += h
            ci += 1
        plt.axis([0, 6, 0, 1])
        plt.xticks('')
        plt.yticks('')
        params = {'legend.fontsize': 10}
        plt.rcParams.update(params)
    
        plt.legend()

    
    imtype = 'pdf'
    plt.savefig('freq_local_hap.%s' % imtype, dpi=None, facecolor='w', edgecolor='w',\
                    orientation='landscape', papertype=None, format=imtype,\
                    transparent=False)

    plt.show()
    
    
    return

###########################


def main():
    import pickle
    args = sys.argv
    try:
        haplo_file = args[1]
    except:
        sys.exit('usage: geo_chart.py haplo_file')
    
    if not os.path.exists(haplo_file):
        sys.exit('sample_directory does not exist')
    
    try:
        dict_file = haplo_file.replace('/', '-') + '.dict'
        dh = open(dict_file)
        freq_dict = pickle.load(dh)
        
    except:
        print >>sys.stderr, 'Running the alignments'
        freq_dict = freq_window(haplo_file)
        
        #   print freq_dict
    plot_pie_hap(freq_dict, haplo_file)

if __name__ == '__main__':
    main()
