#!/usr/bin/env python

import sys
import os
homedir = os.path.expanduser('~/')
sys.path.append(homedir)
import pythonlib

Verbose = False

clones_order = {
    '07-56681': 0, # 50%
    '07-56951': 1, # 25%
    '08-59712': 2, # 12,5%
    '07-54825': 3, # 6,25%
    '08-04134': 4, # 3,125%
    '08-01315': 5, # 1,562%
    '08-55163': 6, # 0,781%
    '08-57881': 7, # 0,39%
    '08-02659': 8, # 0,2%
    '08-04512': 9  # 0,1%
    }

short_clones_order = {
    '56681': 0, # 50%
    '56951': 1, # 25%
    '59712': 2, # 12,5%
    '54825': 3, # 6,25%
    '04134': 4, # 3,125%
    '01315': 5, # 1,562%
    '55163': 6, # 0,781%
    '57881': 7, # 0,39%
    '02659': 8, # 0,2%
    '04512': 9  # 0,1%
    }
order_clones = [None]*10
for c in clones_order:
    order_clones[clones_order[c]] = c


def score_dict(lane):
    import glob
    from cPickle import load
    from heapq import nlargest
    import operator


    score = {}
    files_a = glob.glob('./l_%d-*-a.pck' % lane)
    files_b = glob.glob('./l_%d-*-b.pck' % lane)
    assert len(files_a) == len(clones_order), 'some index has not been read [a]'+str(files_a)
    assert len(files_b) == len(clones_order), 'some index has not been read [b]'
    
    for pair in zip(files_a, files_b):
        a = pair[0]
        b = pair[1]
        ind_id = a.split('-')[2]
        da = load(open(a))
        db = load(open(b))
        
        order = short_clones_order[ind_id]
        
        for ka in da:
            read_id = ka[:-2]
            kb = read_id + '-b'
            if kb in db:
                try:
                    score[read_id]['a'][order] = da[ka]
                except:
                    score[read_id] = {}
                    score[read_id]['a'] = [-1]*10
                    score[read_id]['b'] = [-1]*10
                    score[read_id]['a'][order] = da[ka]
                score[read_id]['b'][order] = db[kb]
        
    handle = open('best-l_%d' % lane, 'w')
    handle.write('sample-lane-x-y-id\ta_clone1\tscore\ta_clone2\tscore\tb_clone1\tscore\tb_clone2\tscore\n')
    for read in score.keys():
        sa = enumerate(score[read]['a'])
        sb = enumerate(score[read]['b'])
        best2_a = nlargest(2, sa, operator.itemgetter(1))
        best2_b = nlargest(2, sb, operator.itemgetter(1))
        ca_1, ca_2 = order_clones[best2_a[0][0]], order_clones[best2_a[1][0]]
        cb_1, cb_2 = order_clones[best2_b[0][0]], order_clones[best2_b[1][0]]
        handle.write('%s\t%s\t%d\t%s\t%d\t' % (read, ca_1, best2_a[0][1], ca_2, best2_a[1][1]))
        handle.write('%s\t%d\t%s\t%d\n' % (cb_1, best2_b[0][1], cb_2, best2_b[1][1]))
    handle.close()


def parse_com_line():
    from optparse import OptionParser
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage = usage)
    
    parser.add_option("-s", "--sample_directory", dest='s', default='../../../data/Sample-174/',
                      metavar="SAMPLE_DIR", help="Illumina files s_x_n_seq.txt are in SAMPLE_DIR [default: %default]"),

    parser.add_option("-c", "--clones", dest='c', default='../../../data/martin_references/martin_clones_nidx/',
                      metavar="CL_DIR", help="clones (haplotypes) file [default: %default]"),
    
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="verbose behaviour [default: %default]")
    
    (options, args) = parser.parse_args()
    
    return options, args

def split2fasta(sample_dir, lane):
    '''
    '''
    import tempfile
    fa = tempfile.NamedTemporaryFile(delete=False)
    fb = tempfile.NamedTemporaryFile(delete=False)
    # fa = open('tmp-%d-a.fas' % lane, 'w')
    # fb = open('tmp-%d-b.fas' % lane, 'w')
    s_files = [ f for f in os.listdir(sample_dir) if f.startswith('s_') and f.endswith('txt') ]
    sample = int(s_files[0].split('_')[1])
    
    r_file = os.path.join(sample_dir, 's_%d_%4.4d_seq.txt' % (sample, lane))
    ih = open(r_file)
    for l in ih:
        lsp = l.strip().split()
        id = 's_%s-l_%s-x_%s-y_%s' % (lsp[0], lsp[1], lsp[2], lsp[3])
        r_len = len(lsp[4])
        fa.write('>%s-a\n' %id)
        fb.write('>%s-b\n' %id)
        fa.write('%s\n' % lsp[4][:r_len/2])
        fb.write('%s\n' % lsp[4][r_len/2:])
    fa.close()
    fb.close()
    
    return fa.name, fb.name

def save_al_score(al_file, d_al_file):
    '''
    Reads alignment file in Novoalign native format and
    returns dictionary of the scores
    '''
    from pythonlib.Novoalign import NovoalignNativeIterator as NNI
    from cPickle import dump
    
    n = NNI(open(al_file))
    score = {}
    
    while True:
        try:
            a = n.next()
        except:
            break
        if a == None:
            break

        if a._status == 'unique':
            score[a._id] = a._score
        
    h = open(d_al_file, 'w')
    dump(score, h)
    h.close()
    
    return
    

def align_lane(al_args):
    '''
    Splits, align with novoalign
    '''
    from pythonlib import Novoalign
    import subprocess
    
    lane, sample_dir, index_file = al_args
    index_suffix = index_file.split('/')[-1]
    al_file_a = 'l_%d-i_%s-a.noa' % (lane, index_suffix)
    d_al_file_a = 'l_%d-i_%s-a.pck' % (lane, index_suffix)
    al_file_b = 'l_%d-i_%s-b.noa' % (lane, index_suffix)
    d_al_file_b = 'l_%d-i_%s-b.pck' % (lane, index_suffix)
    
    if not os.path.exists(d_al_file_a) or not os.path.exists(d_al_file_b):
        fa_name, fb_name = split2fasta(sample_dir, lane)

    #If file doesn' exist, a
    if not os.path.exists(d_al_file_a):
        cml = str(Novoalign.NovoalignCommandline(database=index_file, readfile=fa_name))
        cml += ' > %s' % al_file_a
        print >> sys.stderr, cml
        subprocess.check_call(cml, shell=True)
        save_al_score(al_file_a, d_al_file_a)
    try:
        os.unlink(al_file_a)
    except:
        pass
        
    #If file doesn' exist, b
    if not os.path.exists(d_al_file_b):
        cml = str(Novoalign.NovoalignCommandline(database=index_file, readfile=fb_name))
        cml += ' > %s' % al_file_b
        print >> sys.stderr, cml
        subprocess.check_call(cml, shell=True)
        save_al_score(al_file_b, d_al_file_b)
    try:
        os.unlink(al_file_b)
    except:
        pass
    
    try:
        os.unlink(fa_name)
        os.unlink(fb_name)
    except:
        pass

    return

def main():
    '''
    '''
    
    from multiprocessing import Pool
    
    options, args = parse_com_line()
    sample_dir = options.s
    index_dir = options.c
    
    first_lane, last_lane = 1, 100
    
    for ind_id in clones_order.keys():
        index = os.path.join(index_dir, ind_id)
        lanes = range(first_lane, last_lane+1)
        al_args = [(i, sample_dir, index) for i in lanes]
        pool = Pool()
        it = pool.map(align_lane, al_args)
    
    print >> sys.stderr, 'alignments read'
    pool = Pool()
    score_args = [i for i in range(first_lane, last_lane+1)]
    it = pool.map(score_dict, score_args)
#    score = score_dict(range(first_lane, last_lane+1))

if __name__ == '__main__':
    main()
