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
import os

import logging
import logging.handlers

# Make a global logging object.
declog = logging.getLogger("DecLog")

#################################################
# a common user should not edit above this line #
#################################################
# parameters not controlled by command line options
fasta_length    = 60    # controls line length in fasta files
tolerance       = 0.1   # the portion of in-dels that are tolerated in a single read
go_default      = 6.0   # gap_opening penalty in needle alignment
ge_default      = 3.0   # gap_extension penalty in needle alignment
amb_thresh      = 2     # reads with more than amb_thresh Ns (ambiguous calls) are discarded
win_min_ext     = 0.85  # if the read covers at least win_min_ext of the window, fill it with Ns
verbose         = False # sets verbosity
hist_fraction   = 0.20  # fraction that goes into the history
min_quality     = 0.8   # quality under which discard the correction

# set logging level
declog.setLevel(logging.DEBUG)
# path to diri_sampler program
path_to_shorah  = './'
#################################################
# a common user should not edit below this line #
#################################################

# This handler writes everything to a file.
LOG_FILENAME = './dec.log'
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
declog.addHandler(h)

to_correct = {}
correction = {}
quality = {}

count = {}
count['A'] = 0
count['C'] = 0
count['G'] = 0
count['T'] = 0
count['X'] = 0
count['-'] = 0
clusters  = [ [] ]
untouched = [ [] ]

win_shifts = 3
keep_all_files  = False

def solve_amb(seq_list):
    """Use FASTA format description from NCBI
    to extract a random base from the possible ones
    A --> adenosine           M --> A C (amino)
    C --> cytidine            S --> G C (strong)
    G --> guanine             W --> A T (weak)
    T --> thymidine           B --> G T C
    U --> uridine             D --> G A T
    R --> G A (purine)        H --> A C T
    Y --> T C (pyrimidine)    V --> G C A
    K --> G T (keto)          N --> A G C T (any)
    -  gap of indeterminate length
    """
    
    import random
    solved_list = []
    f_code = {}
    f_code['R'] = ['G', 'A']
    f_code['Y'] = ['T', 'C']
    f_code['K'] = ['G', 'T']
    f_code['M'] = ['A', 'C']
    f_code['S'] = ['G', 'C']
    f_code['W'] = ['A', 'T']
    f_code['B'] = ['G', 'T', 'C']
    f_code['D'] = ['G', 'A', 'T']
    f_code['H'] = ['A', 'C', 'T']
    f_code['V'] = ['G', 'C', 'A']
    f_code['N'] = ['A', 'G', 'C', 'T']

    for base in seq_list:
        if base in f_code.keys():
            solved_list.append( random.choice(f_code[base]) )
        else:
            solved_list.append(base)
    return solved_list



def parse_aligned_reads(reads_file):
    """
    Parse reads from a file with aligned reads in fasta format
    """
    from Bio import SeqIO
    
    format = 'fasta'
    max_read_length = 300
    out_reads = {}
    cp = False
    
    if not os.path.isfile(reads_file):
        declog.error('There should be a file here: ' + reads_file)
        sys.exit('There should be a file here: ' + reads_file)
    else:
        declog.info( 'Using file ' + reads_file + ' of aligned reads')
    
    handle = open(reads_file, 'rU')
    declog.debug('Parsing aligned reads')
    
    for this_read in SeqIO.parse(handle, format):
        try:
            old_length = gen_length
        except NameError:
            gen_length = len(this_read.seq)
            old_length = gen_length
            declog.info('Alignment is ' + str(gen_length) + ' bases long')
            
        gen_length = len(this_read.seq)
        assert gen_length == old_length, 'All reads must have the same length'
        
        mst = this_read.seq.tostring()
        mls = list( mst )
        
        for c in mls:
            if c != '-' and c != 'X':
                mstart = mls.index(c) + 1
                break
        mstop = len ( mst.rstrip('-X') )
        
        # counts the gaps in the read (no flanking gaps)
        gaps_match = mst.strip('-X').count('-')
        '''
        for c in mst.strip('-X'):
            if c =='-':
                gaps_match = gaps_match + 1
        '''
        # match_length = len(mst.strip('-X'))
        # if gaps_match > round (tolerance * match_length):
        # print 'too many indels,', (gaps_query + gaps_match)
        #   continue
        
        name = this_read.id
        out_reads[name] = [ None, None, None, None, [] ]
        out_reads[name][0] = 0 # 
        out_reads[name][1] = gen_length # 
        out_reads[name][2] = mstart # this is
        out_reads[name][3] = mstop  # this too
        
        this_m = list( mst.strip('-X') )
        amb_calls = 0
        
        for i in range(len(this_m)):
            out_reads[name][4].append(this_m[i])
            
            if this_m[i] == 'N':
                amb_calls = amb_calls + 1
                # declog.debug('Found an N in ' + name)
        
        if amb_calls > amb_thresh:
            declog.warning('Read ' + name + ' has too many Ns')
            del out_reads[name]
    '''
    cp = open( './%s.covprof' % reads_file.rstrip('_reads.fas'), 'w' )
    for i in range(1, gen_length):
        cp.write('%i\t%i\n' % (i, cov_prof[i]) )
    cp.close()
    '''
    return out_reads


def print_window(wstart, wend, reads):
    """ 
    Considers only the read overlapping a given window and prints them aligned
    """
    wind_file = open('./w%d-%d.reads.fas' % (wstart, wend), 'w')
    rn = 0
    min_overlap = round( (wend - wstart + 1) * win_min_ext )

    for r in reads:
        rstart = reads[r][2]
        k = len(reads[r][4]) # trailing Ns have already been stripped
        
        if rstart <= wstart and rstart + k - 1 >= wend:
            overlap = wend - wstart + 1
        elif rstart >= wstart and rstart + k - 1 <= wend:
            overlap = k
        elif rstart <= wstart and rstart + k - 1 <= wend:
            overlap = rstart + k - wstart
        elif rstart >= wstart and rstart + k - 1 >= wend:
            overlap = wend - rstart + 1
        
        if overlap > min_overlap:
            rn = rn + 1
            wind_file.write( '>%s %d\n' % (r, reads[r][0]) )
            for i in range(wstart, wend+1):
                if i < rstart or i > rstart + k-1:
                    wind_file.write('N')
                else:
                    wind_file.write( '%s' % reads[r][4][i-rstart] )
                wind_file.flush()
            wind_file.write('\n')
    
    return rn


def run_dpm(run_setting):
    """run the dirichlet process clustering
    """

    import subprocess
    filein, j, a = run_setting
    dn = os.path.dirname(__file__)
    my_prog = os.path.join(dn, 'diri_sampler')
    my_arg =  ' -i %s -j %i -t %i -a %f' % (filein, j, int(j*hist_fraction), a)
    
    try:
        #os.remove('./corrected.tmp' )
        os.remove('./assignment.tmp')
    except:
        pass
    
    # runs the gibbs sampler for the dirichlet process mixture
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            declog.error('%s %s' % (my_prog, my_arg))
            declog.error('Child diri_sampler was terminated by signal %d' % -retcode)
        else:
            declog.debug('run %s finished' % my_arg)
            declog.debug('Child diri_sampler returned %i' % retcode)
    except OSError, ee:
        declog.error('Execution of diri_sampler failed: %s' % ee)
    
    return


def seq_2_fas(in_stem): # To be finished
    """ From Solexa _seq.txt files to fasta format
    """

    fas_handle = open('%s_reads.fas' % in_stem, 'w')
    
    declog.info('Converting files written in %s.seq' % in_stem + ' to fasta format')
    
    # print len( open('%s.seq' % in_stem, 'rU').readlines() ), ' _seq.txt files are to be read'
    
    file_names = open('%s.seq' % in_stem, 'rU')
    for line in file_names:
        this_seq = line.strip()
        
        file_seq = open(this_seq, 'rU')
        index = 0
        for line_seq in file_seq:
            line_entries = line_seq.split()
            if not line_entries:
                continue
            fas_handle.write( '>%s-read_%i\n' % (this_seq, index) )
            index = index + 1
            cc = 0
            read = line_entries[4]
            
            for c in read:
                fas_handle.write(c)
                cc = cc + 1
                if cc % fasta_length == 0:
                    fas_handle.write('\n')
            if cc % fasta_length != 0:
                    fas_handle.write('\n')
    return


def correct_reads(wstart, wend, k):
    """ Parses corrected.tmp (in fasta format) and correct the reads
    """
    #   out_reads[match_rec.id][0] = qstart
    #   out_reads[match_rec.id][1] = qstop
    #   out_reads[match_rec.id][2] = mstart
    #   out_reads[match_rec.id][3] = mstop
    #   out_reads[match_rec.id][4] = Sequence...
    
    from Bio import SeqIO
    import os
    import re
    read_name_rule = re.compile("(.*)\#\d+")
    wlen = wend - wstart + 1
    cor_file = 'w%d-%d.reads.cor.fas' % (wstart, wend)
    try:
        handle = open(cor_file, 'rU')
        for seq_record in SeqIO.parse(handle, 'fasta'):
            assert '\0' not in seq_record.seq.tostring(), 'binary file!!!'
            read_id = seq_record.id # read_name_rule.search(seq_record.id).group(1)
            try:
                correction[read_id][wstart] = list( seq_record.seq.tostring() )
                quality[read_id][wstart] = float(seq_record.description.split('|')[1].split('=')[1])
            except:
                correction[read_id] = {}
                quality[read_id] = {}
                correction[read_id][wstart] = list( seq_record.seq.tostring() )
                quality[read_id][wstart] = float(seq_record.description.split('|')[1].split('=')[1])
        handle.close()
        return
    except IOError:
        declog.warning('No reads in window %d?' % wstart)
        return
    

def get_prop(filename):
    """
    fetch the number of proposed clusters from .dbg file
    """
    h = open(filename)
    for l in h:
        if l.startswith('#made'):
            prop = int(l.split()[1])
            break
    h.close()
    try:
        return prop
    except UnboundLocalError:
        return 'not found'
        

def print_clust_unt(sh, clusters, untouched, in_stem, iterations):
    """
    Print clusters and untouched information
    in two separate files
    """
    del clusters[ len(clusters)-1 ]
    del untouched[ len(untouched)-1 ]
    ch = open('%s-%d.clusters' % (in_stem, sh), 'w')

    for i in range(1, 1 + iterations):
        ch.write('%d' % i)
        try:
            for k in clusters:
                ch.write('\t%s' % k[i])
        except:
            ch.write('\t0')
        
        ch.write('\n')
    uh = open('%s.%d.untouch' % (in_stem, sh), 'w')

    for i in range(1, 1 + iterations):
        uh.write('%d' % i)
        try:
            for k in untouched:
                uh.write('\t%s' % k[i])
        except:
            uh.write('\t0')
        uh.write('\n')
    return

def base_break(baselist):
    """
    """
    import random

    for c1 in count:
        count[c1] = 0
    for c in baselist:
        if c.upper() != 'N':
            count[c.upper()] += 1
        
    max = 0
    out = []
    for b in count:
        if count[b] >= max:
            max = count[b]
    for b in count:
        if count[b] == max:
            out.append(b)

    rc = random.choice(out)
    
    return rc

def main(fas_reads, step, win_shifts, keep_all_files, iters, alpha):
    """
    Performs the error correction analysis, running diri_sampler
    and analyzing the result
    """
    from multiprocessing import Pool
    
    if step % win_shifts != 0:
        sys.exit('Window size must be divisible by win_shifts')
        
    if win_min_ext < 1/float(win_shifts):
        declog.warning('Warning: some bases might not be covered by any window')
        
    assert os.path.isfile(fas_reads), "File '%s' not found" % fas_reads
    
    aligned_reads = parse_aligned_reads(fas_reads)
    r = aligned_reads.keys()[0]
    gen_length = aligned_reads[r][1] - aligned_reads[r][0]
   
    if step > gen_length:
        sys.exit('The window size must be smaller than the genome length')
 
    declog.info('%s reads are being considered' %  len(aligned_reads))
    # declog.info('The others discarded because of too many in-dels or ambiguous calls (Ns)')
    # print '\n'
    
    for k in aligned_reads:
        to_correct[k] = [None, None, None, None, []]
        to_correct[k][0] = aligned_reads[k][0]
        to_correct[k][1] = aligned_reads[k][1]
        to_correct[k][2] = aligned_reads[k][2]
        to_correct[k][3] = aligned_reads[k][3]
        to_correct[k][4] = [] # aligned_reads[k][4][:]

    ############################################
    # Now the windows and the error correction #
    ############################################
    
    # get the new length of the aligned genome (with gaps)
    # maybe there is a more elegant solution....
    gen_length = aligned_reads[aligned_reads.keys()[0]][1]
    
    single_sh = step/win_shifts
    proposed= {}
    runlist = []
    for sh in range(win_shifts):
        shift = sh * single_sh
        ncw = []
        for s in range(1, gen_length+1, step):
            rn = print_window(s+shift, s+shift+step-1, aligned_reads)
            declog.info('%d reads in %d' % (rn, s+shift))
            fst = 'w%d-%d.reads.fas' % (s+shift, s+shift+step-1)
            runlist.append([fst, iters, alpha])
            # if no reads in the window, go on
            if rn <= 1:
                declog.info('No reads in window %d-%d' % (s+shift, s+shift+step-1))
                os.remove(fst)
                continue
            
            declog.info('window %d to %d\t%d reads here' % (s+shift, s+shift+step-1, rn))
            ncw.append(s+shift)
            
    pool = Pool()
    pool.map(run_dpm, runlist) #(fst, options.j, options.a)
    for sh in range(win_shifts):
        shift = sh * single_sh
        for s in range(1, gen_length+1, step):
            correct_reads(s+shift, s+shift+step-1, sh)
            stem = 'w%d-%d' % (s+shift, s+shift+step-1)
            dbg_file = stem + '.reads.dbg'
            proposed[s+shift] = get_prop(dbg_file)
            if not keep_all_files:
                os.remove(dbg_file)
                os.remove(stem + '.reads.sam')
                try:
                    os.remove(stem + '.reads.fas')
                except:
                    pass
                if proposed[s+shift] != 'not found':
                    os.remove(stem + '.reads.cor.fas')
                    os.remove(stem + '-freq.csv')
                    os.remove(stem + '-support.fas')
    
    #################################
    ##  Print the corrected reads  ##
    #################################
    ck_not_found = 0
    for r in to_correct:
        if r not in correction.keys():
            ck_not_found += 1;
            continue
        rlen = len(aligned_reads[r][4])
        rst = aligned_reads[r][2]
        for rpos in range(rlen):
            this = []
            for cst in correction[r]:
                tp = rpos + rst - int(cst)
                if tp >=0 and tp < len(correction[r][cst]) and quality[r][cst] > min_quality:
                    t = correction[r][cst][tp]
                    this.append(t)
            if len(this) > 0:
                cb = base_break(this)
            else:
                cb = 'x'
            to_correct[r][4].append(cb)
            del this

    in_stem = fas_reads.split('.')[0]
    fch = open('%s.cor.fas' % in_stem, 'w')
    print >> sys.stderr, 'not found', ck_not_found, 'reads'
    
    for r in to_correct:
        if to_correct[r][4].count('x') == 0:
            fch.write('>%s %d\n' % (r, to_correct[r][2]) )
            cc = 0
            for c in to_correct[r][4]:
                if c != 'X':
                    fch.write(str(c))
                    fch.flush()
                    cc = cc + 1
                    if cc % fasta_length == 0:
                        fch.write('\n')
                        
            if cc % fasta_length != 0:
                fch.write('\n')
    fch.close()
    
    # write proposed_per_step to file
    ph = open('proposed.dat', 'w')
    ph.write('#base\tproposed_per_step\n')
    for kp in sorted(proposed.iterkeys()):
        if proposed[kp] != 'not found':
            ph.write('%s\t%f\n' % (kp, float(proposed[kp])/iters))
    ph.close()

if __name__ == "__main__":
    
    import optparse
    # parse command line
    optparser = optparse.OptionParser()
    
    optparser.add_option("-f", "--readsfile", help="file with reads <.far format>",
                         default="", type="string", dest="f")
    optparser.add_option("-j", "--iterations", help="iterations in dpm sampling <2000>", default=2000,
                         type="int", dest="j")
    optparser.add_option("-a", "--alpha", help="alpha in dpm sampling <0.01>", default=0.01,
                         type="float", dest="a")
    optparser.add_option("-w", "--windowsize", help="window size <201>", default=201,
                         type="int", dest="w")
    optparser.add_option("-s", "--winshifts", help="number of window shifts <3>", default=3,
                         type="int", dest="s") # window shiftings, such that each base is covered up to win_shifts times
    optparser.add_option("-k","--keep_files",help="keep all intermediate files of diri_sampler <default=False>", default=False,
                         action="store_true", dest="k")
    optparser.add_option("-r","--ref", type="string", default="",
                         dest="ref")
    optparser.add_option("-t","--threshold", help="if similarity is less, throw reads away... <default=0.7>", type="float",
                         dest="threshold", default=0.7)
    optparser.add_option("-d","--delete_s2f_files",help="delete temporary align files of s2f.py <default=False>", action="store_true", default=False,
                         dest="d")
    optparser.add_option("-n","--no_pad_insert",help="do not insert padding gaps in .far file<default=insert>", action="store_false", default=True,
                         dest="pad")

    (options, args) = optparser.parse_args()
    
    supported_formats = {
        'fas': 'fasta_reads',
        'far': 'fasta_aligned_reads'
        # 'seq': 'solexa_genome_analyzer_seq'
        }
    declog.info(' '.join(sys.argv))
    # check the input file is in supported format
    try:
        tmp_filename = os.path.split(options.f)[1]
        [in_stem, in_format]  = [tmp_filename.split('.')[0], tmp_filename.split('.')[-1]]
        t = supported_formats[in_format]
    except IndexError:
        declog.error('The input file must be filestem.format')
        print 'The input file must be filestem.format'
        print 'Supported formats are'
        for sf in supported_formats.iteritems():
            print sf[0], ':', sf[1]
        sys.exit()
    except KeyError:
        declog.error('format unknown')
        print in_format, 'format unknown'
        print 'Supported formats are'
        for sf in supported_formats.iteritems():
            print sf[0], ':', sf[1]
        sys.exit()

    if in_format != 'far':
        import s2f
        ref_file = options.ref
        out_file = os.path.join(os.getcwd(), in_stem+'.far')
        thresh = options.threshold
        pad_insert = options.pad
        declog.debug('running s2f.py')
        s2f.main(ref_file, options.f, out_file, thresh, pad_insert, keep_all_files)
        fas_reads = out_file
    else:
        fas_reads = options.f

    keep_all_files = options.k
    step = options.w
    win_shifts = options.s
    iters = options.j
    alpha = options.a
    main(fas_reads, step, win_shifts, keep_all_files, iters, alpha)
