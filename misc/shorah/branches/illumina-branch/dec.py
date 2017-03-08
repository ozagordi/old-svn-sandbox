#!/usr/bin/env python

# Copyright 2007, 2008
# Niko Beerenwinkel,
# Nicholas Eriksson,
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



#################################################
# a common user should not edit above this line #
#################################################
# parameters not controlled by command line options
fasta_length    = 60    # controls line length in fasta files
tolerance       = 0.1   # the portion of in-dels that are tolerated in a single read
go_default      = 6.0   # gap_opening penalty in needle alignment
ge_default      = 3.0   # gap_extension penalty in needle alignment
amb_thresh      = 2     # reads with more than amb_thresh Ns (ambiguous calls) are discarded
win_min_ext     = 0.7   # if the read covers at least win_min_ext of the window, fill it with Ns
verbose         = False # sets verbosity
hist_fraction   = 0.20  # fraction that goes into the history
min_quality     = 0.8   # quality under which discard the correction

# path to diri_sampler program
path_to_shorah  = './'

#################################################
# a common user should not edit below this line #
#################################################

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
    import time
    
    format = 'fasta'
    max_read_length = 300
    out_reads = {}
    cp = False
    
    if not os.path.isfile(reads_file):
        print 'There should be a file here:', reads_file
        sys.exit()
    else:
        print 'Using file', reads_file, 'of aligned reads'
    
    handle = open(reads_file, 'rU')
    print 'Parsing aligned reads'
    
    for this_read in SeqIO.parse(handle, format) :

        try:
            old_length = gen_length
        except:
            gen_length = len(this_read.seq)
            old_length = gen_length
        gen_length = len(this_read.seq)
        
        assert gen_length == old_length, 'All reads must have the same length'
        
        if not cp: # only the first time allocate cov_prof
            cov_prof = [0]*(2*gen_length + max_read_length)
            cp = True
        #print 'Alignment is', gen_length, 'bases long'
        
        mst = this_read.seq.tostring()
        
        mls = list( mst )
        
        for c in mls:
            if c != '-' and c != 'X':
                mstart = mls.index(c) + 1
                break
        mstop = len ( mst.rstrip('-X') )
#        if mst.startswith('-') and mstart > 20:# or mst.startswith('-'):
 #           print mst[mstart-10: mstop+10]
  #          sys.exit()
        
        # counts the gaps in the read (no flanking gaps)
        gaps_match = 0
        for c in mst.strip('-X'):
            if c =='-':
                gaps_match = gaps_match + 1
        
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
        
        for i in range(mstart, mstop + 1):
            try:
                cov_prof[i] = cov_prof[i] + 1
            except IndexError:
                print 'out of coverage', i
        this_m = list( mst.strip('-X') )
                
        amb_calls = 0
        
        for i in range(len(this_m)):
            out_reads[name][4].append(this_m[i])
            
            if this_m[i] == 'N':
                amb_calls = amb_calls + 1
                if verbose:
                    print >> sys.stderr, 'Found an N in', name
        
        if amb_calls > amb_thresh:
            if verbose:
                print 'Read', name, 'has too many Ns'
            del out_reads[name]
    cp = open( './%s.covprof' % reads_file.rstrip('_reads.fas'), 'w' )
    for i in range(1, gen_length):
        cp.write('%i\t%i\n' % (i, cov_prof[i]) )
    cp.close()
    
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


def run_dpm(filein, j, a):
    """run the dirichlet process clustering
    """

    import subprocess
    
    my_prog = "%s/diri_sampler" % path_to_shorah
    my_arg =  " -i %s -j %i -t %i -a %f > ds-out" % (filein, j, int(j*hist_fraction), a)
    
    try:
        #os.remove('./corrected.tmp' )
        os.remove('./assignment.tmp')
    except:
        pass
    
    # runs the gibbs sampler for the dirichlet process mixture
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            print >>sys.stderr, "'%s %s'"%(my_prog,my_arg)
            print >>sys.stderr, "Child diri_sampler was terminated by signal", -retcode
        else:
            print >>sys.stdout, "Child diri_sampler returned %i" % retcode
    except OSError, ee:
        print >>sys.stderr, "Execution of diri_sampler failed:", ee
    
    assert os.path.isfile('ds-out'), 'File ds-out not found'
    new_clusters = 0
    k = len(clusters)-1
    dsh = open('ds-out', 'rU')
    for line in dsh:
        if line.startswith('iteration'):
            lspl = line.split()
            it = lspl[1]
            clus = lspl[2]
            untouch = lspl[3]
            clusters[k].append(clus)
            untouched[k].append(untouch)
        elif line.startswith('#made'):
            new_clusters = int(line.split()[1])
    assert len(clusters[k]) > 1, 'Clusters not appended'
    assert len(untouched[k]) > 1, 'Untouched not appended'
    clusters.append([])
    untouched.append([])
    return new_clusters



def seq_2_fas(in_stem): # To be finished
    """ From Solexa _seq.txt files to fasta format
    """

    fas_handle = open('%s_reads.fas' % in_stem, 'w')
    
    print 'Converting files written in', '%s.seq' % in_stem, 'to fasta format'
    
    print len( open('%s.seq' % in_stem, 'rU').readlines() ), ' _seq.txt files are to be read'
    
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
    try:
        handle = open('corrected.tmp', 'rU')
        for seq_record in SeqIO.parse(handle, "fasta"):
            assert '\0' not in seq_record.seq.tostring(), 'binary file!!!'
            read_id = seq_record.id # read_name_rule.search(seq_record.id).group(1)
            try:
                correction[read_id][wstart] = list( seq_record.seq.tostring() )
                quality[read_id][wstart] = float(seq_record.description.split('|')[1])
            except:
                correction[read_id] = {}
                quality[read_id] = {}
                correction[read_id][wstart] = list( seq_record.seq.tostring() )
                quality[read_id][wstart] = float(seq_record.description.split('|')[1])
        
        handle.close()
        #os.remove('corrected.tmp')
        return
    except IOError:
        print 'No reads in window %d?' % wstart
        return
    

def print_prop(sh, ncw, nc, in_stem):
    """
    print a single file with proposed new clusters
    """
    nch = open('%s-%d.proposed' % (in_stem, sh), 'w')
    nch.write('#Proposed new clusters in each read\n')
    for nck in ncw:
        nch.write('%i\t%i\n' % (nck, nc[nck]))
    nch.close()
    return


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

def main():
    """ Only called if run not interactively
    """
    import optparse
     
    supported_formats = {
        'far': 'fasta_aligned_reads'
        # 'seq': 'solexa_genome_analyzer_seq'
        }
    
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
    optparser.add_option("-k","--keep_all_files",help="keep all files of the diri_sampler output. useful for detailed analysis. <default=delete>",
                         action="store_true",dest="k")
    
    (options, args) = optparser.parse_args()
    
    # check the input file is in supported format
    try:
        tmp_filename = os.path.split(options.f)[1]
        [in_stem, in_format]  = [tmp_filename.split('.')[0], tmp_filename.split('.')[-1]]
        t = supported_formats[in_format]
    except IndexError:
        print 'The input file must be filestem.format'
        print 'Supported formats are'
        for sf in supported_formats:
            print sf+': ', supported_formats[sf]
        sys.exit()
    except KeyError:
        print in_format, 'format unknown'
        print 'Supported formats are'
        for sf in supported_formats:
            print sf+': ', supported_formats[sf]
        sys.exit()
    cwd = os.getcwd()

    if not os.path.isfile('%s/diri_sampler' % path_to_shorah):
        print "found no 'diri_sampler'-program in the directory '%s'" % path_to_shorah
        sys.exit()

    keep_all_files = False
    try:
        if options.k == True:
            keep_all_files = True
    except:
        pass

    step = options.w
    win_shifts = options.s
    if step % win_shifts != 0:
        sys.exit('Window size must be divisible by win_shifts')
    
    if win_min_ext < 1/float(win_shifts):
        print 'Warning: some bases might not be covered by any window'
    
    assert os.path.isfile(options.f), "File '%s' not found" % options.f
    
    fas_reads = options.f


    if in_format == 'far':
        aligned_reads = parse_aligned_reads(fas_reads)
        r = aligned_reads.keys()[0]
        gen_length = aligned_reads[r][1] - aligned_reads[r][0]
    else:
        sys.exit('Input file must be a .far file (see the README)')
   
    if options.w >= gen_length:
        sys.exit('The window size must be smaller than the genome length')
 
    print len(aligned_reads), 'reads are being considered'
    print 'The others discarded because of too many in-dels or ambiguous calls (Ns)'
    print '\n'
    
    for k in aligned_reads:
        to_correct[k] = [None, None, None, None, []]
        to_correct[k][0] = aligned_reads[k][0]
        to_correct[k][1] = aligned_reads[k][1]
        to_correct[k][2] = aligned_reads[k][2]
        to_correct[k][3] = aligned_reads[k][3]
        to_correct[k][4] = [] # aligned_reads[k][4][:]     
    del k

    ############################################
    # Now the windows and the error correction #
    ############################################

    # get the new length of the aligned genome (with gaps)
    # maybe there is a more elegant solution....
    gen_length = aligned_reads[aligned_reads.keys()[0]][1]


    single_sh = step/win_shifts
    proposed= {}
    for sh in range(win_shifts):
        shift = sh * single_sh
        ncw = []
        for s in range(1, gen_length+1, step):
            rn = print_window(s+shift, s+shift+step-1, aligned_reads)
            fst = 'w%d-%d.reads.fas' % (s+shift, s+shift+step-1)
        
        # print '\x1B[1A\x1B[2K doing from %d to %d    %d reads in this window\n' % (s, s+step, rn)
        
        # if no reads in the window, go on
            if rn <= 1:
                print >>sys.stdout, '\n No reads in window', s+shift, '-', s+shift+step-1
                os.remove(fst)
                continue
        
            print >>sys.stdout, '\n window from %d to %d    %d reads in this window' % (s+shift, s+shift+step-1, rn)

            ncw.append(s+shift)
            proposed[s+shift] = run_dpm(fst, options.j, options.a)
            correct_reads(s+shift, s+shift+step-1, sh)
            if keep_all_files:
                os.rename('ds-out','w%d-%d.ds-out'%(s+shift,s+shift+step-1))
                os.rename('corrected.tmp','w%d-%d.corrected.fas'%(s+shift,s+shift+step-1))
                os.rename('haplo.out','w%d-%d.haplo.out'%(s+shift,s+shift+step-1))
            else:
                os.remove('ds-out')
                os.remove('./%s' % fst)
                os.remove('corrected.tmp')
                os.remove('haplo.out')

        try:
            os.remove('./assignment.tmp')
        except:
            pass
        
        

    #print_prop(sh, ncw, nc, in_stem)
        #del nc
        #del ncw
        #print_clust_unt(sh, clusters, untouched, in_stem, options.j)

    #################################
    ##  Print the corrected reads  ##
    #################################
    
    for r in to_correct:
        if r not in correction.keys():
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
            
    fch = open('%s.cor.fas' % in_stem, 'w')
    
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
    
    
    #################
    # gnuplot files #
    #################
    
    # clusters

    prop_keys = proposed.keys()
    prop_keys.sort()
    pf = open ('%s.proposed' % in_stem, 'w')
    for k in prop_keys:
        pf.write('%d\t%f\n' % (k, float(proposed[k])/options.j))
    pf.close()

    set_prop = 'set title \'proposed per step, sample %s\'\nset yra [0:]\n' % (in_stem)
    title = ''
    gf = open('proposed.gp', 'w')
    
    gf.write(set_prop)
    gf.write( 'plot \'%s.proposed\' u 1:2 w l title \'%s\'\n' % (in_stem, title) )
    gf.close()


if __name__ == "__main__":
    main()
