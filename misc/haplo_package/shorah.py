#!/usr/bin/env python

# Copyright 2007, 2008
# Niko Beerenwinkel,
# Nicholas Eriksson,
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

import subprocess
import os
import sys

def run_dec(filein, fileref, j, a, w, wsh):
    """
    1. local error correction
    output stem.cor.fas
    """
    my_prog = "./dec.py"# % cwd
    my_arg  =  " -f %s -r %s -j %i -a %f -w %i -s %i" % (filein, fileref, j, a, w, wsh)
    assert os.path.isfile(filein), 'File %s not found' % filein
    
    # runs the local error correction
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child %s was terminated by signal" % my_prog, -retcode
        else:
            print >>sys.stderr, "Child %s returned %i" % (my_prog, retcode)
    except OSError, ee:
        print >>sys.stderr, "Execution of %s failed:" % my_prog, ee
        sys.exit()
    
    return retcode


def run_f2r(filein):
    """
    2. translate to read format
    output: filein.read
    """
    my_prog = "./fas2read.pl"
    my_arg  = " -f %s" % filein
    assert os.path.isfile(filein), 'File %s not found' % filein

    # runs the fasta to read conversion
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child %s was terminated by signal" % my_prog, -retcode
        else:
            print >>sys.stderr, "Child %s returned %i" % (my_prog, retcode)
    except OSError, ee:
        print >>sys.stderr, "Execution of %s failed:" % my_prog, ee
    
    return retcode


def run_contain(filein):
    """
    3. eliminate redundant reads
    output: filein.rest
    """
    my_prog = "./contain"
    my_arg  = " -f %s" % filein
    assert os.path.isfile('%s.fas' % filein), 'File %s not found' % filein

    #eliminates redundant reads
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child %s was terminated by signal" % my_prog, -retcode
        else:
            print >>sys.stderr, "Child %s returned %i" % (my_prog, retcode)
    except OSError, ee:
        print >>sys.stderr, "Execution of %s failed:" % my_prog, ee
    
    return retcode


def run_mm(filein, max_hap=10):
    """
    4. run maximum matching, output up to 200 haplotypes
    output: sample.$nr.geno
    """
    #my_prog = "python2.5 -m cProfile mm.py"
    my_prog = "python mm.py"
    my_arg  = " %s %i" % ('%s.rest' % filein, max_hap)

    assert os.path.isfile('%s.rest' % filein), 'File %s not found' % filein

    #eliminates redundant reads
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child %s was terminated by signal" % my_prog, -retcode
        else:
            print >>sys.stderr, "Child %s returned %i" % (my_prog, retcode)
    except OSError, ee:
        print >>sys.stderr, "Execution of %s failed:" % my_prog, ee

    return retcode


def run_freqEst(filein):
    """
    5. run EM
    output: sample.$nr.popl
    """
    my_prog = "./freqEst"
    my_arg  = " -f %s" % filein
    assert os.path.isfile('%s.rest' % filein), 'File %s not found' % filein

    #eliminates redundant reads
    try:
        retcode = subprocess.call(my_prog + my_arg, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child %s was terminated by signal" % my_prog, -retcode
        else:
            print >>sys.stderr, "Child %s returned %i" % (my_prog, retcode)
    except OSError, ee:
        print >>sys.stderr, "Execution of %s failed:" % my_prog, ee

    return retcode


def main():
    """ Only called if run not interactively
    """
    import optparse
    import shutil
    
    # parse command line
    optparser = optparse.OptionParser()
    
    optparser.add_option("-f", "--readsfile", help="file with reads <.sff or .fas format>",
                         default="", type="string", dest="f")
    optparser.add_option("-r", "--refgenome", help="file with reference genome <ref_genome.fas>",
                         default="ref_genome.fas", type="string", dest="r")
    optparser.add_option("-j", "--iterations", help="iterations in dpm sampling <1000>", default=1000,
                         type="int", dest="j")
    optparser.add_option("-a", "--alpha", help="alpha in dpm sampling <0.01>", default=0.01,
                         type="float", dest="a")
    optparser.add_option("-w", "--windowsize", help="window size in <99>", default=99,
                         type="int", dest="w")
    optparser.add_option("-s", "--winshifts", help="number of window shifts <3>", default=3,
                         type="int", dest="s") # window shiftings, such that each base is covered up to win_shifts times
    
    (options, args) = optparser.parse_args()
    
    try:
        [in_stem, in_format]  = [options.f.split('.')[0], options.f.split('.')[1]]
    except IndexError:
        print 'The input file must be filestem.format'
        sys.exit()

    retcode = run_dec(options.f, options.r, options.j, options.a, options.w, options.s)
    if retcode is not 0:
        sys.exit()
    shutil.copy('%s.cor.fas' % in_stem, '%s_cor.fas' % in_stem)
    
    retcode = run_f2r('%s_cor.fas' % in_stem)
    if retcode is not 0:
        sys.exit()
    
    retcode = run_contain('%s_cor' % in_stem)
    if retcode is not 0:
        sys.exit()
    
    retcode = run_mm('%s_cor' % in_stem)
    if retcode is not 0:
        sys.exit()
    
    retcode = run_freqEst('%s_cor' % in_stem)
    if retcode is not 0:
        sys.exit()

if __name__ == "__main__":
    main()
