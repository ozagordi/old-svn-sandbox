#!/usr/bin/env python

__author__ = "Osvaldo Zagordi"
__version__ = "$Revision: 0.1 $"
__date__ = "$Date: 2011/02/18$"
__copyright__ = ""
__license__ = ""
# Copyright 2011 by Osvaldo Zagordi

"""Command line wrapper for the read aligner Mosaik
http://code.google.com/p/mosaik-aligner
by MarthLab http://bioinformatics.bc.edu/marthlab/

Citations:

Last checked against version: 1.1.0021
"""
import sys
import logging
import logging.handlers

LOG_FILENAME = './mosaik-aligner.log'
# Make a global logging object
x = logging.getLogger("logfun")
x.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler('./log', 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s\t%(asctime)s\n\t%(message)s\n")
h.setFormatter(f)
x.addHandler(h)
logfun = logging.getLogger("logfun")

def generic_application(cml):
    import subprocess
    logfun.info('Running: %s' % cml)
    
    p = subprocess.Popen(cml, shell='/bin/bash',# bufsize=bufsize,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    p.wait()
    if p.returncode < 0: logfun.critical('Child terminated')
    oo = unicode(p.stdout.read(), errors='replace')
    logfun.info(oo)
    logfun.info('*'*80)

def mosaik_align(options):
    '''From fastq file to sam alignment
    '''
    import os.path
    
    ref_file = options.ref    
    #Runs MosaikBuild on the reference, if necessary
    if ref_file.endswith('.dat'):
        ref_dat = ref_file
    else:
        ref_dat = os.path.split(ref_file)[1].split('.')[-2] + '.dat'
        cml = 'MosaikBuild -fr %s -oa %s' % (ref_file, ref_dat)
        generic_application(cml)

    try:
        fastq_file = options.q
        fastq_stem = fastq_file.split('.')[-2]
    except:
        sff_file = options.s
        fastq_stem = os.path.split(sff_file)[1].split('.')[-2]
        fastq_file = fastq_stem + '.fastq'
        cml = 'sff2fastq -o %s %s' % (fastq_file, sff_file)
        generic_application(cml)


    
    #Runs MosaikBuild on the fastq
    seq_build = fastq_stem  + '.mkb'
    cml = 'MosaikBuild -q %s -st %s -out %s' % (fastq_file, options.st, seq_build)
    generic_application(cml)
    
    #Runs MosaikAlign
    al_build = fastq_stem + '.mka'
    cml = 'MosaikAligner -in %s -ia %s -out %s' % (seq_build, ref_dat, al_build)
    cml += ' -p %d' % options.p
    cml += ' -bw %d' % options.bw
    if options.mmp: cml += ' -mmp %f' % options.mmp
    if options.mmal: cml += ' -mmal'
    generic_application(cml)
    
    #Runs Sort
    al_sorted = fastq_stem + '.mks'
    cml = 'MosaikSort -in %s -out %s' % (al_build, al_sorted)
    generic_application(cml)
    
    #Runs Text (convert to sam)
    cml ='MosaikText -in %s -sam %s.sam' % (al_sorted, fastq_stem)
    generic_application(cml)

if __name__ == '__main__':
    import optparse
    # parse command line
    optparser = optparse.OptionParser()
    #input files
    optparser.add_option("-f", "--fasta_file", help="file with reads in fastA format",
                         default="", type="string", dest="f")
    optparser.add_option("-q", "--fastq_file", help="file with reads in fastQ format",
                         default="", type="string", dest="q")
    optparser.add_option("-s", "--sff_file", help="file with reads in sff format (calls sff2fastq to convert",
                         default="", type="string", dest="s")    
    optparser.add_option("-r","--ref", help="reference genome in fasta or .dat format (Mosaik database)",
                         type="string", default="", dest="ref")
    optparser.add_option("--st","--platform", help="specifies sequencing technology\
                            <454(default), helicos, illumina, sanger, solid>",
                         type="string", default="454", dest="st")
    
    optparser.add_option("-m", "--mode", help="alignment mode <unique, all(default)>",
                         default="all", type="string", dest="m")
    optparser.add_option("--hs", "--hashsize", help="hash size [4-32] (default 15)",
                         default=15, type="int", dest="hs")
    optparser.add_option("-p", "--processors", help=" number of processors (default 1)",
                         default=1, type="int", dest="p")
    optparser.add_option("--bw", "--bandwidth", help="banded Smith-Waterman for better performance",
                         default=51, type="int", dest="bw")
    
    optparser.add_option("--mm","--mismatches", help="allowed mismatches (default 4)",
                         dest="k", default=4)
    optparser.add_option("--mmp","--mismatch_percentage", help="in fraction rather than count [0.0-1.0] (default 0.1)",
                         default=0.1, type="float", dest="mmp")
    optparser.add_option("--mmal","--mismatch_on_align", help="count the mismatches only on the aligned portion",
                         action="store_true", default=True, dest="mmal")
    optparser.add_option("--minp","--min_percent_aligned",help="specifies min fraction of the read that should be aligned [0.0-1.0]",
                         dest="minp")
                         
    (options, args) = optparser.parse_args()
    
    al_build = mosaik_align(options)
    
    