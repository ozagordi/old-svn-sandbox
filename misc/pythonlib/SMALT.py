#!/usr/bin/env python

__author__ = "Osvaldo Zagordi"
__version__ = "$Revision: 0.1 $"
__date__ = "$Date: 2011/04/13$"
__copyright__ = ""
__license__ = ""
# Copyright 2011 by Osvaldo Zagordi

"""Command line wrapper for the read aligner SMALT
http://www.sanger.ac.uk/resources/software/smalt/
by Hannes Ponstingl @ Sanger Institute

Citations:
ftp://ftp.sanger.ac.uk/pub4/resources/software/smalt/smalt-manual-0.5.2.pdf

Last checked against version: 0.5.2
"""
import sys
import logging
import logging.handlers

LOG_FILENAME = './smalt-aligner.log'
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

def smalt_align(options):
    '''From fastq file to sam alignment
    '''
    import os.path
    
    if os.path.exists('/Users/ozagordi/Work/tools/smalt_MacOSX_i386'):
        smalt_exe = '/Users/ozagordi/Work/tools/smalt_MacOSX_i386'
    else:
        smalt_exe = 'smalt'
    
    # build index options
    ref_file = options.ref
    index_name = options.ind
    wordlen = options.k
    skipstep = options.s
    index_options = '-k %d -s %d' % (wordlen, skipstep)
    
    # build map options
    map_options = ' -f %s' % options.f
    if options.a: map_options += ' -a'
    if options.c: map_options += ' -c %f' % options.c
    if options.d: map_options += ' -d %d' % options.d
    if options.i: map_options += ' -i %d'% options.i
    if options.j: map_options += ' -j %d'% options.j
    if options.m: map_options += ' -m %d'% options.m
    if options.n: map_options += ' -n %d'% options.n
    if options.o: map_options += ' -o %d'% options.o
    if options.p: map_options += ' -p'
    if options.q: map_options += ' -q %d'% options.q
    if options.r: map_options += ' -r %d'% options.r
    if options.y: map_options += ' -y %f'% options.y
    if options.x: map_options += ' -x'
    if options.w: map_options += ' -w'
    
    # Runs smalt index, if necessary
    if not os.path.exists(index_name+'.sma') or not os.path.exists(index_name+'.smi'):
        cml = '%s index %s %s %s' % (smalt_exe, index_options, index_name, ref_file)
        generic_application(cml)

    # converts sff to fastq
    try:
        fastq_file = options.fq
        fastq_stem = ''.join(fastq_file.split('.')[:-1])
    except:
        sff_file = options.sff
        fastq_stem = ''.join(os.path.split(sff_file)[1].split('.')[:-1])
        fastq_file = fastq_stem + '.fastq'
        cml = 'sff2fastq -o %s %s' % (fastq_file, sff_file)
        generic_application(cml)

    # map reads
    try:
        read_file = options.fq
    except:
        read_file = options.ff
    cml = '%s map %s %s %s' % (smalt_exe, map_options, options.ind, read_file)
    generic_application(cml)
    
if __name__ == '__main__':
    import optparse
    # parse command line
    optparser = optparse.OptionParser()
    #input files
    optparser.add_option("--ff", "--fasta_file", help="file with reads in fastA format",
                         default="", type="string", dest="ff")
    optparser.add_option("--fq", "--fastq_file", help="file with reads in fastQ format",
                         default="", type="string", dest="fq")
    optparser.add_option("--sff", "--sff_file", help="file with reads in sff format (calls sff2fastq to convert)",
                         default="", type="string", dest="sff")
    optparser.add_option("--ref", "--ref_file", help="reference genome in fasta format",
                         type="string", default="", dest="ref")
    optparser.add_option("--ind", "--index_name", help="basename for index files",
                         type="string", default="", dest="ind")
    
    optparser.add_option("-k", "--wordlen", help="Sets the length of the hashed words, 2 < wordlen <= 20 [default=13]",
                         default=13, type="int", dest="k")
    
    optparser.add_option("-s", "--skipstep", help="sampling step size [default=2]",
                        default=2, type="int", dest="s")
    '''the distance between successive words that are hashed along the genomic reference sequence.
    With the option -s 1 every word is hashed, with -s 2 every second word, with -s 3 very third etc.
    By default skipstep is set equal to wordlen.'''
    
    optparser.add_option("-a", "--explicit", help="output explicit alignments",
                          action="store_true", default=False, dest="a")
    '''When this flag is set, explicit alignments are output along with the mappings.'''
    
    optparser.add_option("-c", "--mincover", help='''only consider mappings where the k-mer word seeds cover
                                                      the query read to a minimum extent''',
                          type="float", default=None, dest="c")
    '''mincover Only consider mappings where the k-mer word seeds cover the query read to a minimum extent.
    If mincover is an integer or floating point value > 1.0, at least this many bases of the read must be
    covered by k-mer word seeds. If mincover is a floating point value <= 1.0, it specifies the fraction of
    the query read length that must be covered by k-mer word seeds.
    '''
    
    optparser.add_option("-d", "--scorediff", help="threshold of the align score relative to maximum",
                          type="int", default=0, dest="d")
    '''scorediff Set a threshold of the Smith-Waterman alignment score relative to the maximum score. This
    option affects the way alignments are reported in single-read mode. For each read all alignments resulting
    in Smith- Waterman scores within scorediff of the maximum score are reported. Mappings with scores lower
    than this value are skipped. If scordiff (an integer) is set to a value < 0, all alignments are reported
    with scores above the threshold set by the -m minscor option. If set to 0 (default) only mappings with
    the best score are reported. Reads with multiple best mappings are reported as unmapped. This is also
    how read pairs are reported irrespective of the value of scorediff.
    '''
    
    optparser.add_option("-f", "--format", help="output format: cigar sam samsoft [default]",
                          type="string", default="samsoft", dest="f")
    '''format Specifies the output format. format can be one of the following strings:
    cigar (default) Compact Idiosyncratic Gapped Alignment Report (see http://www.sanger.ac.uk/resources/software/ssaha2)
    sam SequenceAlignment/Mapformat(http://samtools.sourceforge.net) with hard clipped sequences.
    samsoft like sam but using soft clipping ssaha native output format of the SSAHA2 software package
    (http://www.sanger.ac.uk/resources/software/ssaha2)
    '''
    
    optparser.add_option("-i", "--insertmax", help="max insert size for paired-end [default=500]",
                          default=None, type="int", dest="i")
    '''insertmax Maximum insert size for paired-end reads. insertmax is a positive integer (default 500).'''
    
    optparser.add_option("-j", "--insertmin", help="min insert size default[0]",
                          default=None, type="int", dest="j")
    '''insertmin Minimum insert size for paired-end reads insertmax is a positive integer (default 0).'''
    
    optparser.add_option("-m", "--minscor", help="absoult align score threshold to report alignments [default minscor = wordlen + skipstep - 1]",
                          type="int", default=None, dest="m")
    '''minscor Sets an absolute threshold of the Smith-Waterman scores. Mappings with scores below that
    threshold will not be reported. minscor is a positive integer (default minscor = wordlen + skipstep - 1).
    '''
    
    optparser.add_option("-n", "--nthreads", help="threads forked, 8 maximum",
                          type="int", dest="n")
    '''nthreads Run SMALT using multiple threads. nthread is the number of threads forked including the main thread.
    A maximum of 8 threads can be forked.
    '''
    
    optparser.add_option("-o", "--outfile", help="writes mapping to file [default to stdout]",
                          type="string", dest="o")
    '''oufilnam Write mapping output (e.g. SAM lines) to a separate file named oufilnam. If this option is not specified,
    mappings are written to standard output together with other messages.
    '''
    
    optparser.add_option("-p", "--partial", help="report partial alignments",
                          action="store_true", default=False, dest="p")
    '''Report partial alignments if they are complementary on the query read (split or chimeric reads). A maximum of
    two partial alignments are output per read. The second alignment is labelled P (-f ssaha or -f cigar formats)
    or has the secondary alignment bit-flag (0x100) of the SAM FLAG field raised (-f sam or -f samsoft).
    '''
    
    optparser.add_option("-q", "--minbasq", help="base quality threshold [default=0]",
                          type="int", default=0, dest="q")
    '''minbasq Sets a base quality threshold 0 <= minbasq <= 10 (default minbasq = 0). k-mer words of the read with base
    pairs that have a base quality below this threshold are not looked up in the hash index.
    '''
    
    optparser.add_option("-r", "--seed", help="same best alignment score is picked at random",
                          type="int", default=None, dest="r")
    '''seed If the there are multiple mappings with the same best alignment score report one one picked at random. This is
    relevant only in paired-end mode or with the option -d 0 (the default). seed is a positive integer used to to seed the
    pseudo-random generator. With seed = 0 a seed is derrived from the current calendar time. Without this option (default)
    reads with multiple best mappings are reported as unmapped.
    '''
    
    optparser.add_option("-y", "--minid", help="min exact matching nucleotides (fraction if <=1.0 or count otherwise)",
                          type="float", dest="y")
    '''minid Filters output alignment by a threshold in the number of exactly matching nucleotides. minid is a positive
    integer or a floating point number <= 1.0 specifying a fraction of the read length.
    '''
    
    optparser.add_option("-x", "--exhaustive", help="more exhaustive search for alignments",
                          action="store_true", default=True, dest="x")
    '''more exhaustive at the cost of decreased speed. In paired-end mode each mate is mapped independently.
    By default the mate with fewer hits in the hash index is mapped first and the vicinity is searched for its mate.
    '''
    
    optparser.add_option("-w", "--weighted", help="output complexity weighted Smith-Waterman scores",
                          action="store_true", default=False, dest="w")
                         
    (options, args) = optparser.parse_args()
    
    al_build = smalt_align(options)
    
    