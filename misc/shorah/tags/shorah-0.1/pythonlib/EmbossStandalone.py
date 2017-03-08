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


import MyAlignIO
import subprocess
import os
#import sys


def needle(needlecmd, aseqfile, bseqfile, **keywds):
    """  Standard (Mandatory) qualifiers:
  [-asequence]         sequence   Sequence filename and optional format, or
                                  reference (input USA)
  [-bsequence]         seqall     Sequence(s) filename and optional format, or
                                  reference (input USA)
   -gapopen            float      [10.0 for any sequence] The gap open penalty
                                  is the score taken away when a gap is
                                  created. The best value depends on the
                                  choice of comparison matrix. The default
                                  value assumes you are using the EBLOSUM62
                                  matrix for protein sequences, and the
                                  EDNAFULL matrix for nucleotide sequences.
                                  (Number from 0.000 to 100.000)
   -gapextend          float      [0.5 for any sequence] The gap extension
                                  penalty is added to the standard gap penalty
                                  for each base or residue in the gap. This
                                  is how long gaps are penalized. Usually you
                                  will expect a few long gaps rather than many
                                  short gaps, so the gap extension penalty
                                  should be lower than the gap penalty. An
                                  exception is where one or both sequences are
                                  single reads with possible sequencing
                                  errors in which case you would expect many
                                  single base gaps. You can get this result by
                                  setting the gap open penalty to zero (or
                                  very low) and using the gap extension
                                  penalty to control gap scoring. (Number from
                                  0.000 to 10.000)
  [-outfile]           align      [*.water] Output alignment file name

   Additional (Optional) qualifiers:
   -datafile           matrixf    [EBLOSUM62 for protein, EDNAFULL for DNA]
                                  This is the scoring matrix file used when
                                  comparing sequences. By default it is the
                                  file 'EBLOSUM62' (for proteins) or the file
                                  'EDNAFULL' (for nucleic sequences). These
                                  files are found in the 'data' directory of
                                  the EMBOSS installation.

   Advanced (Unprompted) qualifiers:
   -[no]brief          boolean    [Y] Brief identity and similarity

   Associated qualifiers:

   "-asequence" associated qualifiers
   -sbegin1            integer    Start of the sequence to be used
   -send1              integer    End of the sequence to be used
   -sreverse1          boolean    Reverse (if DNA)
   -sask1              boolean    Ask for begin/end/reverse
   -snucleotide1       boolean    Sequence is nucleotide
   -sprotein1          boolean    Sequence is protein
   -slower1            boolean    Make lower case
   -supper1            boolean    Make upper case
   -sformat1           string     Input sequence format
   -sdbname1           string     Database name
   -sid1               string     Entryname
   -ufo1               string     UFO features
   -fformat1           string     Features format
   -fopenfile1         string     Features file name

   "-bsequence" associated qualifiers
   -sbegin2            integer    Start of each sequence to be used
   -send2              integer    End of each sequence to be used
   -sreverse2          boolean    Reverse (if DNA)
   -sask2              boolean    Ask for begin/end/reverse
   -snucleotide2       boolean    Sequence is nucleotide
   -sprotein2          boolean    Sequence is protein
   -slower2            boolean    Make lower case
   -supper2            boolean    Make upper case
   -sformat2           string     Input sequence format
   -sdbname2           string     Database name
   -sid2               string     Entryname
   -ufo2               string     UFO features
   -fformat2           string     Features format
   -fopenfile2         string     Features file name

   "-outfile" associated qualifiers
   -aformat3           string     Alignment format
   -aextension3        string     File name extension
   -adirectory3        string     Output directory
   -aname3             string     Base file name
   -awidth3            integer    Alignment width
   -aaccshow3          boolean    Show accession number in the header
   -adesshow3          boolean    Show description in the header
   -ausashow3          boolean    Show the full USA in the alignment
   -aglobal3           boolean    Show the full sequence in alignment

   General qualifiers:
   -auto               boolean    Turn off prompts
   -stdout             boolean    Write standard output
   -filter             boolean    Read standard input, write standard output
   -options            boolean    Prompt for standard and additional values
   -debug              boolean    Write debug output to program.dbg
   -verbose            boolean    Report some/full command line options
   -help               boolean    Report command line options. More
                                  information on associated and general
                                  qualifiers can be found with -help -verbose
   -warning            boolean    Report warnings
   -error              boolean    Report errors
   -fatal              boolean    Report fatal errors
   -die                boolean    Report dying program messages

    """
    att2param = {
        'seq1'     : '-asequence',
        'seq2'     : '-bsequence',
        
        'gapopen'  : '-gapopen',
        'gapext'   : '-gapextend',
        'datafile' : '-datafile',
        
        'sb1'      : '-sbegin1',
        'se1'      : '-send1',

        'sb2'      : '-sbegin2',
        'se2'      : '-send2',
        
        'out'      : '-outfile',
        'stdout'   : '-stdout', # TODO, doesn't work yet
        'outform'  : '-aformat3',
        'aglobal3' : '-aglobal3',

        'usashow'  : '-usashow'
        # and more to come
        }

    if not os.path.exists(needlecmd):
        raise ValueError, "needle does not exist at %s" % needlecmd
    
    params = []

    params.extend([att2param['seq1'], aseqfile])
    params.extend([att2param['seq2'], bseqfile])
    params.extend([att2param['gapopen'], str(10.0)])
    params.extend([att2param['gapext'], str(0.5)])
#    params.extend([att2param['out'], aseqfile.rstrip('.fas') + '.water'])
    params.extend([att2param['outform'], 'markx10'])

    for attr in keywds.keys():
        params.extend([att2param[attr], str(keywds[attr])])

    bufsize = -1
    ind  = None
    outd = None
    errd = open('/dev/null', 'w')

    p = subprocess.Popen(" ".join([needlecmd] + params), shell=True, bufsize=bufsize, stdin=ind, stdout=outd, stderr=errd, close_fds=True)
    p.wait()
    [w, r, e] = [p.stdin, p.stdout, p.stderr]
    
    print " ".join([needlecmd] + params)
    return r


def water(watercmd, aseqfile, bseqfile, **keywds):
    """  Standard (Mandatory) qualifiers:
  [-asequence]         sequence   Sequence filename and optional format, or
                                  reference (input USA)
  [-bsequence]         seqall     Sequence(s) filename and optional format, or
                                  reference (input USA)
   -gapopen            float      [10.0 for any sequence] The gap open penalty
                                  is the score taken away when a gap is
                                  created. The best value depends on the
                                  choice of comparison matrix. The default
                                  value assumes you are using the EBLOSUM62
                                  matrix for protein sequences, and the
                                  EDNAFULL matrix for nucleotide sequences.
                                  (Number from 0.000 to 100.000)
   -gapextend          float      [0.5 for any sequence] The gap extension
                                  penalty is added to the standard gap penalty
                                  for each base or residue in the gap. This
                                  is how long gaps are penalized. Usually you
                                  will expect a few long gaps rather than many
                                  short gaps, so the gap extension penalty
                                  should be lower than the gap penalty. An
                                  exception is where one or both sequences are
                                  single reads with possible sequencing
                                  errors in which case you would expect many
                                  single base gaps. You can get this result by
                                  setting the gap open penalty to zero (or
                                  very low) and using the gap extension
                                  penalty to control gap scoring. (Number from
                                  0.000 to 10.000)
  [-outfile]           align      [*.water] Output alignment file name

   Additional (Optional) qualifiers:
   -datafile           matrixf    [EBLOSUM62 for protein, EDNAFULL for DNA]
                                  This is the scoring matrix file used when
                                  comparing sequences. By default it is the
                                  file 'EBLOSUM62' (for proteins) or the file
                                  'EDNAFULL' (for nucleic sequences). These
                                  files are found in the 'data' directory of
                                  the EMBOSS installation.

   Advanced (Unprompted) qualifiers:
   -[no]brief          boolean    [Y] Brief identity and similarity

   Associated qualifiers:

   "-asequence" associated qualifiers
   -sbegin1            integer    Start of the sequence to be used
   -send1              integer    End of the sequence to be used
   -sreverse1          boolean    Reverse (if DNA)
   -sask1              boolean    Ask for begin/end/reverse
   -snucleotide1       boolean    Sequence is nucleotide
   -sprotein1          boolean    Sequence is protein
   -slower1            boolean    Make lower case
   -supper1            boolean    Make upper case
   -sformat1           string     Input sequence format
   -sdbname1           string     Database name
   -sid1               string     Entryname
   -ufo1               string     UFO features
   -fformat1           string     Features format
   -fopenfile1         string     Features file name

   "-bsequence" associated qualifiers
   -sbegin2            integer    Start of each sequence to be used
   -send2              integer    End of each sequence to be used
   -sreverse2          boolean    Reverse (if DNA)
   -sask2              boolean    Ask for begin/end/reverse
   -snucleotide2       boolean    Sequence is nucleotide
   -sprotein2          boolean    Sequence is protein
   -slower2            boolean    Make lower case
   -supper2            boolean    Make upper case
   -sformat2           string     Input sequence format
   -sdbname2           string     Database name
   -sid2               string     Entryname
   -ufo2               string     UFO features
   -fformat2           string     Features format
   -fopenfile2         string     Features file name

   "-outfile" associated qualifiers
   -aformat3           string     Alignment format
   -aextension3        string     File name extension
   -adirectory3        string     Output directory
   -aname3             string     Base file name
   -awidth3            integer    Alignment width
   -aaccshow3          boolean    Show accession number in the header
   -adesshow3          boolean    Show description in the header
   -ausashow3          boolean    Show the full USA in the alignment
   -aglobal3           boolean    Show the full sequence in alignment

   General qualifiers:
   -auto               boolean    Turn off prompts
   -stdout             boolean    Write standard output
   -filter             boolean    Read standard input, write standard output
   -options            boolean    Prompt for standard and additional values
   -debug              boolean    Write debug output to program.dbg
   -verbose            boolean    Report some/full command line options
   -help               boolean    Report command line options. More
                                  information on associated and general
                                  qualifiers can be found with -help -verbose
   -warning            boolean    Report warnings
   -error              boolean    Report errors
   -fatal              boolean    Report fatal errors
   -die                boolean    Report dying program messages

    """
    att2param = {
        'seq1'     : '-asequence',
        'seq2'     : '-bsequence',
        
        'gapopen'  : '-gapopen',
        'gapext'   : '-gapextend',
        'datafile' : '-datafile',
        
        'sb1'      : '-sbegin1',
        'se1'      : '-send1',

        'sb2'      : '-sbegin2',
        'se2'      : '-send2',
        
        'out'      : '-outfile',
        'stdout'   : '-stdout', # TODO, doesn't work yet
        'outform'  : '-aformat3',

        'usashow'  : '-usashow'
        # and more to come
        }

    if not os.path.exists(watercmd):
        raise ValueError, "water does not exist at %s" % watercmd
    
    params = []

    params.extend([att2param['seq1'], aseqfile])
    params.extend([att2param['seq2'], bseqfile])
    params.extend([att2param['gapopen'], str(10.0)])
    params.extend([att2param['gapext'], str(0.5)])
#    params.extend([att2param['out'], aseqfile.rstrip('.fas') + '.water'])
    params.extend([att2param['outform'], 'markx10'])

    for attr in keywds.keys():
        params.extend([att2param[attr], str(keywds[attr])])

    bufsize = -1
    ind  = None
    outd = None
    errd = open('/dev/null', 'w')

    p = subprocess.Popen(" ".join([watercmd] + params), shell=True, bufsize=bufsize, stdin=ind, stdout=outd, stderr=errd, close_fds=True)
    p.wait()
    [w, r, e] = [p.stdin, p.stdout, p.stderr]
    
    print " ".join([watercmd] + params)
    return r

def stretcher(stretchercmd, aseqfile, bseqfile, **keywds):
    """Standard (Mandatory) qualifiers:
    [-asequence]         sequence   Sequence filename and optional format, or
                                      reference (input USA)
    [-bsequence]         sequence   Sequence filename and optional format, or
                                  reference (input USA)
    [-outfile]           align      [*.stretcher] Output alignment file name

   Additional (Optional) qualifiers:
   -datafile           matrix     [EBLOSUM62 for protein, EDNAFULL for DNA]
                                  This is the scoring matrix file used when
                                  comparing sequences. By default it is the
                                  file 'EBLOSUM62' (for proteins) or the file
                                  'EDNAFULL' (for nucleic sequences). These
                                  files are found in the 'data' directory of
                                  the EMBOSS installation.
   -gapopen            integer    [12 for protein, 16 for nucleic] Gap penalty
                                  (Positive integer)
   -gapextend          integer    [2 for protein, 4 for nucleic] Gap length
                                  penalty (Positive integer)

   Advanced (Unprompted) qualifiers: (none)
   Associated qualifiers:

   "-asequence" associated qualifiers
   -sbegin1            integer    Start of the sequence to be used
   -send1              integer    End of the sequence to be used
   -sreverse1          boolean    Reverse (if DNA)
   -sask1              boolean    Ask for begin/end/reverse
   -snucleotide1       boolean    Sequence is nucleotide
   -sprotein1          boolean    Sequence is protein
   -slower1            boolean    Make lower case
   -supper1            boolean    Make upper case
   -sformat1           string     Input sequence format
   -sdbname1           string     Database name
   -sid1               string     Entryname
   -ufo1               string     UFO features
   -fformat1           string     Features format
   -fopenfile1         string     Features file name

   "-bsequence" associated qualifiers
   -sbegin2            integer    Start of the sequence to be used
   -send2              integer    End of the sequence to be used
   -sreverse2          boolean    Reverse (if DNA)
   -sask2              boolean    Ask for begin/end/reverse
   -snucleotide2       boolean    Sequence is nucleotide
   -sprotein2          boolean    Sequence is protein
   -slower2            boolean    Make lower case
   -supper2            boolean    Make upper case
   -sformat2           string     Input sequence format
   -sdbname2           string     Database name
   -sid2               string     Entryname
   -ufo2               string     UFO features
   -fformat2           string     Features format
   -fopenfile2         string     Features file name

   "-outfile" associated qualifiers
   -aformat3           string     Alignment format
   -aextension3        string     File name extension
   -adirectory3        string     Output directory
   -aname3             string     Base file name
   -awidth3            integer    Alignment width
   -aaccshow3          boolean    Show accession number in the header
   -adesshow3          boolean    Show description in the header
   -ausashow3          boolean    Show the full USA in the alignment
   -aglobal3           boolean    Show the full sequence in alignment

   General qualifiers:
   -auto               boolean    Turn off prompts
   -stdout             boolean    Write standard output
   -filter             boolean    Read standard input, write standard output
   -options            boolean    Prompt for standard and additional values
   -debug              boolean    Write debug output to program.dbg
   -verbose            boolean    Report some/full command line options
   -help               boolean    Report command line options. More
                                  information on associated and general
                                  qualifiers can be found with -help -verbose
   -warning            boolean    Report warnings
   -error              boolean    Report errors
   -fatal              boolean    Report fatal errors
   -die                boolean    Report dying program messages
    """

    att2param = {
        'seq1'     : '-asequence',
        'seq2'     : '-bsequence',
        
        'gapopen'  : '-gapopen',
        'gapext'   : '-gapextend',
        'datafile' : '-datafile',
        
        'sb1'      : '-sbegin1',
        'se1'      : '-send1',

        'sb2'      : '-sbegin2',
        'se2'      : '-send2',
        
        'out'      : '-outfile',
        'stdout'   : '-stdout', # TODO, doesn't work yet
        'outform'  : '-aformat3',

        'usashow'  : '-usashow'
        # and more to come
        }

    if not os.path.exists(stretchercmd):
        raise ValueError, "stretcher does not exist at %s" % stretchercmd
    
    params = []

    params.extend([att2param['seq1'], aseqfile])
    params.extend([att2param['seq2'], bseqfile])
    #    params.extend([att2param['gapopen'], str(16)])
    #    params.extend([att2param['gapext'], str(4)])
    params.extend([att2param['outform'], 'markx10'])

    try:
        keywds['stdout'] == True
        params.extend([att2param['stdout']])
        print 'going to stdout'
    except KeyError:
        if not keywds.has_key('out'):
            params.extend([att2param['out'], aseqfile.rstrip('.fas') + '.water'])
            print 'going to', aseqfile.rstrip('.fas') + '.water'
    

    for attr in keywds.keys():
        if attr != 'stdout':
            params.extend([att2param[attr], str(keywds[attr])])
    print params
    bufsize = -1
    ind  = None
    outd = None #open('./OUTSTREAM', 'w')
    errd = open('/dev/null', 'w')

    p = subprocess.Popen(" ".join([stretchercmd] + params), shell=True, bufsize=bufsize, stdin=ind, stdout=outd, stderr=errd, close_fds=True)
    #p = subprocess.Popen(" ".join([stretchercmd] + params), shell=True)#, bufsize=bufsize, stdin=ind, stdout=outd, stderr=errd, close_fds=True)
    w, r, e = (p.stdin, p.stdout, p.stderr)
    p.wait()
#    print " ".join([stretchercmd] + params)
    return r

if __name__ == '__main__':
    
    my_water_exe = '/Users/ozagordi/Work/tools/EMBOSS-5.0.0/emboss/water'   

    s1 = 'short1.fas' # asequence is the reference
    s2 = 'short2.fas' # bsequence are the reads
    go = 10.0
    ge = 0.5
    of = 'out'
    sb2 = 0
    se2 = 10
    af = 'markx10'

    water(my_water_exe, 'short1.fas', 'short2.fas', sb2=1, se2=0, out=of)
    handle = open(of)
    for alin in MyAlignIO.parse(handle, "markx10") :
            assert len(alin.get_all_seqs()) == 2, "Should be pairwise!"
            print '%s' % '-'*30
            print "Alignment length is %i" % alin.get_alignment_length()
            print 'all seq %s' % alin.get_all_seqs()
#            for k in alin._annotations:
#                print k, '=', alin._annotations[k]
#            print '%s' % '-'*30                
            for record in alin :
#                print '%s' % '-'*30
#                for k in record.annotations:
#                    print k,'=', record.annotations[k]
                print record.seq, record.name, record.id, record.annotations['al_start'], record.annotations['al_stop']
