#!/usr/bin/env python

from __future__ import print_function
import sys
import pysam
assert pysam.__version__ == '0.4.1'

def main():
    import sys
    import os.path
    import subprocess

    from Bio import SeqIO
    
    file_in = sys.argv[1]
    ref_name = sys.argv[2]
    start, end = sys.argv[3:5]
    if not os.path.exists(file_in):
        sys.exit('file %s not found' % file_in)
    #ref_file = sys.argv[2]
    file_stem = os.path.split(file_in)[1].split('.')[0]
    samfile = pysam.Samfile(file_in, "r" )
    #SAM -> BAM if not there
    if not os.path.exists(file_stem+'.bam'):
        print('SAM -> BAM', file=sys.stderr)
        cml = 'samtools view -b -S -o %s %s' % (file_stem+'.bam', file_in)
        print(cml, file=sys.stderr)
        subprocess.call(cml, shell='/bin/bash')
    else:
        print('BAM file was there, using it', file=sys.stderr)
        
    #sort BAM file
    if not os.path.exists(file_stem+'_sorted.bam'):
        print('Sorting', file=sys.stderr)
        pysam.sort(file_stem+'.bam', file_stem+'_sorted')
    else:
        print('BAM file already sorted, using it', file=sys.stderr)
        
    #index the sorted bam in any case    
    pysam.index(file_stem+'_sorted.bam')
    
    sorted_bam = pysam.Samfile(file_stem+'_sorted.bam', 'rb')
    #fastafile = pysam.Fastafile(ref_file)
    ref_name = sorted_bam.header['SQ'][0]['SN']
    print(ref_name, 'was given', file=sys.stderr)
    print(sorted_bam.header['SQ'][0]['SN'], 'was found', file=sys.stderr)
    
    al_reads = [a for a in sorted_bam.fetch( ref_name, 6815, 7367)]
    print(len(al_reads))
    print(al_reads[0])
            
if __name__ == '__main__':
    main()
