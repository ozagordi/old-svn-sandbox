#!/usr/bin/env python

#from __future__ import print_function
import pysam
assert pysam.__version__ == '0.4.1'

def my_fetch_callback( alignment ):
    print str(alignment)

def main():
    import sys
    import os.path
    import subprocess
    from textwrap import fill
    
    file_in = sys.argv[1]
    if not os.path.exists(file_in):
        sys.exit('file %s not found' % file_in)
    
    #file_stem = ''.join(os.path.split(file_in)[1].split('.')[:-2])
    if file_in.endswith('.sam'):
        sys.exit('input must be sorted and indexed BAM file')
    elif file_in.endswith('.bam'):
        samfile = pysam.Samfile(file_in, "rb" )
    else:
        sys.exit('file must have bam extension!')
    
    sam_ref = map(str, samfile.references)
    i = 0
    for ref_name in sam_ref:
        print ref_name
        i += 1
        h = open('../reads_per_ampl/reads_ampl_%d.fasta' % i, 'w')
        for j, read in enumerate(samfile.fetch(ref_name)):
            assert ref_name == samfile.references[read.tid] , ref_name + samfile.references[read.tid]
            h.write('>%s | ref %s pos %d qstart %d\n' % (read.qname, ref_name, read.pos, read.qstart))
            h.write(fill(read.seq, 80))
            h.write('\n')
        h.close()
        print j+1
        #break
        #print >> sys.stderr, count
        

if __name__ == '__main__':
    main()
