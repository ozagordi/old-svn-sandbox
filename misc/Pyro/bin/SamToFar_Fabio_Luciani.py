#!/usr/bin/python

import sys
import os
import getopt

def rmEnter(inFile):
    """Removes newlines from quality info in sam."""
    out = open(inFile+'b','w')
    input = open(inFile)
    first = input.readline()
    while first[0]=='@':
        first = input.readline()
    init = first[:6]
    input.close()
    input = open(inFile)
    for i in input:
        a=i.rstrip()
        if a[:6]!=init:
            out.write(a)
        else:
            out.write('\n'+a)
    input.close()
    out.close()
    return inFile+'b'

def splitCig(s):
    """Splits a cigar string into parts."""
    a=0
    i=0
    l=[]
    while i<len(s):
        if s[i].isdigit()!=True:
            l.append((int(s[a:i]),s[i]))
            a=i+1
        i+=1
    return l

def getSam(line):
    """Parses a sam read."""
    parts=line.split('\t')
    title=parts[0]
    pos = parts[3]
    cigar = parts[5]
    seq = parts[9]
    return title, int(pos), cigar, seq

def expCig(ref, pos, seq, cigL, begin, end):
    """Expands a cigar string to fasta format."""
    lSeqPos = 0
    rSeqPos = 0 
    lRefPos = pos-1
    rRefPos = pos-1
    newSeq=''
    for i in cigL:
        if i[1]=='D':
            rRefPos+=i[0]
            if int(i[0])%3!=0:
                newSeq+=ref[lRefPos:rRefPos]
            else:
                newSeq+=('-'*int(i[0]))
            lRefPos=rRefPos
        elif i[1]=='I':
            rSeqPos+=(i[0])
            lSeqPos+=(i[0])
        else:
            rSeqPos+=i[0]
            newSeq+=seq[lSeqPos:rSeqPos]
            lSeqPos=rSeqPos
            lRefPos+=i[0]
            rRefPos+=i[0]
    return newSeq

def SamToFar(InRef,InSam,OutFar, begin, end):
    """Converts a slice of sam file to a far file."""
    x=0
    ref = open(str(InRef),'r')
    sam = open(InSam,'r')
    out = open(str(OutFar),'w')
    refSeq = ''
    next=ref.readline()
    if next[0]=='>':
        next=ref.readline()
    while next:
        refSeq+=next.rstrip()
        next = ref.readline()
    next=sam.readline()
    while next[0]=='@':
        next = sam.readline()
    while next:
        title, pos, cigar, seq = getSam(next)
        if cigar!='*':
            x+=1
            cigL=splitCig(cigar)
            read = expCig(refSeq,pos,seq,cigL,begin,end)
            bp=begin-pos
            epe=end-(pos+len(read))+1
            ep=end-pos +1
            pb=pos-begin
            if (begin<=pos and pos<= end) or (begin<=pos+len(read) and pos+len(read)<=end):
                if bp>0:
                    read=read[bp:]+epe*'-'
                elif ep>0 and epe<0:
                    read=pb*'-'+read[:ep]
                else:
                    read=pb*'-'+read[:]+epe*'-' 			
                seg=0
                if ('A' in read) or ('T' in read) or ('C' in read) or ('G' in read):
                    s=read.rstrip('-')
		    s=s.lstrip('-')
		    if len(s)>=30:			
                        out.write('>read_'+title+'_'+str(x)+' '+title+'_'+str(x)+'\n')
                        while seg < (end-begin-79):
                            out.write(read[seg:(seg+80)]+'\n')
                            seg+=80
                        out.write(read[seg:(end-begin+1)]+'\n')     
        next = sam.readline()
    out.close()

def usage():
    """Describes how to run this file."""
    print "\n\nSamToFar.py:\n"
    print "Takes a slice from a .sam file and creates a"
    print ".far file for input into shorah. Any reads"
    print "starting within a given beginning and end"
    print "position (inclusive) are written to the .far file."
    print "\nArguments:\n"
    print "     -h: prints these instructions."
    print "     -s: .sam input file."
    print "     -r: reference file in fasta format."
    print "     -o: desired output file name."
    print "     -b: begin position for slice."
    print "     -e: end position for slice.\n\n"

def main(argv):
    """Parses keyboard input and runs program."""
    inSam=''
    inRef=''
    outFar=''
    beg=''
    end=''
    try:
        opts, args = getopt.getopt(argv, "hs:r:o:b:e:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-s':
            inSam=arg
        elif opt == '-r':
            inRef=arg
        elif opt == '-o':
            outFar=arg
        elif opt == '-b':
            beg=int(arg)
        elif opt == '-e':
            end=int(arg)
    if (inSam=='' or inRef=='' or outFar=='' or beg=='' or end==''):
        usage()
        sys.exit()
    input = rmEnter(inSam)
    SamToFar(inRef, input, outFar, beg, end)
    os.remove(input) 
    
if __name__=="__main__":
    main(sys.argv[1:])
