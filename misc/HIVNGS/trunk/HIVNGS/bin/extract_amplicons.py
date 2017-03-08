#!/usr/bin/env python
'''See also https://wiki-bsse.ethz.ch/display/HIVNGS/References
'''
import sys
import csv
import os
import socket
from Bio import SeqIO

hostname = socket.gethostname().split('.')[0]
if hostname == 'bs-submit01':
    al_exe = '/usr/local/bsse/EMBOSS/bin/needle'
else:
    al_exe = '/usr/local/bsse/EMBOSS/bin/needle'

homedir = os.path.expanduser('~/')
filein = os.path.join(homedir, 'Work/HIVNGS/data/FLS1_amplicons.csv')
#hxb2file = os.path.join(homedir, 'References/HIV-HXB2.fasta')
hxb2file = os.path.join(homedir, 'References/HIV1-HXB2-RNA.fasta')
print hxb2file
try:
    h = open(hxb2file)
except IOError:
    hxb2file = os.path.join(homedir, 'server_home/References/HIV1-HXB2-RNA.fasta')
    h = open(hxb2file)
HXB2_seq = list(SeqIO.parse(h, 'fasta'))[0]
h.close()

def extract_amplicon(refs, clone_file=hxb2file):
    import subprocess
    from Bio import AlignIO
    n, fw_primer, rev_primer = refs
    
    cml = '%s %s asis:%s -auto -aformat fasta -outfile ppp' % (al_exe, clone_file, fw_primer)
    p = subprocess.call(cml, shell='/bin/bash')
    al = AlignIO.read('ppp', 'fasta')
    os.unlink('ppp')
    seq_a, seq_b = al[0].seq, al[1].seq
    f_ampl_start = len(seq_b) - len(seq_b.lstrip('-'))
    f_ampl_stop = len(seq_b.rstrip('-'))
    assert seq_b[f_ampl_start] != '-'
    if f_ampl_start > 0: assert seq_b[f_ampl_start-1] == '-'
    assert seq_b[f_ampl_stop-1] != '-'
    if f_ampl_stop < len(seq_b): assert seq_b[f_ampl_stop] == '-'
    
    cml = '%s %s asis:%s -auto -aformat fasta -outfile ppp' % (al_exe, clone_file, rev_primer)
    p = subprocess.call(cml, shell='/bin/bash')
    al = AlignIO.read('ppp', 'fasta')
    os.unlink('ppp')
    seq_a, seq_b = al[0].seq, al[1].seq
    r_ampl_start = len(seq_b) - len(seq_b.lstrip('-'))
    r_ampl_stop = len(seq_b.rstrip('-'))
    assert seq_b[r_ampl_start] != '-'
    if r_ampl_start > 0: assert seq_b[r_ampl_start-1] == '-'
    assert seq_b[r_ampl_stop-1] != '-'
    if r_ampl_stop < len(seq_b): assert seq_b[r_ampl_stop] == '-'
    
    ampl_coord = f_ampl_start, f_ampl_stop, r_ampl_start, r_ampl_stop
    return ampl_coord


if __name__ == '__main__':
    from Bio.SeqRecord import SeqRecord
    ampl_list = []
    myreader = csv.reader(open(filein, 'rU'), delimiter=',')
    myreader.next()
    for row in myreader:
        ampl_coord = extract_amplicon((row[0], row[2], row[4]))
        f_ampl_start, f_ampl_stop, r_ampl_start, r_ampl_stop = ampl_coord
        this_rec = SeqRecord(HXB2_seq.seq[f_ampl_start:r_ampl_stop],
                             id='FLS1_amplicon_%2.2d-%d-%d' % (int(row[0]), f_ampl_start+1, r_ampl_stop),
                             description='| HIV RNA amplicon from FLS1')
        ampl_list.append(this_rec)
        if len(this_rec) == 0:
            print row[0], f_ampl_start, r_ampl_stop
    SeqIO.write(ampl_list, 'FLS1_all_ampl.fasta', 'fasta')
    print 'Done'
