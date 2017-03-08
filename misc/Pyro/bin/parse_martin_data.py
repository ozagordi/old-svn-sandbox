#!/usr/bin/env python
import user
import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

'''
Parses Martin's data according to the barcode. See mail-martin-2009-02-05.txt
'''


#>FP3CAUT03DRHMI rank=0000006 x=1426.0 y=616.0 length=73
#ACGAGTGCGTTATCCTTAGCTTCCTCAATCACTCTTGGCAACGACCCTTAGTCACAATAA
#GATAGGGNGGACA

def make_record(record, ad1, desc=''):
    ret_rec = SeqRecord(seq=record.seq[ad1:], \
                            id = record.id,\
                            description = desc)
    return ret_rec

def make_rc_record(record, ad1, desc=''):
    ret_rec = SeqRecord(seq = record.seq.reverse_complement()[ad1:], \
                            id = record.id,\
                            description = desc)
    return ret_rec

args = sys.argv

try:
    sample = int(args[1])
except:
    sys.exit('usage: parse_martin_data.py sample_number')
if sample == 3:
    tag1 = 'ACGAGTGCGT'
    tag2 = 'ACGCTCGACA'
if sample == 4:
    tag1 = 'AGACGCACTC'
    tag2 = 'AGCACTGTAG'

length = []
after_len_a = []
after_len_b = []
home_dir = os.path.expanduser('~/')
h = open(os.path.join(home_dir ,'Work/Pyro/data/Martin_2009-02-05-raw/fna-files/%d.TCA.454Reads.fna' % sample))

tmp_h = open('tmp.fas', 'w')
for line in h:
    if line.startswith('>'):
        lsp = line.split()
        tmp_h.write(lsp[0] + '\n')
        r4 = lsp[4]
        length.append(int(r4.partition('=')[2]))
    else:
        tmp_h.write(line)

nseq = len(length)
print '    Found', nseq, 'reads'
tmp_o = open('tmp.fas')

t1 = []
t2 = []
array = [0]*5
count_tag = {}
thresh = 1 # allowed errors in the tags
seqs = SeqIO.parse(tmp_o, 'fasta')

while True:
    try:
        s = seqs.next()
    except:
        break
    
    check_tag = 0
    tst = s.seq.tostring()[:len(tag1)]
    tst_rev = s.seq.reverse_complement().tostring()[:len(tag1)]
    
    ed1 = 0
    ed2 = 0
    ed3 = 0
    ed4 = 0
    
    for i in range(len(tst)): #1; seq <-> tag1
        ed1 += (tst[i] != tag1[i])
        
    for i in range(len(tst_rev)): #2; seq_rev <-> tag1
        ed2 += (tst_rev[i] != tag1[i])

    for i in range(len(tst)): #3; seq <-> tag2
        ed3 += (tst[i] != tag2[i])
        
    for i in range(len(tst_rev)): #4; seq_rev <-> tag2
        ed4 += (tst_rev[i] != tag2[i])
    
    if ed1 <= thresh:
        t1.append(make_record(s, len(tag1)))
        check_tag = 1
        array[check_tag] += 1
        after_len_a.append(len(s))
        continue
    
    if ed2 <= thresh:
#        t1.append(make_rc_record(s, len(tag1)))
        check_tag = 2
        array[check_tag] += 1
        continue
    
    if ed3 <= thresh:
        t2.append(make_record(s, len(tag2)))
        check_tag = 3
        array[check_tag] += 1
        after_len_b.append(len(s))
        continue
    
    if ed4 <= thresh:
#        t2.append(make_rc_record(s, len(tag2)))
        check_tag = 4
        array[check_tag] += 1
        continue
    
    try:
        count_tag[tst] += 1
    except:
        count_tag[tst] = 1

print '   sample a, selected reads', len(after_len_a)
print '   sample b, selected reads', len(after_len_b)
#for k in count_tag:
#    if count_tag[k] > 100:
#        print k, count_tag[k]

if sample == 3:
    h31 = open(os.path.join(home_dir, 'Work/Pyro/data/Martin_2009-02-05-parsed/clone_mix/reads.fas'), 'w')
    h32 = open(os.path.join(home_dir, 'Work/Pyro/data/Martin_2009-02-05-parsed/PCR_clone_mix/reads.fas'), 'w')
    SeqIO.write(t1, h31, 'fasta')
    SeqIO.write(t2, h32, 'fasta')

if sample == 4:
    h41 = open(os.path.join(home_dir, 'Work/Pyro/data/Martin_2009-02-05-parsed/patient_1/reads.fas'), 'w')
    h42 = open(os.path.join(home_dir, 'Work/Pyro/data/Martin_2009-02-05-parsed/patient_2/reads.fas'), 'w')
    SeqIO.write(t1, h41, 'fasta')
    SeqIO.write(t2, h42, 'fasta')
    
try:
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
except:
    print 'exit'
    sys.exit()

#plt.subplots_adjust(hspace=0.1)

#plt.grid(True)

# first subplot
plt.subplot(211)
n, bins, patches = plt.hist(after_len_a, 20, normed=0, facecolor='green', alpha=0.75)
plt.xlabel(r'read length distribution')

# second subplot
plt.subplot(212)
n2, bins2, patches2 = plt.hist(after_len_b, 20, normed=0, facecolor='green', alpha=0.75)

plt.xlabel(r'read length distribution')

#plt.xlabel('cycle')
#plt.ylabel('entropy')

#plt.xlabel('time')
#plt.ylabel('times')
#xt =  ('7:30', '8:00', '8:30', '9:00', '9:30', '10:00', '10:30', '11:00')
imtype = 'pdf'
plt.savefig(os.path.join(home_dir, 'Work/Pyro/data/Martin_2009-02-05-parsed/fig%d.%s' % (sample, imtype)),
            dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=imtype,
            transparent=False)

plt.show()
