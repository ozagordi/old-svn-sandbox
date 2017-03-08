#!/usr/bin/env python

import os
import sys
import logging
import logging.handlers
import socket
from Bio import SeqIO

homedir = os.path.expanduser('~/')
sys.path.append(homedir)

our_hosts = ['bs-mbp08', 'bs-dsvr07', 'bs-dsvr24']
hostname = socket.gethostname().split('.')[0]
min_length = 200
amb_thresh = 2

LOG_FILENAME = './decontaminate.log'
# Make a global logging object.
x_log = logging.getLogger('log_inst')
x_log.setLevel(logging.DEBUG)
# This handler writes everything to a file.
#h = logging.FileHandler('./log', 'w')
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter('%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s')
h.setFormatter(f)
x_log.addHandler(h)

def check_reference(fas_reference):
    ''' check the reference genome
    '''
    from Bio import SeqIO
    
    rh = open(fas_reference, 'rU')
    ref = list(SeqIO.parse(rh, 'fasta') )
    assert len(ref) == 1, 'One and only one sequence must be in the reference file'
    gen_length = len(ref[0].seq)
    if gen_length > 20:
        gen_ref_start = ref[0].seq.tostring()[:20]
    else:
        gen_ref_start = ref[0].seq.tostring()

    assert 'N' not in ref[0].seq, "Found an ambiguous position in the reference '%s'" % ref[0].id
    x_log.info('The reference genome length is %d' % gen_length )
    
    return gen_ref_start


def prepare_reads(f_fasta_filename):
    '''
    '''
    from Bio import SeqIO
    import math
    f_fasta = open(f_fasta_filename)
    seqlist_raw = list(SeqIO.parse(f_fasta, 'fasta'))
    
    seqlist = [ s for s in seqlist_raw if len(s) >= min_length ]
    countreads = len(seqlist)
    x_log.info('removed %d reads because shorter than %d' % (len(seqlist_raw) - countreads, min_length))
    
    # forward...
    f_fasta_forward_filename = 'tmp_reads_f.fas'
    f_fasta_forward = open(f_fasta_forward_filename, 'w')
    SeqIO.write(seqlist, f_fasta_forward, 'fasta')
    f_fasta_forward.close()
    
    # ...and reverse
    for seq in seqlist:
        seq.seq = seq.seq.reverse_complement()
    f_fasta.close()
    f_fasta_reverse_filename = 'tmp_reads_r.fas'
    f_fasta_reverse = open(f_fasta_reverse_filename,'w')
    SeqIO.write(seqlist,f_fasta_reverse,'fasta')
    f_fasta_reverse.close()
    
    x_log.info('Remaining %d reads', countreads)

    length = 0.
    length2 = 0.
    n = 0.

    for read in seqlist:
        if read.seq.count('N') < amb_thresh:
            len_seq = len(read.seq)
            length += float(len_seq)
            length2 += float(len_seq*len_seq)
            # readdict[read.] = [seq,len_seq]
            n += 1.
    meanlr = length/n
    stdlr = math.sqrt((n*length2 - length*length)/(n*n - n))
    #allowed_length = [meanlr - acclength * stdlr, meanlr + (1+acclength) * stdlr]
    x_log.info('mean length=%d, stderr=%d' % (meanlr, stdlr))
    
    return



def align_strand(al_info):
    '''
    '''
    from Bio.Emboss.Applications import NeedleCommandline
    import subprocess
    
    ref_file = al_info['ref_file']
    in_file = al_info['in_file']
    out_file = al_info['out_file']
    
    if os.path.exists(out_file):
        return
    
    cline = NeedleCommandline(gapopen=6.0, gapextend=3.0)
    cline.asequence = ref_file
    cline.bsequence = in_file
    cline.outfile = out_file
    cline.aformat = 'markx10'
    cml = str(cline) + ' -adesshow3 -auto'
    x_log.info(cml)
    try:
        retcode = subprocess.call(cml, shell=True)
        if retcode < 0:
            x_log.info('Child diri_sampler was terminated by signal %d' -retcode)
        else:
            x_log.info('Child diri_sampler returned %d' % retcode)
    except OSError, ee:
        x_log.exception('Execution of diri_sampler failed:' + ee)
    
    return

def parse_alignments(ref_name):
    '''
    '''
    
    from pythonlib.MarkxIO import Markx10Iterator

    count_forward = 0
    count_reverse = 0
    score = {}
    
    f_forward = open('tmp_align_f_%s.needle' % ref_name)
    f_reverse = open('tmp_align_r_%s.needle' % ref_name)
    
    forwardaligniter = Markx10Iterator(f_forward)
    reversealigniter = Markx10Iterator(f_reverse)
    
    x_log.info('parsing the alignments %s' % ref_name)
    
    # iterates through the alignments
    while True:
        
        try:
            f_align = forwardaligniter.next()
            r_align = reversealigniter.next()
        except:
            break
        if f_align is None or r_align is None:
            break
        
        assert f_align.get_all_seqs()[1].id == r_align.get_all_seqs()[1].id, 'same seq back and forward'
        descr = f_align.get_all_seqs()[1].id
        
        if float(f_align._annotations['sw_score']) > float(r_align._annotations['sw_score']):
            score[descr] = float(f_align._annotations['sw_score'])
            count_forward += 1
        else:
            score[descr] = float(r_align._annotations['sw_score'])
            count_reverse += 1
    
    return score, count_forward, count_reverse

def main():
    '''
    '''
    
    from multiprocessing import Pool
    import os.path
    pad_insert = False
    score = {}
    args = sys.argv

    x_log.info('\n\n')
    x_log.info('running with:')
    x_log.info(' '.join(sys.argv) + '\n')
    
    try:
        f_fasta_filename = args[1]
        references = args[2:]
    except:
        x_log.exception('wrong usage')
        sys.exit('usage: decontaminate.py reads.fas reference_1 reference_2...')

    gen_ref_start = [check_reference(ref) for ref in references]
    
    prepare_reads(f_fasta_filename)

    strands = []
    for ref in references:
        name = ref.split('/')[-1].split('.')[0]
        
        af = {
            'in_file': 'tmp_reads_f.fas',
            'ref_file': ref,
            'out_file': 'tmp_align_f_%s.needle' % name
            }
        strands.append(af)
        ar = {
            'in_file': 'tmp_reads_r.fas',
            'ref_file': ref,
            'out_file': 'tmp_align_r_%s.needle' % name
            }
        strands.append(ar)
    
    pool = Pool()
    pool.map(align_strand, strands)

    refs = [ i[1]['ref_file'].split('/')[-1].split('.')[0] for i in enumerate(strands) if i[0]%2 ]
    
    
    for r in refs:
        score[r], cf, cr = parse_alignments(r)
        x_log.info('forward=%d\treverse=%d\n' % (cf, cr))

    
    assignment = {}
    for read in score[refs[0]]:
        ts = [ (score[r][read], r) for r in refs ]
        tss = sorted(ts, reverse=True)
        try:
            assignment[tss[0][1]].append(read)
        except:
            assignment[tss[0][1]] = []
            assignment[tss[0][1]].append(read)

    
    read_dict = SeqIO.to_dict(SeqIO.parse(open(f_fasta_filename), 'fasta'))
    
    for ref in assignment:
        oh = open('closer_to_%s.fas' % ref, 'w')
        for read in assignment[ref]:
            or_read = read.split('#')[0]
            oh.write('>%s\n' % or_read)
            oh.write('%s\n' % read_dict[or_read].seq.tostring())
        oh.close()

if __name__ == '__main__':
    main()
