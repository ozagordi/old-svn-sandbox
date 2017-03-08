#!/usr/bin/env python

import sys
import os

import tempfile
import pickle
import socket
from Bio import SeqIO

import logging
import logging.handlers

reclog = logging.getLogger("reccolog") # Make a global logging object.
reclog.setLevel(logging.DEBUG) # set logging level
# This handler writes everything to a file.
LOG_FILENAME = './recco.log'
h = logging.handlers.RotatingFileHandler(LOG_FILENAME, 'w', maxBytes=100000, backupCount=5)
f = logging.Formatter("%(levelname)s\t%(asctime)s\n\t%(message)s\n")
h.setFormatter(f)
reclog.addHandler(h)

savings_thresh = 5
results = {}

homedir = os.path.expanduser('~/')
sys.path.append(homedir)
hostname = socket.gethostname().split('.')[0]
if hostname == 'bs-mbp115':
    clustalw_bin = '/opt/local/sbin/clustalw2'
    muscle_exe = '/Users/ozagordi/Work/tools/muscle'
    recco_jar = '/Users/ozagordi/Work/Pyro/bin/recco/Recco.jar'
    needle_exe = '/usr/local/bin/needle'
    water_exe = '/usr/local/bin/water'
    nucmx = '/Users/ozagordi/Work/tools/nucmx'
else:
    clustalw_bin = '/usr/local/bewi/clustalw-2.0.12/clustalw-2.0.12-linux-i686-libcppstatic/clustalw2'
    recco_jar = '/nas/ozagordi/bin/recco/Recco.jar'
    needle_exe = 'needle'
    muscle_exe = 'muscle'
    nucmx = '/home/ozagordi/Work/tools/nucmx'
    tempfile.tempdir = '/local0/tmp/'
        
aligner = 'muscle'

def parse_com_line():
    from optparse import OptionParser
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage = usage)
    parser.add_option("-r", "--reads_file", dest='r', type="string",
                      metavar="R_FILE", help="reads are contained in R_FILE"),
    #parser.add_option("-c", "--cons_file", dest='c', type="string",
    #                  metavar="C_FILE", help="consensus sequence is in C_FILE"),
    parser.add_option("-l", "--clones_file", dest='l', type="string",
                      metavar="CL_FILE", help="clones (haplotypes) are contained in CL_FILE"),
    parser.add_option("-m", "--min_length", dest='m', type="int",
                      metavar="MIN_LENGTH", help="min length of reads to analyse with recco"),
    parser.add_option("-v", "--verbose", help="verbose behaviour [default]",
                      action="store_true", dest="verbose", default=False)
    
    (options, args) = parser.parse_args()
    
    return options, args


def align_read2clones(scc):
    '''
    '''
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import AlignIO
    from pythonlib.MarkxIO import Markx10Iterator
    from pythonlib import Alignment
    import subprocess
    import hashlib
    import copy
    
    delta_start = 50
    delta_stop = 50
    
    seq_inst, cl_rec = scc
    m = hashlib.md5()
    m.update(seq_inst.id)
    seq_id = m.hexdigest()
    
    # read and clones
    all_recs = copy.copy(cl_rec)
    all_recs.append(SeqRecord(Seq(seq_inst.seq.tostring().upper().strip('N')), id=seq_id, description=seq_inst.id))
    n_seqs = len(all_recs)
    assert n_seqs == len(cl_rec) + 1
    handle = tempfile.NamedTemporaryFile(mode='w', delete=False)    
    SeqIO.write(all_recs, handle.name, 'fasta')
    handle.close()
    
    r_file_name = handle.name
    out_clustal = tempfile.NamedTemporaryFile(mode='w', delete=False)
    align_out = out_clustal.name
    
    if aligner == 'clustalw':
        # run clustalw
        com_line = '%s -INFILE=%s -OUTFILE=%s -output=PIR &> /dev/null' % (clustalw_bin, tmp_align_name, align_out)
        child_process = subprocess.call(str(com_line), shell=True)
        try:
            os.unlink(tmp_align_name+'.dnd')
        except:
            reclog.warning(tmp_align_name+'.dnd not found')
        try:
            os.unlink(tmp_align_name)
        except:
            reclog.warning(tmp_align_name+' not found either')
        
        #convert 'tmp_align.pir' to 'read_clones_al.fas'
        input_handle = open(align_out , "rU")
        #output_handle = open("read_clones_al.fas", "w")
        output_handle = tempfile.NamedTemporaryFile(mode='w', delete=False)
        read_clones_al = output_handle.name
        sequences = SeqIO.parse(input_handle, 'pir')
        count = SeqIO.write(sequences, output_handle, 'fasta')
        output_handle.close()
        input_handle.close()
        os.unlink(align_out)
        #END of convert
        return seq_inst, read_clones_al, seq_id
        
    elif aligner == 'muscle':
        #com_line = '%s -maxiters 1 -matrix %s -seqtype protein -in %s -out %s &> /dev/null' % (muscle_exe, nucmx, r_file_name, align_out)
        com_line = '%s -maxiters 1 -in %s -out %s &> /dev/null' % (muscle_exe, r_file_name, align_out)
        child_process = subprocess.call(str(com_line), shell=True)
        alignment = AlignIO.read(align_out, 'fasta')
        os.unlink(align_out)
        for s in alignment._records:
            if s.id == seq_id:
                read_start = len(s.seq) - len(s.seq.lstrip('-'))
                read_stop = len(s.seq.rstrip('-'))
                assert s[read_start] != '-'
                if read_start > 0: assert s[read_start-1] == '-'
                assert s[read_stop-1] != '-'
                if read_stop < len(s): assert s[read_stop] == '-'
        
        new_al = alignment[:, read_start:read_stop-1]
        SeqIO.write(new_al, align_out, 'fasta')
        
        return seq_inst, align_out, seq_id
    

def recomb_analysis(rc):
    '''
    This uses recco
    '''
    import subprocess
    read, read_clones_al, read_id = rc
    nperm=100
    read_id = read_id.split('#')[0]
    recco_out = tempfile.NamedTemporaryFile(mode='w', delete=False)
    recco_out.close()
    
    recco_rout = tempfile.NamedTemporaryFile(mode='w', delete=False)
    recco_rout.close()
    
    com_line = 'java -jar %s -c DNA -nperm %d -rec %s -out %s -rout %s -o -v %s' % \
               (recco_jar, nperm , read_id, recco_out.name, recco_rout.name, read_clones_al)
    
    child_process = subprocess.call(str(com_line), shell=True)
    os.unlink(read_clones_al)
    
    return recco_out.name, recco_rout.name


def parse_recco_output(ror):
    read, read_id, recco_out, recco_rout = ror
    oh = open(recco_out, "r")
    flines = oh.readlines()
    
    i = 0
    recco_dict = {}
    for line in flines:
        i += 1
        if line.startswith("END SequencePValues"):
           PVLine = flines[i-2].split()
           id_line = flines[i-3].split()
           break
        
    for i, j in enumerate(id_line):
        if j == read_id:
            rpv = PVLine[i]
            
    rh = open(recco_rout)
    
    
    idline = rh.next()
    ids = idline.replace(' pv', '-pv').replace('n C', 'n-C').replace('Seq ', 'Seq-').split()
    max_sav = -1000
    while True:
        try:
            values = rh.next().split()
        except:
            break
        assert values[0] == read_id.split('#')[0], '%s %s' % (values[0], read_id)
        this_results = {}
        for k, v in zip(ids[1:], values[1:]):
            this_results[k] = float(v)
        if this_results['Savings'] > max_sav:
            max_sav = this_results['Savings']
            res2ret = this_results    
        
        #reclog.info('rout file %s has no results' % rh.name)
        #print recco_rout, recco_out
    rh.close()
    oh.close()
    os.unlink(recco_out)
    os.unlink(recco_rout)
    
    if max_sav == -1000:
        return read_id, read.seq.tostring().upper().strip('N'), {'Savings': -1000.0}
    else:
        return read_id, read.seq.tostring().upper().strip('N'), res2ret



def main():
    '''
    '''
    from multiprocessing import Pool
    from textwrap import fill
    
    n_processes = 7
    print >> sys.stderr, "\n" + "-" * 55
    print >> sys.stderr, " WARNING: READS MUST BE IN THE CORRECT ORIENTATION"
    print >> sys.stderr, " THE PROGRAM WILL NOT TRY TO FIND THE BEST ORIENTATION"
    print >> sys.stderr, "-" * 55 + "\n"
    min_length = 200
    sav_thresh = 4
    options, args = parse_com_line()
    reclog.info(options)
    reclog.info(' '.join(sys.argv))
    in_stem = os.path.split(options.r)[1].split('.')[0]
    pv_pck = in_stem + '.pck'
    
    if not os.path.exists(pv_pck):
        clones_rec = list(SeqIO.parse(open(options.l), 'fasta'))
        # loop on all reads in reads_file (options.r)
        reads_to_run = [(s, clones_rec) for s in SeqIO.parse(open(options.r), 'fasta') if len(s) >= min_length]#[:20]
        print >> sys.stderr, 'Reads to run are', len(reads_to_run)
        reclog.info('Reads to run are %d' % len(reads_to_run))
        
        reclog.info('Now the MSA')
        pool = Pool(processes=n_processes)
        al_reads = pool.map(align_read2clones, reads_to_run)
        reclog.info('MSA have been computed')
        # align_read2clones returns read, read_clones_al, read_id_hashed
        print >> sys.stderr, al_reads[-1]
        
        reclog.info('Now the recombination analysis')
        pool = Pool(processes=n_processes)
        or_files = pool.map(recomb_analysis, al_reads)
        reclog.info('Recombination analysis completed')
        print >> sys.stderr, or_files[-1]
         
        ror = [(a[0][0], a[0][2], a[1][0], a[1][1]) for a in zip(al_reads, or_files)]
        reclog.info('Parsing recco output')
        recco_out = map(parse_recco_output, ror)
        print >> sys.stderr, recco_out[-1]
        
        h = open(pv_pck, 'w')
        pickle.dump(recco_out, h)
        h.close()
    else:
        recco_out = pickle.load(open(pv_pck))
        
    print >> sys.stderr, '%d reads analysed' % len(recco_out)
    
    pot_rec = 0
    count_sav = {}
    count = 0
    hdict = {}
    for ar in recco_out:
        sav_here = int(round(ar[2]['Savings']))
        try:
            hdict[sav_here].write('>%s\n' % ar[0])
            hdict[sav_here].write('%s\n' % fill(ar[1], 80))
        except KeyError:
            hdict[sav_here] = open('tmp-%d-savings.fasta' % sav_here, 'w')
            hdict[sav_here].write('>%s\n' % ar[0])
            hdict[sav_here].write('%s\n' % fill(ar[1], 80))
        if sav_here - int(ar[2]['Savings']) != 0.0:
            count += 1
        count_sav[sav_here] = count_sav.get(sav_here, 0) + 1
    #print 'sum is', sum(count_sav.values())
    #print sorted(count_sav.keys())
    print >> sys.stderr, count, 'reads have non-integer Savings'
    reclog.info('%d reads have non-integer Savings' % count)
    maxsav = int(max(count_sav.keys()))
    print "\"%s\"" % in_stem
    print '\"savings\", \"reads\", \"%\"'
    for sav in range(maxsav+1):
        try:
            print '%d, %d, %4.2f' % (sav, count_sav[sav], 100.0*count_sav[sav]/len(recco_out))
            #print >> sys.stderr, '-> %4.2f %% of the total\n' % (100.0*count_sav[sav]/len(recco_out))
        except KeyError:
            print '%d, 0, 0.0' % sav
            #print >> sys.stderr, 'no reads have %d savings\n' % (sav)

if __name__ == '__main__':
    main()
