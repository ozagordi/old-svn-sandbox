#!/usr/bin/env python
'''mainly used to submit amplicon jobs to the grid'''

import socket
import logging
import logging.handlers

hostname = socket.gethostname().split('.')[0]
LOG_FILENAME = './submission.log'
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
    
def compile_string(exe=None, exe_args=None, slot_range=8, wd=None, err_file=None, binary=False, name='noname'):
    
    queue = None#'regular.q'
    
    resource = 'background'

    cml = 'qsub'
    #cml += ' -pe make %d' % slot_range
    cml += ' -m e -M ozagordi@bsse.ethz.ch'
    #cml += ' -v OMP_NUM_THREADS' # variable list
    #cml += ' -o %s' % out_file
    if err_file: cml += ' -e %s' % err_file
    if resource: cml += ' -l %s' % resource
    elif queue: cml += ' -q %s' % queue
    cml += ' -V' # export all env variables
    if wd: cml += ' -wd %s' % wd
    cml += ' -N %s' % name
    if binary: cml += ' -b y' #command as binary
    cml += ' %s ' % exe
    cml += exe_args
    return cml

def main():
    import sys
    import os
    
    args = sys.argv
    job_type = args[1]
    amplicon = int(args[2])
    analyses = ['s2f', 'diri', 'recco', 'error', 'sim_rec', 'recco_sim', 'pr_rec']
    homedir = os.path.expanduser('~/')
    if job_type not in analyses: sys.exit('only possible one of: [ ' + ' | '.join(analyses) + ' ]')
    
    if job_type == 'pr_rec':
        #rd = os.path.join(homedir, 'Work/HIVNGS/experiments/2011-05-19-recombination_run2_again/')
        rd = os.path.join(homedir, 'Work/HIVNGS/experiments/2011-06-28-recombination_run_5_new_ref/')
        os.chdir(rd)
        amp_dir = 'PR_%d' % amplicon
        try:
            os.mkdir(amp_dir)
        except:
            pass
        os.chdir(amp_dir)
        clones_file = os.path.join(homedir, 'References/viral_mix_5_seqs_amplicon.fasta')
        reads_file = '/home/ozagordi/Work/HIVNGS/data/PR_31-PR_36/PR_%d.far' % amplicon
        exe = '/home/ozagordi/server_home/Work/Pyro/bin/run_recco.py'
        exe_args = ' -r %s -l %s' % (reads_file, clones_file)
        cml = compile_string(exe=exe, exe_args=exe_args, err_file='queue_err', wd=os.path.join(rd, amp_dir), name='recco-PR_%d' %  amplicon)
        #print cml
        generic_application(cml)
    
    if job_type == 'recco_sim':
        rd = os.path.join(homedir, 'Work/HIVNGS/experiments/2011-05-16-sim_rec_per_ampl-rec_0.0/')
        os.chdir(rd)
        amp_dir = 'amplicon_%2.2d' % amplicon
        os.chdir(amp_dir)
        clones_file = 'FLS1_all_clones_ampl_%2.2d.fasta' % amplicon
        reads_file = 'sim.fasta'
        exe = '/home/ozagordi/server_home/Work/Pyro/bin/run_recco.py'
        exe_args = ' -r %s -l %s' % (reads_file, clones_file)
        cml = compile_string(exe=exe, exe_args=exe_args, err_file='queue_err', wd=os.path.join(rd, amp_dir), name='rcs-a-%d' %  amplicon)
        #print cml
        generic_application(cml)
    
    if job_type == 'recco':
        clones_file = '/home/ozagordi/References/HIV1_clones.fasta'
        reads_file = '/home/ozagordi/Work/HIVNGS/experiments/2011-04-11-separate_amplicons/reads_per_ampl/reads_ampl_%d.fasta' % amplicon
        rd = os.path.join(homedir, 'Work/HIVNGS/experiments/2011-04-18-recco_per_amplicon')
        os.chdir(rd)
        amp_dir = 'amplicon_%2.2d' % amplicon
        try:
            os.mkdir(amp_dir)
            logfun.info('creating dir for amplicon %d' % amplicon)
        except OSError:
            logfun.info('directory exists')
        os.chdir(amp_dir)
        
        exe = '/home/ozagordi/server_home/Work/Pyro/bin/run_recco.py'
        exe_args = ' -r %s -l %s' % (reads_file, clones_file)
        cml = compile_string(exe=exe, exe_args=exe_args, wd=os.path.join(rd, amp_dir), name='%s-a-%d' % (job_type, amplicon))
        #print cml
        generic_application(cml)

    if job_type == 's2f':
        from Bio import SeqIO
        all_ampl = os.path.join(homedir, 'server_home/References/FLS1_all_ampl.fasta')
        ampl_recs = SeqIO.parse(all_ampl, 'fasta')
        for sr in ampl_recs:
            if sr.id.startswith('FLS1_amplicon_%2.2d' % amplicon):
                ref_rec = sr
        fileref = os.path.join(homedir, 'server_home/References/HIV1-HXB2-RNA.fasta')
        reads_file = '/home/ozagordi/server_home/Work/HIVNGS/experiments/2011-04-11-separate_amplicons/reads_per_ampl/reads_ampl_%d.fasta' % amplicon
        rd = os.path.join(homedir, 'server_home/Work/HIVNGS/experiments/2011-04-19-diri_per_ampl_all_reads') 
        os.chdir(rd)
        amp_dir = 'amplicon_%2.2d' % amplicon
        try:
            os.mkdir(amp_dir)
            logfun.info('creating dir for amplicon %d' % amplicon)
        except OSError:
            logfun.info('directory exists')
        os.chdir(amp_dir)
        ref_ampl = 'ref_ampl_%2.2d.fasta' % amplicon
        SeqIO.write(ref_rec, ref_ampl, 'fasta')
        exe = os.path.join(homedir, 'server_home/bin/shorah/trunk/s2f.py')
        exe_args = ' -r %s -f %s -o reads_ampl_%2.2d.far' % (ref_ampl, reads_file, amplicon)
        cml = compile_string(exe=exe, exe_args=exe_args, slot_range=2, wd=os.path.join(rd, amp_dir), name='%s-a-%d' % (job_type, amplicon))
        #print cml
        generic_application(cml)
        
    if job_type == 'diri':
        n_iter = 10000
        alpha = 10.0
        history = 2000
        
        amp_dir = 'amplicon_%2.2d' % amplicon
        rd = os.path.join(homedir, 'server_home/Work/HIVNGS/experiments/2011-04-19-diri_per_ampl_all_reads', amp_dir) 
        os.chdir(rd)
        
        #cuts the high coverage part
        exe = os.path.join(homedir, 'server_home/Work/HIVNGS/bin/cut_high_coverage.py')
        cml = exe + ' reads_ampl_%2.2d.far' % amplicon
        generic_application(cml)
        
        #runs diri_sampler
        exe = os.path.join(homedir, 'server_home/bin/shorah/trunk/diri_sampler_grid')
        exe_args = ' -i reads_ampl_%2.2d_high_cov.far -j %d -a %f -K 20 -t %d' % (amplicon, n_iter, alpha, history)
        cml = compile_string(exe=exe, exe_args=exe_args, slot_range=1, binary=True, wd=rd, name='%s-a-%d' % (job_type, amplicon))
        #print cml
        generic_application(cml)


    if job_type == 'sim_rec':
        if hostname == 'bs-submit01':
            rd = os.path.join(homedir, 'server_home/Work/HIVNGS/experiments/2011-05-16-sim_rec_per_ampl-rec_0.0')
        else:
            rd = os.path.join(homedir, 'Work/HIVNGS/experiments/2011-05-16-sim_rec_per_ampl-rec_0.0')
        os.chdir(rd)
        amp_dir = 'amplicon_%2.2d' % amplicon
        try:
            os.mkdir(amp_dir)
            logfun.info('creating dir for amplicon %d' % amplicon)
        except OSError:
            logfun.info('directory exists')
        os.chdir(amp_dir)
        
        #extract all amplicons
#        if hostname == 'bs-submit01':
#            exe = os.path.join(homedir, 'server_home/Work/HIVNGS/bin/extract_all_amplicons.py')
#        else:
#            exe = os.path.join(homedir, 'Work/HIVNGS/bin/extract_all_amplicons.py')
#        cml = exe + ' %d' % amplicon
        #print cml
#        generic_application(cml)
        
        #simulate recombination
        if hostname == 'bs-submit01':
            exe = os.path.join(homedir, 'server_home/Work/Pyro/bin/simulate_recomb.py')
        else:
            exe = os.path.join(homedir, 'Work/Pyro/bin/simulate_recomb.py')
        cml = exe + ' FLS1_all_clones_ampl_%2.2d.fasta' % amplicon
        #print cml
        generic_application(cml)
        
    
    
if __name__ == '__main__':
    main()
