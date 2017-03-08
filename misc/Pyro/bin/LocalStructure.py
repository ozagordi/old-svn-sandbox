#!/usr/bin/env python

from __future__ import print_function
import sys
import os.path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

hdir = os.path.expanduser('~/')#/pythonlib')
sys.path.append(hdir)
from pythonlib.HXB2_data import gene_coord

def str_num_compare(x, y):
    return int(x) - int(y)
    
class LocalStructure:
    '''This class    
    '''
    
    def __init__(self, sup_file, ref_file='References/HIV-HXB2.fasta', gene='protease'):
        '''
        '''
        import os.path
        
        full_ref = os.path.join(hdir, ref_file)
        full_sup = os.path.join(hdir, sup_file)
        
        self.sup_file = full_sup
        self.ref_file = full_ref
        
        s_head, self.name = os.path.split(full_sup)
        h = open(full_sup)
        self.seq_obj = SeqIO.parse(h, 'fasta')
        start, stop = gene_coord[gene]
        try:
            h = open(full_ref)
            sr = list(SeqIO.parse(h, 'fasta'))[0].seq.tostring()[start:stop].upper()
            self.ref = Seq(sr, IUPAC.unambiguous_dna)
        except IOError:
            self.ref = None
        
        self.cons = Seq(self.get_cons(), IUPAC.unambiguous_dna)
        
        self.frame = None
        # self.offset = None
        self.get_offset()
        sup_far = os.path.join(s_head, self.name.split('-')[0] + '.far')
        ns = SeqIO.parse(open(sup_far), 'fasta')
        self.n_reads = len([s for s in ns])
        
        self.dna_seqs = []
        self.prot_seqs = []
        
    def alignedvariants(self, threshold=0.9):
        import subprocess
        import re
        import itertools
        import hashlib
        from Bio.Emboss.Applications import NeedleCommandline
        from pythonlib import Alignment
        
        files = []
        var_dict = {}
        for i, s in enumerate(self.seq_obj):
            m_obj = re.search('posterior=(.*)\s*ave_reads=(.*)', s.description)
            post, ave_reads = map(float, (m_obj.group(1), m_obj.group(2)))
            if post < threshold or ave_reads < 1.:
                continue
            if post > 1.0:
                print('WARNING: posterior=', post, file=sys.stderr)
            outfile = 'tmp%d.needle' % i
            files.append(outfile)
            needle_cline = NeedleCommandline(asequence='asis:%s' % self.ref, bsequence='asis:%s' % s.seq.tostring().strip('-'), \
                                   outfile=outfile, gapopen=10.0, gapextend=0.5, aformat='markx10')
            needle_cline.auto = True
            
            try:
                retcode = subprocess.call(str(needle_cline), shell=True)
                if retcode < 0:
                    sys.exit('Child needle was terminated by signal %d' % -retcode)
 #               else:
 #                   print >> sys.stderr, 'Child needle returned %i' % retcode
            except OSError:
                sys.exit('Execution of needle failed: %s' % ee)
                pass
            
            tal = Alignment.alignfile2dict([outfile], 'support_seqs%d' % i, 10.0, 0.5, Verbose=False)
            os.remove(outfile)
            ka = tal.keys()[0]
            this = tal[ka]['asis']
            it_pair = itertools.izip(this.seq_a, this.seq_b)
            #this.summary()
            #start, stop = this.start, this.stop
            #it_pair = itertools.izip(this.seq_a[start-1:stop], this.seq_b[start-1:stop])
            
            this_seq = []
            while True:
                try:
                    p = it_pair.next()
                except StopIteration:
                    break
                if p is None:
                    break
                if p[1] == '-':
                    assert p[0] != '-', 'gap-gap?'
                    this_seq.append(p[0])
                elif p[0] != '-':
                    this_seq.append(p[1])
            ws = ''.join(this_seq)
            var_dict[ws] = var_dict.get(ws, 0) + ave_reads
            
        for k, v in var_dict.items():
            ts = Seq(k, IUPAC.unambiguous_dna)
            tsr = SeqRecord(ts, id = hashlib.sha224(k).hexdigest(), \
                            name='Reconstructed local hap')
            tsr.description = 'ave_reads=%f' % v
            self.dna_seqs.append(tsr)
        print('%d haplotypes have support >=%f'\
              % (len(files), threshold), file=sys.stderr)
        return self.dna_seqs
    

    def get_cons(self, plurality=0.1):
        '''Consensus by running EMBOSS cons
        '''
        import subprocess
        import os
        import itertools
        from pythonlib import Alignment
        
        cline = 'cons -sequence %s -stdout -auto' % self.sup_file
        cline += ' -plurality %f' % plurality
        
        p = subprocess.Popen(cline, shell=True, bufsize=1024, \
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
                             close_fds=True)
        
        sc =  list(SeqIO.parse(p.stdout, 'fasta'))[0].seq.tostring().upper()
        strcons = sc.replace('N', '')
        
        outfile = 'tmp.tmp'
        Alignment.needle_align(self.ref_file, 'asis:%s' % strcons, \
                                   outfile, go=10.0, ge=0.5)
        tal = Alignment.alignfile2dict([outfile], 'ref_cons_alignment', 10.0, 0.5)
        os.remove(outfile)
        ka = tal.keys()[0]
        this = tal[ka]['asis']
        it_pair = itertools.izip(this.seq_a, this.seq_b)
        
        this_seq = []
        while True:
            try:
                p = it_pair.next()
            except StopIteration:
                break
            if p is None:
                break
            if p[1] == '-':
                assert p[0] != '-', 'gap-gap?'
                this_seq.append(p[0])
            elif p[0] != '-':
                this_seq.append(p[1])
        ws = ''.join(this_seq)
        
        return ws
    

    def mytranslate(self):
        '''
        '''
        import re
        import hashlib
        from Bio.Seq import translate
        
        prot_dict = {}
        minstop = 10000
        for fr in range(3):
            s1 = self.dna_seqs[1].seq[fr:].translate()
            if minstop >=  s1.count('*'):
                self.frame = fr
                minstop = s1.count('*')
                
        print('Frame is', self.frame+1, file=sys.stderr)
        
        for s in self.dna_seqs:
            cod = Seq(s.seq[self.frame:].tostring())
            aas = translate(cod).tostring()
            ave_reads = re.search('ave_reads=(.*)', s.description).group(1)
            prot_dict[aas] = prot_dict.get(aas, 0) + float(ave_reads)
            
        for k, v in prot_dict.items():
            ts = Seq(k, IUPAC.protein)
            tsr = SeqRecord(ts, id = hashlib.sha224(k).hexdigest(), \
                            name='Reconstructed local hap: translated')
            tsr.description = 'ave_reads=%f' % v
            self.prot_seqs.append(tsr)
            
        return self.prot_seqs

    
    def mutations(self, wh='DNA', out_format='csv', out_file=sys.stdout):
        '''
        '''
        from operator import itemgetter
        from pythonlib.Alignment import dna_code
        
        if wh == 'DNA': seqs = self.dna_seqs
        elif wh == 'aa': seqs = self.prot_seqs
        
        print('#' * 60, file=sys.stderr)
        print(str(' Now %s variants ' % wh).center(60, '#'), file=sys.stderr)
        print('#' * 60, file=sys.stderr)
        dna_offset = self.offset
        aa_offset = (self.offset+1)/3
        if wh == 'DNA': print('DNA offset is', dna_offset, file=sys.stderr)
        elif wh == 'aa': print('aa offset is', aa_offset, file=sys.stderr)
        all_mut = []
        mut_info = []
        # ref_cons_mut = []
        trans_ref = self.ref[self.frame:].translate()#[aa_offset:]
        print('Translated reference is:', file=sys.stderr)
        print(trans_ref, '\n', file=sys.stderr)
        # parse the mutation in the haplotypes
        topr = True
        for s in seqs:
            mut =[]
            if wh == 'DNA':
                for i, p in enumerate(zip(self.ref, s)):
                    if dna_code[p[0].upper()] & dna_code[p[1].upper()] == set([]):
                        this_mut = '%s%d%s' % (p[0].upper(), i+1+dna_offset, p[1].upper())
                        if this_mut not in all_mut:
                            all_mut.append(this_mut)
                        mut.append(this_mut)
                        
            if wh == 'aa':
                if topr:
                    print('This is to chech that reference matches query', file=sys.stderr)
                    print(trans_ref[:10], file=sys.stderr)
                    print(s[:10].seq, file=sys.stderr)
                    print('', file=sys.stderr)
                    topr = False
                for i, p in enumerate(zip(trans_ref, s)):
                    if p[0].upper() != p[1].upper():
                        this_mut = '%s%d%s' % (p[0].upper(), i+1+aa_offset, p[1].upper())
                        if this_mut not in all_mut:
                            all_mut.append(this_mut)
                        mut.append(this_mut)
            mut_info.append([s.name, float(s.description.split('=')[1]), mut])
            
        mut_info = sorted(mut_info, key=itemgetter(1), reverse=True)
        all_mut = sorted(all_mut, key=itemgetter(slice(1,-1)), cmp=str_num_compare)
        tr_reads = sum([m[1] for m in mut_info])
        print('\nHaps account for %8.1f reads' % tr_reads, file=sys.stderr)
        print('That is %f%% of the total\n' % (100 * tr_reads/self.n_reads), file=sys.stderr)

        if out_format == 'human':
            for m in mut_info:
                print(m[0], 'ave_reads = %6.2f%%' % (100 * m[1]/tr_reads), file=sys.stderr)
                for mm in m[2]:
                    print(mm, file=sys.stderr)
                print('\n', file=sys.stderr)

        elif out_format == 'csv':
            import csv
            fh = open(out_file, 'wb')
            #writer = csv.writer(fh, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer = csv.writer(fh, dialect='excel')#delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['freq\\muts'] + all_mut)
            for m in mut_info:
                to_write =[('%4.2f' % (100 * m[1]/tr_reads))]
                for mm in all_mut:
                    if mm in m[2]:
                        to_write.append('X')
                    else:
                        to_write.append(' ')
                writer.writerow(to_write)
            
        elif out_format == 'tex':
            out_file.write('\\begin{sidewaystable}\n')
            out_file.write('\\centering\n')
            out_file.write('\\rowcolors{1}{Apricot}{cyan}\n')
            out_file.write('\\begin{tabular}{c*{%d}{|c}}\n' % len(all_mut))
            out_file.write('\\%& ')
            out_file.write(' &'.join(all_mut))
            out_file.write('\\\\\n')
            for m in mut_info:
                to_write =[('%4.2f' % (100 * m[1]/tr_reads))]
                for mm in all_mut:
                    if mm in m[2]:
                        to_write.append('X')
                    else:
                        to_write.append(' ')
                out_file.write(' &'.join(to_write))
                out_file.write('\\\\\n')
            out_file.write('\\end{tabular}\n')
            out_file.write(
            '\\caption{Patient PR: haplotypes account for %f \\%% of the reads.}\n'
            %  (100 * tr_reads/self.n_reads))
            out_file.write('\\end{sidewaystable}\n')
            
        return tr_reads, 100 * tr_reads/self.n_reads
    
    
    def get_offset(self, ref_file='~/References/HIV-HXB2.fasta', gene='protease'):
        '''
        '''
        from pythonlib import Alignment
        import os
        
        outfile = 'ppp.tmp'
        start, stop = gene_coord[gene]
        usa_seq = ref_file + '[%d:%d]' % (start, stop)
        Alignment.needle_align(usa_seq, 'asis:%s' % self.cons, outfile, go=10.0, ge=0.5)
        tal = Alignment.alignfile2dict([outfile], 'get_offset', 10.0, 0.5)
        os.remove(outfile)
        ka = tal.keys()[0]
        this = tal[ka]['asis']
        this.summary()
        self.offset = this.start
        print('Offset consensus w.r.t', ref_file, 'is', self.offset, file=sys.stderr)
        return

    
if __name__ == '__main__':
    
    try:
        sup_file = sys.argv[1]
    except ValueError:
        sys.exit('usage: %s sup_file gene[default=protease]' % sys.argv[0])
        
    gene_name = 'protease'
    
    try:
        gene_name = sys.argv[2]
    except:
        print('Allowed gene names', gene_coord.keys(), file=sys.stderr)
        print('Now using default ', gene_name, file=sys.stderr)
    else:
        if gene_name not in gene_coord.keys():
            print('Allowed gene names2', gene_coord.keys(), file=sys.stderr)
    
    sample_ls = LocalStructure(sup_file=sup_file, gene=gene_name)#, ref_file=ref_file)
    
    vd = sample_ls.alignedvariants(threshold=0.9)
    print('There are ', len(vd), ' DNA variants', file=sys.stderr)
    #SeqIO.write(vd, sys.stdout, 'fasta')
    
    pvd = sample_ls.mytranslate()
    print('There are ', len(pvd), ' protein variants', file=sys.stderr)
    
    sample_ls.mutations(wh='DNA', out_format='csv', out_file='mutations_DNA.csv')
    sample_ls.mutations(wh='aa', out_format='csv', out_file='mutations_aa.csv')
    #sample_ls.mutations(wh='DNA', out_format='human')
    
'''
n = 1
while len(qualities) != len(nucleotides):
    nucleotides, nsub = re.subn(
        '[+-]%i[ATCGNatcgn]{%i}' % (n, n), '', nucleotides)
    n += 1
    if n > 100 or len(nucleotides) < len(qualities):
        sys.exit("Something is wrong: Less nucleotides than \
        qualities and/or insertions longer than 100nt!") #TRAP!
'''
