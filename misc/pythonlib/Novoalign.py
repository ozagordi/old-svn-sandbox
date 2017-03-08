__author__ = "Osvaldo Zagordi"
__version__ = "$Revision: 0.1 $"
__date__ = "$Date: 2009/08/23$"
__copyright__ = ""
__license__ = ""
# Copyright 2009 by Osvaldo Zagordi
# Modification of MuscleCommandline: class implemented in Biopython
# by Cymon J. Cox. 2009

"""Command line wrapper for the short read aligner Novoalign by Novocraft (www.novocraft.com)

Citations:

Last checked against version: 2.05.04
"""
import types
import sys
from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline
DEBUG = False

class NovoalignCommandline(AbstractCommandline):
    """Command line wrapper for the short read alignment program novoalign by novocraft."""
    def __init__(self, cmd="novoalign", **kwargs):
        
        READ_FORMAT = ['FA', 'SLXFQ', 'STDFQ', 'ILMFQ', 'PRB', 'PRBnSEQ']
        REPORT_FORMAT = ['Native', 'Pairwise', 'SAM']
        REPEAT_METHOD = ['None', 'Random', 'All', 'Exhaustive', '0.99']
        
        self.parameters = \
           [
            _Option(["-d", "database"], ["input", "file"],
                    None, 0, "database filename",
                    0),
            _Option(["-f", "readfile"], ["input", "file"],
                    None, 0, "read file",
                    0),
            _Option(["-F", "format"], ["input", "option"],
                    lambda x: x in READ_FORMAT,
                    0, "Format of read files (checked anyway)",
                    0),
            
            # Alignment scoring options
            _Option(["-t", "threshold"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Threshold for alignment score",
                    0),
            _Option(["-g", "gap_open"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Gap opening penalty [default: 40]",
                    0),
            _Option(["-x", "gap_extend"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Gap extend penalty [default: 15]",
                    0),
            _Option(["-u", "unconverted"], ["input"],
                    lambda x: isinstance(x, types.IntType), 0,
                    "Experimental: unconverted cytosines penalty in bisulfite mode [default: no penalty]",
                    0),
            
            # Quality control and read filtering
            _Option(["-l", "good_bases"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Minimum number of good quality bases [default: log(N_g, 4) + 5]",
                    0),
            _Option(["-h", "homopolymer"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Homopolymer read filter [default: 20; disable: negative value]",
                    0),
            
            # Read preprocessing options
            _Option(["-a", "adapter"], ["input"],
                    lambda x: isinstance(x, types.StringType),
                    0, "Strips a 3' adapter sequence prior to alignment. With paired ends two adapters can be specified",
                    0),
            _Option(["-n", "truncate"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Truncate to specific length before alignment",
                    0),
            _Option(["-s", "trimming"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "If fail to align, trim by s bases until they map or become shorter than l [default: 2]",
                    0),
            _Option(["-5", "adapter_5"], ["input"],
                    lambda x: isinstance(x, types.StringType),
                    0, "Strips a 5' adapter sequence. Similar to -a, but on the 5'",
                    0),
            # Reporting options
            _Option(["-o", "report"], ["input"],
                    lambda x: x in REPORT_FORMAT,
                    0, "Specifies the report format [default: Native]",
                    0),
            _Option(["-Q", "quality"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Lower threshold for an alignment to be reported [default: 0]",
                    0),
            _Option(["-R", "repeats"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "If score difference is higher, report repeats, o.w. '-r method' applies [default: 5]",
                    0),
            _Option(["-r", "r_method"], ["input"],
                    lambda x: x.split()[0] in REPEAT_METHOD,
                    0, "Methods to report reads with multiple matches. 'All' and 'Exhaustive' accept limits",
                    0),
            _Option(["-e", "recorded"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Alignments recorded with score equal to the best [default: 1000 in default r_method, o.w. no limit]",
                    0),
            _Option(["-q", "qual_decimal"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Decimal digits for quality scores [default: 0]",
                    0),

            # Paired end options
            _Option(["-i", "fragment"], ["input"],
                    lambda x: len(x.split()) == 2,
                    0, "Fragment length (2 reads + insert) and standard deviation [default: 250 30]",
                    0),
            _Option(["-v", "variation"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Structural variation penalty [default: 70]",
                    0),
            
            # miRNA mode
            _Option(["-m", "miRNA"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Sets miRNA mode and optionally sets a value for the region scanned [default: off]",
                    0),
            
            # Multithreading
            _Option(["-c", "cores"], ["input"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Number of threads, disabled on free versions [default: number of cores]",
                    0),
            
            # Quality calibrations
            _Option(["-k", "read_cal"], ["input"],
                    lambda x: isinstance(x, types.StringType),
                    0, "Read quality calibration from file (mismatch counts)",
                    0),
            _Option(["-K", "write_cal"], ["input"],
                    lambda x: isinstance(x, types.StringType),
                    0, "Accumulate mismatch counts and write to file",
                    0)
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class NovoindexCommandline(AbstractCommandline):
    """Command line wrapper for the short read alignment program novoalign by novocraft."""
    def __init__(self, cmd="novoindex", **kwargs):
        
        self.parameters = \
            [
            _Option(["-k", "kmer"], ["input", "option"],
                    lambda x: isinstance(x, types.IntType),
                    0, "k-mer length used for the index, typically 14 [default: set by the program]",
                    0),
            _Option(["-s", "step"], ["input", "option"],
                    lambda x: isinstance(x, types.IntType),
                    0, "Step size used for the index, typically from 1 to 3 [default: set by the program]",
                    0),
            _Option(["-n", "name"], ["input", "option"],
                    lambda x: isinstance(x, types.StringType),
                    0, "Internal name for the reference sequence [default: indexfile nam]",
                    0),
            _Argument(["-i", "indexfile"], ["input", "option"],
                      lambda x: isinstance(x, types.StringType),
                      0, "Indexed reference sequence generated by novoindex"),
            _Argument(["-f", "sequencefiles"], ["input", "option"],
                      lambda x: isinstance(x, types.StringType),
                      0, "List of sequence files to include in the index"),
            _Switch(["-m", "masking"], ["input"],
                    "Lower case masking, if included lower case sequences are not indexed"),
            _Switch(["-b", "bisulphite"], ["input"],
                    "Turns on bisulphite mode, creating index based on C->T and G->A conversion")
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)



#from Bio.Align.Generic import Alignment
from Bio.AlignIO.Interfaces import AlignmentIterator
#from Interfaces import AlignmentIterator
#from Bio.Alphabet import single_letter_alphabet, generic_dna, generic_protein
#from Bio.Alphabet import Gapped

class NovoalignNativeIterator(AlignmentIterator) :
    """
    Alignment iterator for the novoalign (novocraft) tool's pairwise native alignment output.
    """

    def next(self) :
        '''
        Reads from the handle to construct and return the next alignment.
        '''
        handle = self.handle
        
        try :
            #Header we saved from when we were parsing
            #the previous alignment.
            line = self._header
            del self._header
        except AttributeError:      
            line = handle.readline()
        if not line:
            return None
        
        
        line = line.strip()
        if line.startswith("# novoalign") :
            #Skip the file header before the alignments
            line = self._parse_file_header(line)
            
        if line.startswith('#') and 'Read' in line:
            #Parse the footer after the alignments
            line = self._parse_file_footer(line)
        
        if line.startswith('>'):
            #Parse the footer after the alignments
            al = self._parse_alignment_line(line)
            # line, al = self._parse_alignment_line(line)
            return al

    def _parse_alignment_line(self, line):
        al_mode = {'S': 'single', 'L': 'first_file', 'R': 'second_file'}

        #Parse the alignment lines
        lsp = line.strip().split()
        
        self._id = lsp[0][1:]
        self._al_mode = al_mode[lsp[1]]
        self._query = lsp[2]
        
        if lsp[3] == 'NM':
            self._status = 'no_match'
        
        if lsp[3] == 'U':
            self._status = 'unique'
            self._score, self._quality = map(int, [lsp[4], lsp[5]])
            self._match = lsp[6][1:]
                  
        if lsp[3] == 'QC':
            self._status = 'low_quality'
            
        al = self._id

        return self

    def _parse_file_header(self, line) :
        """Helper function for the main parsing code.

        parse file header region.
        """
        #e.g. This region:
        """
        # novoalign (2.05.04) - short read aligner with qualities.
        # (C) 2008 NovoCraft
        # Licensed for evaluation, educational, and not-for-profit use only.
        #  novoalign -d ../../data/references/martin_clones_nidx/07-54825 -f /tmp/tmpwxBNCD 
        # Interpreting input files as FASTA.
        # Index Build Version: 2.5
        # Hash length: 11
        # Step size: 1
        """
        
        self._file_header_annotation = {}
#        self._file_descr = ""
        
        # line = self.handle.readline()

        assert line.startswith('#'), line+'ERROR 111'
        
        while line.startswith('#'):
            if not line :
                raise ValueError("Premature end of file!")
            if line.split()[1].startswith('novoalign'):
                self._file_header_annotation['com_line'] = ' '.join(line.split()[1:])
            if line.split()[1].startswith('Hash'):
                self._file_header_annotation['hash'] = int(line.split()[3])
            if line.split()[1].startswith('Step'):
                self._file_header_annotation['step'] = int(line.split()[3])
                
            line = self.handle.readline()
        if DEBUG:
            for k in self._file_header_annotation:
                print k+':', self._file_header_annotation[k]
            
        return line



    def _parse_file_footer(self, line):
        '''
        parse this section
        
        #     Read Sequences:   104715
        #            Aligned:    65408
        #   Unique Alignment:    65408
        #   Gapped Alignment:     1191
        #     Quality Filter:      233
        # Homopolymer Filter:        5
        #       Elapsed Time: 18,784s
        # Done.
        '''
        self._file_footer_annotation = {}

        if not line.startswith('#'):
            raise ValueError("Expected line starting '# '")
        while line:
            if line.split()[1].startswith('Read'):
                self._file_footer_annotation['read_sequences'] = int(line.split()[3])
            if line.split()[1].startswith('Aligned'):
                self._file_footer_annotation['aligned'] = int(line.split()[2])
            if line.split()[1].startswith('Unique'):
                self._file_footer_annotation['unique'] = int(line.split()[3])
            if line.split()[1].startswith('Gapped'):
                self._file_footer_annotation['gapped'] = int(line.split()[3])
            if line.split()[1].startswith('Quality'):
                self._file_footer_annotation['quality_filter'] = int(line.split()[3])
            if line.split()[1].startswith('Homopolymer'):
                self._file_footer_annotation['homopolymer_filter'] = int(line.split()[3])
            if line.split()[1].startswith('Quality'):
                self._file_footer_annotation['quality_filter'] = int(line.split()[3])
            if line.split()[1].startswith('Elapsed'):
                t = float(line.split()[3][:-1].replace(',', '.'))
                self._file_footer_annotation['elapsed_time'] = t

            line = self.handle.readline()

        if DEBUG:
            for k in self._file_footer_annotation:
                print k+':', self._file_footer_annotation[k]

        return line


if __name__ == '__main__':
    cml = NovoalignCommandline(database='~/some_dir/some_db',
                               readfile='~/some_dir/some_seq.txt')
    cml.format = 'PRBnSEQ'
    cml.r_method='0.99'
    cml.fragment = '250 20' # must be given as a string
    cml.miRNA = 100
#    print cml
#    subprocess.call(str(cml), shell=True)
    cml = NovoindexCommandline()
    cml.masking = True
    cml.indexfile='~/some_dir/some_index'
    cml.sequencefiles='~/some_dir/some_seq'
#    print cml
    h = open('./small.noa')
    aliter = NovoalignNativeIterator(h)
    alignments = list(aliter)
    print dir(alignments[1])
    print alignments[1]._file_footer_annotation is alignments[2]._file_footer_annotation

    print aliter._file_header_annotation
    print aliter._file_footer_annotation
