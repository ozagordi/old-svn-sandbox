# Copyright 2008 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
This module contains a parser for the pairwise alignments produced by Bill
Pearson's FASTA tool, for use from the Bio.AlignIO interface where it is
refered to as the "fasta-m10" file format (as we only support the machine
readable output format selected with the -m 10 command line option).

This module does NOT cover the generic "fasta" file format originally
developed as an input format to the FASTA tools.  The Bio.AlignIO and
Bio.SeqIO both use the Bio.SeqIO.FastaIO module to deal with these files,
which can also be used to store a multiple sequence alignments.
"""

from Bio.Align.Generic import Alignment
from Interfaces import AlignmentIterator

class FastaM10Iterator(AlignmentIterator) :
    """Alignment iterator for the FASTA tool's pairwise alignment output.

    This is for reading the pairwise alignments output by Bill Pearson's
    FASTA program when called with the -m 10 command line option for machine
    readable output.  For more details about the FASTA tools, see the website
    http://fasta.bioch.virginia.edu/ and the paper:

         W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

    This class is intended to be used via the Bio.AlignIO.parse() function
    by specifying the format as "fasta-m10" as shown in the following code:

        from Bio import AlignIO
        handle = ...
        for a in AlignIO.parse(handle, "fasta-m10") :
            assert len(a.get_all_seqs()) == 2, "Should be pairwise!"
            print "Alignment length %i" % a.get_alignment_length()
            for record in a :
                print record.seq, record.name, record.id

    Note that this is not a full blown parser for all the information
    in the FASTA output - for example, most of the header and all of the
    footer is ignored.  Also, the alignments are not batched according to
    the input queries.

    Also note that there can be up to about 30 letters of flanking region
    included in the raw FASTA output as contextual information.  This is NOT
    part of the alignment itself, and is not included in the resulting
    Alignment objects returned.
    """
    
    def next(self) :
        """Reads from the handle to construct and return the next alignment.

        This returns the pairwise alignment of query and match/library
        sequences as an Alignment object containing two rows."""
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

        if line.startswith("#") :
            #Skip the file header before the alignments.  e.g.
            line = self._skip_file_header(line)
        if ">>>" in line and not line.startswith(">>>") :
            #Moved onto the next query sequence!
            self._query_descr = ""
            self._query_header_annotation = {}
            #Read in the query header
            line = self._parse_query_header(line)
        if not line :
            #End of file
            return None
        if ">>><<<" in line :
            #Reached the end of the alignments, no need to read the footer...
            return None

            
        assert line.startswith(">>") and not line.startswith(">>>"), line

        query_seq_parts, match_seq_parts = [], []
        query_annotation, match_annotation = {}, {}
        match_descr = ""
        alignment_annotation = {}

        #This should be followed by the target match ID line, then more tags.
        #e.g.
        """
        >>gi|152973545|ref|YP_001338596.1| putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
        ; fa_frame: f
        ; fa_initn:  52
        ; fa_init1:  52
        ; fa_opt:  70
        ; fa_z-score: 105.5
        ; fa_bits: 27.5
        ; fa_expect:  0.082
        ; sw_score: 70
        ; sw_ident: 0.279
        ; sw_sim: 0.651
        ; sw_overlap: 43
        """
        if not line.startswith(">>") and not line.startswith(">>>") :
            raise ValueError("Expected target line starting '>>'")
        match_descr = line[2:].strip()
        #Handle the following "alignment hit" tagged data, e.g.
        line = handle.readline()
        line = self._parse_tag_section(line, alignment_annotation)
        assert not line.startswith("; ")

        #Then we have the alignment numbers and sequence for the query
        """
        >gi|10955265| ..
        ; sq_len: 346
        ; sq_offset: 1
        ; sq_type: p
        ; al_start: 197
        ; al_stop: 238
        ; al_display_start: 167
        DFMCSILNMKEIVEQKNKEFNVDIKKETIESELHSKLPKSIDKIHEDIKK
        QLSC-SLIMKKIDVEMEDYSTYCFSALRAIEGFIYQILNDVCNPSSSKNL
        GEYFTENKPKYIIREIHQET
        """
        if not (line.startswith(">") and line.strip().endswith("..")):
            raise ValueError("Expected line starting '>' and ending '..'")
        assert self._query_descr.startswith(line[1:].split()[0])
        
        #Handle the following "query alignment" tagged data
        line = handle.readline()
        line = self._parse_tag_section(line, query_annotation)
        assert not line.startswith("; ")

        #Now should have the aligned query sequence (with leading flanking region)
        while not line.startswith(">") :
            query_seq_parts.append(line.strip())
            line = handle.readline()
        
        #Handle the following "match alignment" data
        """
        >gi|152973545|ref|YP_001338596.1| ..
        ; sq_len: 242
        ; sq_type: p
        ; al_start: 52
        ; al_stop: 94
        ; al_display_start: 22
        IMTVEEARQRGARLPSMPHVRTFLRLLTGCSRINSDVARRIPGIHRDPKD
        RLSSLKQVEEALDMLISSHGEYCPLPLTMDVQAENFPEVLHTRTVRRLKR
        QDFAFTRKMRREARQVEQSW
        """
        #Match identifier
        if not (line.startswith(">") and line.strip().endswith("..")):
            raise ValueError("Expected line starting '>' and ending '..', got '%s'" % repr(line))
        assert match_descr.startswith(line[1:].split()[0])
        
        #Tagged data,
        line = handle.readline()
        line = self._parse_tag_section(line, match_annotation)
        assert not line.startswith("; ")

        #Now should have the aligned query sequence with flanking region...
        while not (line.startswith(">") or ">>>" in line):
            match_seq_parts.append(line.strip())
            line = handle.readline()
        self._header = line

        #We built a list of strings and then joined them because
        #its faster than appending to a string.
        query_seq = "".join(query_seq_parts)
        match_seq = "".join(match_seq_parts)
        del query_seq_parts, match_seq_parts
        #Note, query_seq and match_seq will usually be of different lengths, apparently
        #because in the m10 format leading gaps are added but not trailing gaps!

        #Remove the flanking regions,
        query_align_seq = self._extract_alignment_region(query_seq, query_annotation)
        match_align_seq = self._extract_alignment_region(match_seq, match_annotation)

        #The "sq_offset" values can be specified with the -X command line option.
        #The appear to just shift the origin used in the calculation of the coordinates.
        
        if ("sq_offset" in query_annotation and query_annotation["sq_offset"]<>"1") \
        or ("sq_offset" in match_annotation and match_annotation["sq_offset"]<>"1") :
            #Note that until some point in the v35 series, FASTA always recorded one
            #for the query offset, and ommitted the match offset (even when these were
            #query_seq the -X command line option).
            #TODO - Work out how exactly the use of -X offsets changes things.
            #raise ValueError("Offsets from the -X command line option are not (yet) supported")
            pass

        if len(query_align_seq) <> len(match_align_seq) :
            raise ValueError("Problem parsing the alignment sequence coordinates")
        if "sw_overlap" in alignment_annotation :
            if int(alignment_annotation["sw_overlap"]) <> len(query_align_seq) :
                raise ValueError("Specified sw_overlap = %s does not match expected value %i" \
                                 % (alignment_annotation["sw_overlap"],
                                    len(query_align_seq)))

        #TODO - Look at the "sq_type" to assign a sensible alphabet?
        alignment = Alignment(self.alphabet)

        #TODO - Introduce an annotated alignment class?
        #For now, store the annotation a new private property:
        alignment._annotations = {}
        
        #Want to record both the query header tags, and the alignment tags.
        for key, value in self._query_header_annotation.iteritems() :
            alignment._annotations[key] = value
        for key, value in alignment_annotation.iteritems() :
            alignment._annotations[key] = value
            

        #TODO - Once the alignment object gets an append method, use it.
        #(i.e. an add SeqRecord method)
        alignment.add_sequence(self._query_descr, query_align_seq)
        record = alignment.get_all_seqs()[-1]
        assert record.id == self._query_descr or record.description == self._query_descr
        assert record.seq.tostring() == query_align_seq
        record.id = self._query_descr.split()[0].strip(",")
        record.name = "query"
        record.annotations["original_length"] = int(query_annotation["sq_len"])
        
        alignment.add_sequence(match_descr, match_align_seq)
        record = alignment.get_all_seqs()[-1]
        assert record.id == match_descr or record.description == match_descr
        assert record.seq.tostring() == match_align_seq
        record.id = match_descr.split()[0].strip(",")
        record.name = "match"
        record.annotations["original_length"] = int(match_annotation["sq_len"])

        return alignment

    def _skip_file_header(self, line) :
        """Helper function for the main parsing code.

        Skips over the file header region.
        """
        #e.g. This region:
        """
        # /home/xxx/Downloads/FASTA/fasta-35.3.6/fasta35 -Q -H -E 1 -m 10 -X "-5 -5" NC_002127.faa NC_009649.faa
        FASTA searches a protein or DNA sequence data bank
         version 35.03 Feb. 18, 2008
        Please cite:
         W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

        Query: NC_002127.faa
        """
        #Note that there is no point recording the command line here
        #from the # line, as it is included again in each alignment
        #under the "pg_argv" tag.  Likewise for the program version.
        while ">>>" not in line :
            line = self.handle.readline()
        return line

    def _parse_query_header(self, line) :
        """Helper function for the main parsing code.

        Skips over the free format query header, extracting the tagged parameters.

        If there are no hits for the current query, it is skipped entirely."""
        #e.g. this region (where there is often a histogram too):
        """
          2>>>gi|10955264|ref|NP_052605.1| hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa - 126 aa
        Library: NC_009649.faa   45119 residues in   180 sequences

          45119 residues in   180 sequences
        Statistics: (shuffled [500]) Expectation_n fit: rho(ln(x))= 5.0398+/-0.00968; mu= 2.8364+/- 0.508
         mean_var=44.7937+/-10.479, 0's: 0 Z-trim: 0  B-trim: 0 in 0/32
         Lambda= 0.191631
        Algorithm: FASTA (3.5 Sept 2006) [optimized]
        Parameters: BL50 matrix (15:-5) ktup: 2
         join: 36, opt: 24, open/ext: -10/-2, width:  16
         Scan time:  0.040

        The best scores are:                                      opt bits E(180)
        gi|152973462|ref|YP_001338513.1| hypothetical prot ( 101)   58 23.3    0.22
        gi|152973501|ref|YP_001338552.1| hypothetical prot ( 245)   55 22.5    0.93
        """
        #Sometime have queries with no matches, in which case we continue to the next query block:
        """
          2>>>gi|152973838|ref|YP_001338875.1| hypothetical protein KPN_pKPN7p10263 [Klebsiella pneumoniae subsp. pneumonia 76 aa - 76 aa
         vs  NC_002127.faa library

            579 residues in     3 sequences
         Altschul/Gish params: n0: 76 Lambda: 0.158 K: 0.019 H: 0.100

        FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
         join: 36, opt: 24, open/ext: -10/-2, width:  16
         Scan time:  0.000
        !! No library sequences with E() < 0.5
        """

        self._query_header_annotation = {}
        self._query_descr = ""

        assert ">>>" in line and not line.startswith(">>>")
        #There is nothing useful in this line, the query description is truncated.
        
        line = self.handle.readline()
        #We ignore the free form text...
        while not line.startswith(">>>") :
            #print "Ignoring %s" % line.strip()
            line = self.handle.readline()
            if not line :
                raise ValueError("Premature end of file!")
            if ">>><<<" in line :
                #End of alignments, looks like the last query
                #or queries had no hits.
                return line

        #Now want to parse this section:
        """
        >>>gi|10955264|ref|NP_052605.1|, 126 aa vs NC_009649.faa library
        ; pg_name: /home/pjcock/Downloads/FASTA/fasta-35.3.6/fasta35
        ; pg_ver: 35.03
        ; pg_argv: /home/pjcock/Downloads/FASTA/fasta-35.3.6/fasta35 -Q -H -E 1 -m 10 -X -5 -5 NC_002127.faa NC_009649.faa
        ; pg_name: FASTA
        ; pg_ver: 3.5 Sept 2006
        ; pg_matrix: BL50 (15:-5)
        ; pg_open-ext: -10 -2
        ; pg_ktup: 2
        ; pg_optcut: 24
        ; pg_cgap: 36
        ; mp_extrap: 60000 500
        ; mp_stats: (shuffled [500]) Expectation_n fit: rho(ln(x))= 5.0398+/-0.00968; mu= 2.8364+/- 0.508  mean_var=44.7937+/-10.479, 0's: 0 Z-trim: 0  B-trim: 0 in 0/32  Lambda= 0.191631
        ; mp_KS: -0.0000 (N=1066338402) at  20
        ; mp_Algorithm: FASTA (3.5 Sept 2006) [optimized]
        ; mp_Parameters: BL50 matrix (15:-5) ktup: 2  join: 36, opt: 24, open/ext: -10/-2, width:  16
        """

        assert line.startswith(">>>"), line
        self._query_descr = line[3:].strip()

        #Handle the following "program" tagged data,
        line = self.handle.readline()
        line = self._parse_tag_section(line, self._query_header_annotation)
        assert not line.startswith("; "), line
        assert line.startswith(">>"), line
        return line


    def _extract_alignment_region(self, alignment_seq_with_flanking, annotation) :
        """Helper function for the main parsing code.

        To get the actual pairwise alignment sequences, we must first
        translate the un-gapped sequence based coordinates into positions
        in the gapped sequence (which may have a flanking region shown
        using leading - characters).  To date, I have never seen any
        trailing flanking region shown in the m10 file, but the
        following code should also cope with that.

        Note that this code seems to work fine even when the "sq_offset"
        entries are prsent as a result of using the -X command line option.
        """
        align_stripped = alignment_seq_with_flanking.strip("-")
        start = int(annotation['al_start']) \
              - int(annotation['al_display_start'])
        end   = int(annotation['al_stop']) \
              - int(annotation['al_display_start']) \
              + align_stripped.count("-") + 1
        return align_stripped[start:end]


    def _parse_tag_section(self, line, dictionary) :
        """Helper function for the main parsing code.

        line - supply line just read from the handle containing the start of
               the tagged section.
        dictionary - where to record the tagged values

        Returns a string, the first line following the tagged section."""
        if not line.startswith("; ") :
            raise ValueError("Expected line starting '; '")
        while line.startswith("; ") :
            tag, value = line[2:].strip().split(": ",1)
            #fasta34 and fasta35 will reuse the pg_name and pg_ver tags
            #for the program executable and name, and the program version
            #and the algorithm version, respectively.  This may be a bug.
            #if tag in dictionary :
            #    raise ValueError("Repeated tag '%s' in section" % tag)
            dictionary[tag] = value
            line = self.handle.readline()
        return line
    
if __name__ == "__main__" :
    print "Running a quick self-test"

    #http://emboss.sourceforge.net/docs/themes/alnformats/align.simple
    simple_example = \
"""# /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
FASTA searches a protein or DNA sequence data bank
 version 34.26 January 12, 2007
Please cite:
 W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

Query library NC_002127.faa vs NC_009649.faa library
searching NC_009649.faa library

  1>>>gi|10955263|ref|NP_052604.1| plasmid mobilization [Escherichia coli O157:H7 s 107 aa - 107 aa
 vs  NC_009649.faa library

  45119 residues in   180 sequences
  Expectation_n fit: rho(ln(x))= 6.9146+/-0.0249; mu= -5.7948+/- 1.273
 mean_var=53.6859+/-13.609, 0's: 0 Z-trim: 1  B-trim: 9 in 1/25
 Lambda= 0.175043

FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
 join: 36, opt: 24, open/ext: -10/-2, width:  16
 Scan time:  0.000
The best scores are:                                      opt bits E(180)
gi|152973457|ref|YP_001338508.1| ATPase with chape ( 931)   71 24.9    0.58
gi|152973588|ref|YP_001338639.1| F pilus assembly  ( 459)   63 23.1    0.99

>>>gi|10955263|ref|NP_052604.1|, 107 aa vs NC_009649.faa library
; pg_name: /opt/fasta/fasta34
; pg_ver: 34.26
; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
; pg_name: FASTA
; pg_ver: 3.5 Sept 2006
; pg_matrix: BL50 (15:-5)
; pg_open-ext: -10 -2
; pg_ktup: 2
; pg_optcut: 24
; pg_cgap: 36
; mp_extrap: 60000 180
; mp_stats:  Expectation_n fit: rho(ln(x))= 6.9146+/-0.0249; mu= -5.7948+/- 1.273  mean_var=53.6859+/-13.609, 0's: 0 Z-trim: 1  B-trim: 9 in 1/25  Lambda= 0.175043
; mp_KS: -0.0000 (N=0) at 8159228
>>gi|152973457|ref|YP_001338508.1| ATPase with chaperone activity, ATP-binding subunit [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  65
; fa_init1:  43
; fa_opt:  71
; fa_z-score: 90.3
; fa_bits: 24.9
; fa_expect:   0.58
; sw_score: 71
; sw_ident: 0.250
; sw_sim: 0.574
; sw_overlap: 108
>gi|10955263| ..
; sq_len: 107
; sq_offset: 1
; sq_type: p
; al_start: 5
; al_stop: 103
; al_display_start: 1
--------------------------MTKRSGSNT-RRRAISRPVRLTAE
ED---QEIRKRAAECGKTVSGFLRAAALGKKVNSLTDDRVLKEVM-----
RLGALQKKLFIDGKRVGDREYAEVLIAITEYHRALLSRLMAD
>gi|152973457|ref|YP_001338508.1| ..
; sq_len: 931
; sq_type: p
; al_start: 96
; al_stop: 195
; al_display_start: 66
SDFFRIGDDATPVAADTDDVVDASFGEPAAAGSGAPRRRGSGLASRISEQ
SEALLQEAAKHAAEFGRS------EVDTEHLLLALADSDVVKTILGQFKI
KVDDLKRQIESEAKR-GDKPF-EGEIGVSPRVKDALSRAFVASNELGHSY
VGPEHFLIGLAEEGEGLAANLLRRYGLTPQ
>>gi|152973588|ref|YP_001338639.1| F pilus assembly protein [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  33
; fa_init1:  33
; fa_opt:  63
; fa_z-score: 86.1
; fa_bits: 23.1
; fa_expect:   0.99
; sw_score: 63
; sw_ident: 0.266
; sw_sim: 0.656
; sw_overlap: 64
>gi|10955263| ..
; sq_len: 107
; sq_offset: 1
; sq_type: p
; al_start: 32
; al_stop: 94
; al_display_start: 2
TKRSGSNTRRRAISRPVRLTAEEDQEIRKRAAECGKTVSGFLRAAALGKK
VNSLTDDRVLKEV-MRLGALQKKLFIDGKRVGDREYAEVLIAITEYHRAL
LSRLMAD
>gi|152973588|ref|YP_001338639.1| ..
; sq_len: 459
; sq_type: p
; al_start: 191
; al_stop: 248
; al_display_start: 161
VGGLFPRTQVAQQKVCQDIAGESNIFSDWAASRQGCTVGG--KMDSVQDK
ASDKDKERVMKNINIMWNALSKNRLFDG----NKELKEFIMTLTGTLIFG
ENSEITPLPARTTDQDLIRAMMEGGTAKIYHCNDSDKCLKVVADATVTIT
SNKALKSQISALLSSIQNKAVADEKLTDQE
  2>>>gi|10955264|ref|NP_052605.1| hypothetical protein pOSAK1_02 [Escherichia coli O157:H7 s 126 aa - 126 aa
 vs  NC_009649.faa library

  45119 residues in   180 sequences
  Expectation_n fit: rho(ln(x))= 7.1374+/-0.0246; mu= -7.6540+/- 1.313
 mean_var=51.1189+/-13.171, 0's: 0 Z-trim: 1  B-trim: 8 in 1/25
 Lambda= 0.179384

FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
 join: 36, opt: 24, open/ext: -10/-2, width:  16
 Scan time:  0.000
The best scores are:                                      opt bits E(180)
gi|152973462|ref|YP_001338513.1| hypothetical prot ( 101)   58 22.9    0.29

>>>gi|10955264|ref|NP_052605.1|, 126 aa vs NC_009649.faa library
; pg_name: /opt/fasta/fasta34
; pg_ver: 34.26
; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
; pg_name: FASTA
; pg_ver: 3.5 Sept 2006
; pg_matrix: BL50 (15:-5)
; pg_open-ext: -10 -2
; pg_ktup: 2
; pg_optcut: 24
; pg_cgap: 36
; mp_extrap: 60000 180
; mp_stats:  Expectation_n fit: rho(ln(x))= 7.1374+/-0.0246; mu= -7.6540+/- 1.313  mean_var=51.1189+/-13.171, 0's: 0 Z-trim: 1  B-trim: 8 in 1/25  Lambda= 0.179384
; mp_KS: -0.0000 (N=0) at 8159228
>>gi|152973462|ref|YP_001338513.1| hypothetical protein KPN_pKPN3p05904 [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  50
; fa_init1:  50
; fa_opt:  58
; fa_z-score: 95.8
; fa_bits: 22.9
; fa_expect:   0.29
; sw_score: 58
; sw_ident: 0.289
; sw_sim: 0.632
; sw_overlap: 38
>gi|10955264| ..
; sq_len: 126
; sq_offset: 1
; sq_type: p
; al_start: 1
; al_stop: 38
; al_display_start: 1
------------------------------MKKDKKYQIEAIKNKDKTLF
IVYATDIYSPSEFFSKIESDLKKKKSKGDVFFDLIIPNGGKKDRYVYTSF
NGEKFSSYTLNKVTKTDEYN
>gi|152973462|ref|YP_001338513.1| ..
; sq_len: 101
; sq_type: p
; al_start: 44
; al_stop: 81
; al_display_start: 14
DALLGEIQRLRKQVHQLQLERDILTKANELIKKDLGVSFLKLKNREKTLI
VDALKKKYPVAELLSVLQLARSCYFYQNVCTISMRKYA
  3>>>gi|10955265|ref|NP_052606.1| hypothetical protein pOSAK1_03 [Escherichia coli O157:H7 s 346 aa - 346 aa
 vs  NC_009649.faa library

  45119 residues in   180 sequences
  Expectation_n fit: rho(ln(x))= 6.0276+/-0.0276; mu= 3.0670+/- 1.461
 mean_var=37.1634+/- 8.980, 0's: 0 Z-trim: 1  B-trim: 14 in 1/25
 Lambda= 0.210386

FASTA (3.5 Sept 2006) function [optimized, BL50 matrix (15:-5)] ktup: 2
 join: 37, opt: 25, open/ext: -10/-2, width:  16
 Scan time:  0.020
The best scores are:                                      opt bits E(180)
gi|152973545|ref|YP_001338596.1| putative plasmid  ( 242)   70 27.5   0.082

>>>gi|10955265|ref|NP_052606.1|, 346 aa vs NC_009649.faa library
; pg_name: /opt/fasta/fasta34
; pg_ver: 34.26
; pg_argv: /opt/fasta/fasta34 -Q -H -E 1 -m 10 NC_002127.faa NC_009649.faa
; pg_name: FASTA
; pg_ver: 3.5 Sept 2006
; pg_matrix: BL50 (15:-5)
; pg_open-ext: -10 -2
; pg_ktup: 2
; pg_optcut: 25
; pg_cgap: 37
; mp_extrap: 60000 180
; mp_stats:  Expectation_n fit: rho(ln(x))= 6.0276+/-0.0276; mu= 3.0670+/- 1.461  mean_var=37.1634+/- 8.980, 0's: 0 Z-trim: 1  B-trim: 14 in 1/25  Lambda= 0.210386
; mp_KS: -0.0000 (N=0) at 8159228
>>gi|152973545|ref|YP_001338596.1| putative plasmid SOS inhibition protein A [Klebsiella pneumoniae subsp. pneumoniae MGH 78578]
; fa_frame: f
; fa_initn:  52
; fa_init1:  52
; fa_opt:  70
; fa_z-score: 105.5
; fa_bits: 27.5
; fa_expect:  0.082
; sw_score: 70
; sw_ident: 0.279
; sw_sim: 0.651
; sw_overlap: 43
>gi|10955265| ..
; sq_len: 346
; sq_offset: 1
; sq_type: p
; al_start: 197
; al_stop: 238
; al_display_start: 167
DFMCSILNMKEIVEQKNKEFNVDIKKETIESELHSKLPKSIDKIHEDIKK
QLSC-SLIMKKIDVEMEDYSTYCFSALRAIEGFIYQILNDVCNPSSSKNL
GEYFTENKPKYIIREIHQET
>gi|152973545|ref|YP_001338596.1| ..
; sq_len: 242
; sq_type: p
; al_start: 52
; al_stop: 94
; al_display_start: 22
IMTVEEARQRGARLPSMPHVRTFLRLLTGCSRINSDVARRIPGIHRDPKD
RLSSLKQVEEALDMLISSHGEYCPLPLTMDVQAENFPEVLHTRTVRRLKR
QDFAFTRKMRREARQVEQSW
>>><<<


579 residues in 3 query   sequences
45119 residues in 180 library sequences
 Scomplib [34.26]
 start: Tue May 20 16:38:45 2008 done: Tue May 20 16:38:45 2008
 Total Scan time:  0.020 Total Display time:  0.010

Function used was FASTA [version 34.26 January 12, 2007]

"""                 


    from StringIO import StringIO

    alignments = list(FastaM10Iterator(StringIO(simple_example)))
    assert len(alignments) == 4, len(alignments)
    assert len(alignments[0].get_all_seqs()) == 2
    for a in alignments :
        print "Alignment %i sequences of length %i" \
              % (len(a.get_all_seqs()), a.get_alignment_length())
        for r in a :
            print "%s %s %i" % (r.seq, r.id, r.annotations["original_length"])
        #print a.annotations
    print "Done"

    import os
    path = "./"
    files = [f for f in os.listdir(path) if os.path.splitext(f)[-1] == ".m10"]
    files.sort()
    for filename in files :
        if os.path.splitext(filename)[-1] == ".m10" :
            print
            print filename
            print "="*len(filename)
            for a in FastaM10Iterator(open(os.path.join(path,filename))):
                print a
