# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Alphabet import single_letter_alphabet, Gapped
   
class AlignmentIterator :
    """Base class for building Alignment iterators.

    You should write a next() method to return Aligment
    objects.  You may wish to redefine the __init__
    method as well.
    """
    #TODO - Should the default be Gapped(single_letter_alphabet) instead?
    def __init__(self, handle, count=None,
                 alphabet = single_letter_alphabet) :
        """Create an AlignmentIterator object.

        handle   - input file
        count    - optional, expected number of records per alignment
                   Recommend for fasta file format.
        alphabet - optional, e.g. Bio.Alphabet.generic_protein

        Note when subclassing:
        - there should be a single non-optional argument, the handle,
          and optional count and alphabet IN THAT ORDER.
        - you do not have to require an alphabet (?).
        - you can add additional optional arguments."""
        self.handle = handle
        self.records_per_alignment = count
        self.alphabet = alphabet
        #####################################################
        # You may want to subclass this, for example        #
        # to read through the file to find the first record,#
        # or if additional arguments are required.          #
        #####################################################

    def next(self) :
        """Return the next alignment in the file.
        
        This method should be replaced by any derived class to do something
        useful."""
        raise NotImplementedError, "This object should be subclassed"
        #####################################################
        # You SHOULD subclass this, to split the file up    #
        # into your individual alignments and convert these #
        # into Alignment objects.                           #
        #####################################################

    def __iter__(self):
        """Iterate over the entries as Alignment objects.

        Example usage for (concatenated) PHYLIP files:

        myFile = open("many.phy","r")
        for alignment in PhylipIterator(myFile) :
            print "New alignment:"
            for record in alignment :
                print record.id
                print record.seq
        myFile.close()"""
        return iter(self.next, None)

class AlignmentWriter :
    """Base class for building Alignment writers.
    
    You should write a write_alignment() method.
    You may wish to redefine the __init__ method as well"""

    def __init__(self, handle) :
        self.handle = handle
       
    def write_file(self, alignments) :
        """Use this to write an entire file containing the given alignments.

        alignments - A list or iterator returning Alignment objects

        In general, this method can only be called once per file.
        
        This method should be replaced by any derived class to do something
        useful."""
        raise NotImplementedError, "This object should be subclassed"
        #####################################################
        # You SHOULD subclass this, to write the alignment  #
        # objecta to the file handle                        #
        #####################################################

    def clean(self, text) :
        """Use this to avoid getting newlines in the output."""
        answer = text
        for x in ["\n", "\r"] :
            answer = answer.replace(x, " ")
        return answer.replace("  ", " ")
    
class SequentialAlignmentWriter(AlignmentWriter) :
    """Base class for building Alignment writers.
    
    This assumes each alignment can be simply appended to the file.
    You should write a write_alignment() method.
    You may wish to redefine the __init__ method as well"""

    def __init__(self, handle) :
        self.handle = handle
       
    def write_file(self, alignments) :
        """Use this to write an entire file containing the given alignments.

        alignments - A list or iterator returning Alignment objects

        In general, this method can only be called once per file."""
        self.write_header()
        for alignment in alignments :
            self.write_alignment(alignment)
        self.write_footer()
        
    def write_header(self) :
        """Use this to write any header.
        
        This method should be replaced by any derived class to do something
        useful."""
        pass
    
    def write_footer(self) :
        """Use this to write any footer.
        
        This method should be replaced by any derived class to do something
        useful."""
        pass

    def write_alignment(self, alignment) :
        """Use this to write a single alignment.
        
        This method should be replaced by any derived class to do something
        useful."""
        raise NotImplementedError, "This object should be subclassed"
        #####################################################
        # You SHOULD subclass this, to write the alignment  #
        # objecta to the file handle                        #
        #####################################################
