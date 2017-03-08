# Modified by Osvaldo Zagordi (2009)

# Copyright 2006-2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# This module is for reading and writing QUAL (454) format files as SeqRecord
# objects.  The code is partly inspired  by earlier Biopython modules,
# Bio.Fasta.* and the now deprecated Bio.SeqIO.FASTA

"""Bio.SeqIO support for the "fasta" (aka FastA or Pearson) file format.

You are expected to use this module via the Bio.SeqIO functions."""

from Bio.Alphabet import generic_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter

#This is a generator function!
def QualIterator(handle, alphabet = generic_alphabet, title2ids = None) :
    """Generator function to iterate over Fasta records (as SeqRecord objects).

    handle - input file
    alphabet - optional alphabet
    title2ids - A function that, when given the title of the FASTA
    file (without the beginning >), will return the id, name and
    description (in that order) for the record as a tuple of strings.

    If this is not given, then the entire title line will be used
    as the description, and the first word as the id and name.

    Note that use of title2ids matches that of Bio.Fasta.SequenceParser
    but the defaults are slightly different.
    """
    #Skip any text before the first record (e.g. blank lines, comments)
    while True :
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line[0] == ">" :
            break

    while True :
        if line[0]!=">" :
            raise ValueError("Records in Fasta files should start with '>' character")
        if title2ids :
            id, name, descr = title2ids(line[1:].rstrip())
        else :
            descr = line[1:].rstrip()
            id   = descr.split()[0]
            name = id

        lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">": break
            #Remove trailing whitespace, and any internal spaces
            #(and any embedded \r which are possible in mangled files
            #when not opened in universal read lines mode)
            lines.append(line.rstrip().replace("\r",""))
            line = handle.readline()

        #Return the record and then continue...
        yield SeqRecord(Seq(" ".join(lines), alphabet),
                         id = id, name = name, description = descr)

        if not line : return #StopIteration
    assert False, "Should not reach this line"


if __name__ == "__main__" :
    print "Running quick self test"

    import os
    from Bio.Alphabet import generic_protein, generic_nucleotide

    #Download the files from here:
    #ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/Nanoarchaeum_equitans
    fna_filename = "NC_005213.fna"
    faa_filename = "NC_005213.faa"

    def genbank_name_function(text) :
        text, descr = text.split(None,1)
        id = text.split("|")[3]
        name = id.split(".",1)[0]
        return id, name, descr

    def print_record(record) :
        #See also bug 2057
        #http://bugzilla.open-bio.org/show_bug.cgi?id=2057
        print "ID:" + record.id
        print "Name:" + record.name
        print "Descr:" + record.description
        print record.seq
        for feature in record.annotations :
            print '/%s=%s' % (feature, record.annotations[feature])
        if record.dbxrefs :
            print "Database cross references:"
            for x in record.dbxrefs : print " - %s" % x

    if os.path.isfile(fna_filename) :
        print "--------"
        print "FastaIterator (single sequence)"
        iterator = FastaIterator(open(fna_filename, "r"), alphabet=generic_nucleotide, title2ids=genbank_name_function)
        count=0
        for record in iterator :
            count=count+1
            print_record(record)
        assert count == 1
        print str(record.__class__)

    if os.path.isfile(faa_filename) :
        print "--------"
        print "FastaIterator (multiple sequences)"
        iterator = FastaIterator(open(faa_filename, "r"), alphabet=generic_protein, title2ids=genbank_name_function)
        count=0
        for record in iterator :
            count=count+1
            print_record(record)
            break
        assert count>0
        print str(record.__class__)

    from cStringIO import StringIO
    print "--------"
    print "FastaIterator (empty input file)"
    #Just to make sure no errors happen
    iterator = FastaIterator(StringIO(""))
    count = 0
    for record in iterator :
        count = count+1
    assert count==0

    print "Done"
