import codonTable
import sys
import bx.seq.twobit
class GtfRec(object):
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attributes):
        """represent a Gene Transfer Format File - http://mblab.wustl.edu/GTF2.html"""
        #gtf is one based - but make it zero-based, half-open
        self.seqname=seqname
        self.source=source
        self.feature=feature
        self.start=int(start)-1
        self.end=int(end)
        self.score=score
        self.strand=strand
        self.frame=int(frame)
        self.attributes=attributes

    def getSeqName(self): return self.seqname
    def getSource(self): return self.source
    def getFeature(self): return self.feature
    def getStart(self): return self.start
    def getEnd(self): return self.end
    def getScore(self): return self.score
    def getStrand(self): return self.strand
    def getFrame(self): return self.frame
    def getAttributes(self): return self.attributes

    def getSequence(self, twobit):
        """ return sequence of Gtf object's feature given a bx seq twobit object"""
        sequence=''
        try:
            sequence = twobit[self.seqname][self.start:self.end]
            sequence=sequence.upper()
        except:
            sys.stderr.write("Unable to fetch the sequence\n")
        if self.strand == '-':
            return codonTable.reverse_complement(sequence)
        return sequence


    def walkCds(self, twobit):
        if self.strand == '+':
            self.walkCdsFwd(twobit)
        if self.strand == '-':
            self.walkCdsMinus(twobit)

    def walkCdsFwd(self, twobit):
        """  given a twobit object generator to walk CDS feature and yield codon, aa & codonStart, codonEnd intervals """
       
        #walk the CDS feature in frame, so adjust the counter variable i if needed
        i=self.start
        if self.frame == 1:
            i=i+1
        if self.frame == 2:
            i= i+2
        while (i +3 < self.end):
            sequence =twobit[self.seqname][i:i+3]
            sequence=sequence.upper()
            yield sequence, codonTable.translateCodon(sequence), i, i+3
            i=i+3

    def walkCdsMinus(self, twobit):
        """  given a twobit object generator to walk CDS feature and yield codon, aa & codonStart, codonEnd intervals """
        #walk the CDS feature in frame, so adjust the counter variable i if needed

        i=self.end
        
        if self.frame == 1:
            i=i-1
        if self.frame == 2:
            i=i-2
        while ( i > self.start+1):
            sequence=twobit[self.seqname][i-3:i]
            sequence=sequence.upper()
            sequence=codonTable.reverse_complement(sequence)
            yield sequence, codonTable.translateCodon(sequence), i-3, i
            i=i-3

    def getSequenceInFrame(self, twobit):
        """return a sequence of Gtf object's feature given a twobit object
        return in frame - if frame is 1,2 remove the partial codon
        0 indicates that the first whole codon of the reading frame is located at 5'-most base.
        1 means that there is one extra base before the first codon and
        2 means that there are two extra bases before the first codon """

        sequence=self.getSequence(twobit)
        if self.frame==2:
            sequence=sequence[2::]
        if self.frame==1:
            sequence=sequence[1::]

        return sequence