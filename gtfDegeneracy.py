#!/usr/bin/python
from Gtf import GtfRec
from gtfIO import *
from optparse import OptionParser
from bx.intervals import *
import bx.seq.twobit
import codonTable
import sys


def walkExonSequence(sequence):
    
    while (i < len(sequence) ):
        yield sequence[i:i+3], codonTable.translateCodon(sequence[i:i+3])
        i+=3


def main():
    """ walk the codons of CDS gtf feature and write the bed interval of nsyn, 2-4x degeneracy to appropriate bed file  """

    usage="usage: %prog [options] file.gtf\n\nWalk the codons of CDS gtf feature and write the bed interval of nsyn, 2-4x degeneracy to appropriate bed file\n"
    parser=OptionParser(usage)
    parser.add_option("--twoBitFile", type="string", dest="tbf", help="*.2bit file to  extract sequence from")
    parser.add_option("--debug", action="store_true", dest="debug", default = False, help="print debug statements")
    (options, args)=parser.parse_args()


    onefoldFH = open('gencode.oneFold.bed', 'w')
    twofoldFH = open('gencode.twoFold.bed', 'w')
    threefoldFH = open('gencode.threeFold.bed', 'w')
    fourfoldFH = open('gencode.fourFold.bed', 'w')

    if len(args) != 1:
        parser.error("incorrect number of arguments")
        parser.usage()

    gtf_file=args[0]

    try:
        gtf_fh=open(gtf_file)
    except:
        sys.stderr.write("cannot open gtf file!\n")

    try:
        twobit=bx.seq.twobit.TwoBitFile( open( options.tbf ) ) 
    except:
        sys.stderr.write("Unable to open twobit file!\n")
    sys.stderr.write("walking CDS exons in gtf ....\n")
    #iterate CDS records in gtf
    for cdsgtf in getCds(gtf_fh):

        if cdsgtf.getStrand() == '+':
            for codon, aa, codonStart, codonEnd in cdsgtf.walkCdsFwd( twobit ):
                if options.debug == True: print codon, aa, cdsgtf.getSeqName(), codonStart, codonEnd, codonTable.codon_degeneracy(codon,1), codonTable.codon_degeneracy(codon,2),codonTable.codon_degeneracy(codon,3), cdsgtf.getStrand()

                #determine degeneracy in each position of a codon
                degenThree = codonTable.codon_degeneracy(codon,3)
                degenTwo = codonTable.codon_degeneracy(codon,2)
                degenOne = codonTable.codon_degeneracy(codon,1)
                if degenThree == 1:
                    onefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "nsyn\n"]))
                if  degenThree == 2:
                    twofoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "twoFold\n"]))
                if  degenThree == 3:
                    threefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "threeFold\n"]))
                if  degenThree == 4:
                    fourfoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "fourdFold\n"]))
                
                if degenTwo == 1:
                    onefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "oneFold\n"]))
                if  degenTwo == 2:
                    twofoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "twoFold\n"]))
                if  degenTwo == 3:
                    threefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "threeFold\n"]))
                if  degenTwo == 4:
                    fourfoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "fourdFold\n"]))

                if degenOne == 1:
                    onefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "nsyn\n"]))
                if  degenOne == 2:
                    twofoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "twoFold\n"]))
                if  degenOne == 3:
                    threefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "threeFold\n"]))
                if  degenOne == 4:
                    fourfoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "fourdFold\n"]))

        if cdsgtf.getStrand() == '-':
            for codon, aa, codonStart, codonEnd in cdsgtf.walkCdsMinus( twobit ):
                if options.debug == True: print codon, aa, cdsgtf.getSeqName(), codonStart, codonEnd, codonTable.codon_degeneracy(codon,1),codonTable.codon_degeneracy(codon,2),codonTable.codon_degeneracy(codon,3), cdsgtf.getStrand()
                #determine degeneracy in each position of a codon
                degenThree = codonTable.codon_degeneracy(codon,3)
                degenTwo = codonTable.codon_degeneracy(codon,2)
                degenOne = codonTable.codon_degeneracy(codon,1)
                
                if degenThree == 1:
                    onefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "nsyn\n"]))
                if  degenThree == 2:
                    twofoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "twoFold\n"]))
                if  degenThree == 3:
                    threefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "threeFold\n"]))
                if  degenThree == 4:
                    fourfoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart), str(codonStart+1), "fourdFold\n"]))
                
                if degenTwo == 1:
                    onefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "nsyn\n"]))
                if  degenTwo == 2:
                    twofoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "twoFold\n"]))
                if  degenTwo == 3:
                    threefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "threeFold\n"]))
                if  degenTwo == 4:
                    fourfoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+1), str(codonStart+2), "fourdFold\n"]))

                if degenOne == 1:
                    onefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "nsyn\n"]))
                if  degenOne == 2:
                    twofoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "twoFold\n"]))
                if  degenOne == 3:
                    threefoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "threeFold\n"]))
                if  degenOne == 4:
                    fourfoldFH.write("\t".join([cdsgtf.getSeqName(), str(codonStart+2), str(codonStart+3), "fourdFold\n"]))

if __name__ == "__main__":
    main()
