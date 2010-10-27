import re
import string
GENETIC_CODE = """
TTT (Phe/F)Phenylalanine
TTC (Phe/F)Phenylalanine
TTA (Leu/L)Leucine
TTG (Leu/L)Leucine, Start
TCT (Ser/S)Serine
TCC (Ser/S)Serine
TCA (Ser/S)Serine
TCG (Ser/S)Serine
TAT (Tyr/Y)Tyrosine
TAC (Tyr/Y)Tyrosine
TAA Ochre (Stop)
TAG Amber (Stop)
TGT (Cys/C)Cysteine
TGC (Cys/C)Cysteine
TGA Opal (Stop)
TGG (Trp/W)Tryptophan
CTT (Leu/L)Leucine
CTC (Leu/L)Leucine
CTA (Leu/L)Leucine
CTG (Leu/L)Leucine, Start
CCT (Pro/P)Proline
CCC (Pro/P)Proline
CCA (Pro/P)Proline
CCG (Pro/P)Proline
CAT (His/H)Histidine
CAC (His/H)Histidine
CAA (Gln/Q)Glutamine
CAG (Gln/Q)Glutamine
CGT (Arg/R)Arginine
CGC (Arg/R)Arginine
CGA (Arg/R)Arginine
CGG (Arg/R)Arginine
ATT (Ile/I)Isoleucine, Start2
ATC (Ile/I)Isoleucine
ATA (Ile/I)Isoleucine
ATG (Met/M)Methionine, Start1
ACT (Thr/T)Threonine
ACC (Thr/T)Threonine
ACA (Thr/T)Threonine
ACG (Thr/T)Threonine
AAT (Asn/N)Asparagine
AAC (Asn/N)Asparagine
AAA (Lys/K)Lysine
AAG (Lys/K)Lysine
AGT (Ser/S)Serine
AGC (Ser/S)Serine
AGA (Arg/R)Arginine
AGG (Arg/R)Arginine
GTT (Val/V)Valine
GTC (Val/V)Valine
GTA (Val/V)Valine
GTG (Val/V)Valine, Start2
GCT (Ala/A)Alanine
GCC (Ala/A)Alanine
GCA (Ala/A)Alanine
GCG (Ala/A)Alanine
GAT (Asp/D)Aspartic acid
GAC (Asp/D)Aspartic acid
GAA (Glu/E)Glutamic acid
GAG (Glu/E)Glutamic acid
GGT (Gly/G)Glycine
GGC (Gly/G)Glycine
GGA (Gly/G)Glycine
GGG (Gly/G)Glycine
"""



""" parse the doc string to hash the genetic code"""
GEN_CODE = {}
for line in GENETIC_CODE.split('\n'):
    if line.strip() == '': continue
    f = re.split('\s|\(|\)|\/',line)
    codon = f[0]
    c1,c2,c3 = codon
    aminoacid = f[3]
    if c1 not in GEN_CODE: GEN_CODE[c1] = {}
    if c2 not in GEN_CODE[c1]: GEN_CODE[c1][c2] = {}

    GEN_CODE[c1][c2][c3] = aminoacid



def translate( codon, genetic_code):
    c1,c2,c3 = codon
    return genetic_code[c1][c2][c3]

def dumpGeneTable():
    print GEN_CODE

REVMAP = string.maketrans("ACGTacgt","TGCAtgca")
def revComp(seq):
    return seq[::-1].translate(REVMAP)

def Comp(seq):
    return seq.translate(REVMAP)


def codon_degeneracy( codon, position):
    all = ['A','C','G','T']
    aa = translate( codon, GEN_CODE )
    if position==1:
        degeneracy = [GEN_CODE[ k ][ codon[1] ][ codon[2] ] for k in all].count(aa)
    elif position==2:
        degeneracy = [GEN_CODE[ codon[0] ][ k ][ codon[2] ] for k in all].count(aa)
    elif position==3:
        degeneracy = GEN_CODE[ codon[0] ][ codon[1] ].values().count(aa)

    return degeneracy





def codon_syn( codon, position=3 ):
    aa = translate( codon, GEN_CODE )
    all = ['A','C','G','T']
    possible_aa = []
    if position==1:
        all = filter(lambda l: codon[0] not in l, all)
        possible_aa = [GEN_CODE[ k ][ codon[1] ][ codon[2] ] for k in all]
        possible_aa = filter(lambda l: 'Stop' not in l, possible_aa)
        degeneracy = float([GEN_CODE[ k ][ codon[1] ][ codon[2] ] for k in all].count(aa))
        degeneracy=degeneracy/len(possible_aa)
    elif position==2:
        all = filter(lambda l: codon[1] not in l, all)
        possible_aa = [GEN_CODE[ codon[0] ][ k ][ codon[2] ] for k in all]
        possible_aa = filter(lambda l: 'Stop' not in l, possible_aa)
        degeneracy = float([GEN_CODE[ codon[0] ][ k ][ codon[2] ] for k in all].count(aa))
        degeneracy=degeneracy/len(possible_aa)
    elif position==3:
        all = filter(lambda l: codon[2] not in l, all)
        possible_aa = [GEN_CODE[ codon[0] ][ codon[1] ][ k ] for k in all]
        possible_aa = filter(lambda l: 'Stop' not in l, possible_aa)
        degeneracy = float([GEN_CODE[ codon[0] ][ codon[1] ][ k ] for k in all].count(aa))
        degeneracy=float(degeneracy/len(possible_aa) )
    #print possible_aa, codon, aa, position, all, degeneracy
    return  degeneracy



def  translateCodon(codon):
    codonTable = { 'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C','TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C','TTA': 'L', 'TCA': 'S', 'TAA': 'X', 'TGA': 'X', 'TTG': 'L', 'TCG': 'S', 'TAG': 'X', 'TGG': 'W', 'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R','CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S', 'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S','ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G', 'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}
    if codon not in codonTable:
        #print "Error! ", codon,  "not in codonTable!\n"
        return None
    else:
        return codonTable[codon]


def reverse_complement ( s ):
    complement_dna = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "N":"N", "n":"n" }
    reversed_s = []
    for i in s:
        reversed_s.append( complement_dna[i] )
    reversed_s.reverse()
    return "".join( reversed_s )
