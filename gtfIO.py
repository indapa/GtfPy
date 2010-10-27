from Gtf import GtfRec
def getCds(fh):
    """ python generator that yields Gtf object for CDS feature """
    for line in fh:
        if '#' not in line:
            (name,src,feat,begin, end, score, strand, frame, attr) = line.strip().split('\t')
            if 'CDS' == feat:
                cdsGtf=GtfRec(name,src,feat,begin, end, score, strand, frame, attr)
                yield cdsGtf

def getStartCodon(fh):
    """ python generator that yields Gtf object for start_codon feature """
    for line in fh:
        if '#' not in line:
            (name,src,feat,begin, end, score, strand, frame, attr) = line.strip().split('\t')
            if 'start_codon' == feat:
                start_codonGtf=GtfRec(name,src,feat,begin, end, score, strand, frame, attr)
                yield start_codonGtf

def getStopCodon(fh):
    """ python generator that yields Gtf object for stop_codon feature """
    for line in fh:
        if '#' not in line:
            (name,src,feat,begin, end, score, strand, frame, attr) = line.strip().split('\t')
            if 'stop_codon' == feat:
                stop_codonGtf=GtfRec(name,src,feat,begin, end, score, strand, frame, attr)
                yield stop_codonGtf


def getUtr(fh):
    """ python generator that yields Gtf object for utr feature """
    for line in fh:
        if '#' not in line:
            (name,src,feat,begin, end, score, strand, frame, attr) = line.strip().split('\t')
            if 'UTR' == feat:
                featutreGtf=GtfRec(name,src,feat,begin, end, score, strand, frame, attr)
                yield featureGtf


def getExon(fh):
    """ python generator that yields Gtf object for exon feature """
    for line in fh:
        if '#' not in line:
            (name,src,feat,begin, end, score, strand, frame, attr) = line.strip().split('\t')
            if 'exon' == feat:
                featutreGtf=GtfRec(name,src,feat,begin, end, score, strand, frame, attr)
                yield featureGtf


def getExon(fh):
    """ python generator that yields Gtf object for gene feature """
    for line in fh:
        if '#' not in line:
            (name,src,feat,begin, end, score, strand, frame, attr) = line.strip().split('\t')
            if 'gene' in feat:
                featutreGtf=GtfRec(name,src,feat,begin, end, score, strand, frame, attr)
                yield featureGtf

def bed_get_position(fh):
    """ get the first three fields of a bed file """
    for line in fh:
            if '#' in line: continue
            datafields=line.strip().split('\t')
            yield datafields[0], int(datafields[1]), int(datafields[2])