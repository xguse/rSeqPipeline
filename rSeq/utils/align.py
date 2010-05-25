import sys
from rSeq.utils.errors import *
from rSeq.utils.externals import runExternalApp


def exonerateCigar2BEDline(cigarLine,rgb="0,0,0"):
    """Takes a cigar line from an exonerate output file and returns
    a BED formated line representing the alignment blocks of the query
    as aligned to the target. """
    
    line  = cigarLine.strip('\n').split() # its a space delimited line
    
    # exonerate data
    query   = line[1]   # query id
    qStart  = line[2]   # ! will be lower number than qEnd if qStrand is '-'
    qEnd    = line[3]   # ! will be higher number than qStart if qStrand is '-'
    qStrand = line[4]   #
    target  = line[5]   # target id
    tStart  = line[6]   # ! will be lower number than tEnd if tStrand is '-'
    tEnd    = line[7]   # ! will be higher number than tStart if tStrand is '-'
    tStrand = line[8]   #
    score   = line[9]   # "sum of transistion scores used in the dynamic programming."
    cigInfo = line[10:]  # actual CIGAR encoding
    
    # Format cigar info for easy use:
    cigar = []
    if not len(cigInfo)%2==0:
        raise UnexpectedValueError('cigInfo var is not a list with an even number of indexes.')
    while cigInfo:
        cigar.append((cigInfo.pop(0),cigInfo.pop(0)))
    if "I" in [x[0] for x in cigar]:
        sys.stderr.write('warn: just encountered a cigar containing "I".  Skipping.\n')
        return '' # for now just do nothing

        #raise UnexpectedValueError('cigar info includeds "I" and I have not been taught how to deal with this!')
    
    # If match is to '-' strand of target, reverse the cigar info
    if tStrand == '-':
        cigar.reverse()
        
    # Bed Terms
    chrom       = target
    chromStart  = min(int(tStart),int(tEnd))
    chromEnd    = max(int(tStart),int(tEnd))
    name        = query
    score       = 0
    strand      = tStrand
    thickStart  = chromStart
    thickEnd    = chromEnd
    itemRgb     = rgb
    blockCount  = len([x[0] for x in cigar if x[0]=="M"])
    blockSizes  = ','.join([x[1] for x in cigar if x[0]=="M"])
    blockStarts = ['0']  # for now.  See below.
    
    # Calculate blockStarts
    i = 0
    while 1:
        if i < len(cigar)-1:
            blockStarts.append(str(int(cigar[i][1])+int(cigar[i+1][1])))
            i += 2
        else:
            break
    blockStarts = ','.join(blockStarts)
    
    # Formulate and return BedLine
    return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
           (chrom,chromStart,chromEnd,name,score,strand,
            thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts)
    
def bowtieMakeIndex(optionsHere):
    """Create bowtie indexes from new fasta set."""
    
def bowtieAlign(optionsHere):
    """Run alignment of fastQ to bowtie index."""
    
    