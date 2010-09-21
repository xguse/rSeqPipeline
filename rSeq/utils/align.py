import sys
from rSeq.utils.errors import *
from rSeq.utils.externals import runExternalApp



#def samCH2supCntg(samFile,outFile,cnvsnDict):
    #from gusPyCode.defs.mosqData.AaCntgCvrt import supContigConvert# FIX THIS!
    #"""Convert CHxxx to superContig1.xxx"""
    
    #outFile = open(outFile, 'w')
    
    #for line in open(samFile, 'rU'):
        #if line.startswith('@SQ'):
            #line = line.split('\t')
            #line[1] = line[1].split(':')[0],':',cnvsnDict[line[1].split(':')[1]]

def cigarStr2AlignCoords(cigarString,minCoord,maxCoord,alignStrand,cigKind='ensembl',intron='I'):
    """Returns List of coordinate tuples for each block in form of:
    [
    (blkMinCoord_1,blkMaxCoord_1),
    (blkMinCoord_2,blkMaxCoord_2),
    ...,
    (blkMinCoord_n,blkMaxCoord_n)
    ]
    
    NOTE: Assumes 0-referenced "in-between" coordinate system as used in BED/exonserate."""
    
    cigInfo = parseCigarString(cigarString,kind=cigKind)
    
    # Reverse cigData if alignStrand is neg
    if alignStrand == ('-' or '-1'):
        cigInfo.reverse()
    elif alignStrand == ('+' or '1'):
        pass
    else:
        raise InvalidOptionError("cigarStr2AlignCoords: valid alignStrand values are %s. You gave: %s" \
                                 % (['-','+','-1','1'],alignStrand))
    if intron == 'I':
        ignore = 'D'
    elif intron == 'D':
        ignore = 'I'
    else:
        raise InvalidOptionError("cigarStr2AlignCoords: valid intron values are %s. You gave: %s" \
                                 % (['I','D'],intron)) 
    
    # 1: remove all 'ignore' info since these seem to relate to insertions in the aligned seq
    # and do not affect coords on the chrm seq.
    cigNoIgnores = []
    for elem in cigInfo:
        if elem[0] == ignore:
            continue
        else:
            cigNoIgnores.append(elem)
            
    # 2: combine any 'M' info that are adjacent to each other after cleaning out the 'ignore' info.
    cigCmbnMs = []
    while 1:
        if not len(cigNoIgnores) > 0:
            break
        Ms = []
        while 1:
            try:
                if cigNoIgnores[0][0] == 'M':
                    Ms.append(cigNoIgnores.pop(0))
                else:
                    break
            except IndexError:
                break
        cigCmbnMs.append(('M',sum([x[1] for x in Ms])))
        if len(cigNoIgnores) > 0:
            cigCmbnMs.append(cigNoIgnores.pop(0))
        else:
            break
    
    # 3: walk the cig, keeping track of where on Chrm you are and logging boundaries.
    blkCoordsRelative = []
    i = 0
    while i < len(cigCmbnMs):
        if cigCmbnMs[i][0] == 'M':
            end   = sum([x[1] for x in cigCmbnMs[:i+1]]) 
            start = end - cigCmbnMs[i][1]
            blkCoordsRelative.append((start,end))
            i += 1
        else:
            i += 1
    blkCoordsAbsolute = []
    for coord in blkCoordsRelative:
        blkCoordsAbsolute.append((coord[0]+minCoord,coord[1]+minCoord))
    
    # Confirm that calculated max coord == maxCoord
    calcMax = max([x[1] for x in blkCoordsAbsolute])
    if calcMax != maxCoord:
        raise UnexpectedValueError("cigarStr2AlignCoords: calculated max coord != maxCoord.  My code may not be correct, but check that you used the correct 'intron' value first.")
    else:
        return blkCoordsAbsolute
    
def parseCigarString(cigarString,kind='ensembl'):
    """Returns cigarInfo from cigarString in form:
    [
    [operationLetter_1,run_1],
    ...,
    [operationLetter_n,run_n]
    ]
    
    Validates operation letters before returning."""
    
    validKinds  = ['ensembl','exonerate']
    validCigOps = ['M','D','I']
    if kind not in validKinds:
        raise InvalidOptionError('parseCigarString: valid kinds are: %s.  You gave: %s' % (validKinds,kind))
    if not type(cigarString) == type(''):
        raise InvalidOptionError('parseCigarString: type(cigarString) != type(""): %s' % (kind))
    
    cigData = None
    if kind == validKinds[0]:
        cigData = parseEnsemblCigar(cigarString)
    else:
        cigData = parseExonerateCigar(cigarString)
    
    # Do some validation of cigar operation letters:
    cigOps = list(set([x[0] for x in cigData]))
    for op in cigOps:
        if not op in validCigOps:
            raise InvalidOptionError("parseCigarString: valid cigar operation letters are %s.  malformed cigarString (%s)." \
                                     % (validCigOps,cigarString))
    return cigData
    

def parseExonerateCigar(cigarString):
    """Given exonerate style cigar string, return list of lists representing
    [[letter,number],[letter,number],...]."""
    cigarList = cigarString.split()
    
    # Format cigar info for easy use:
    cigar = []
    if not len(cigarList)%2==0:
        raise InvalidOptionError('parseExonerateCigar: cigarList var is not a list with an even number of indexes.')
    while cigarList:
        cigar.append((cigarList.pop(0),int(cigarList.pop(0))))
    return cigar

def parseEnsemblCigar(cigarString):
    """Given ensembl style cigar line, return list of lists representing
    [[letter,number],[letter,number],...]."""
    def getNumEnd(cigString,index):
        try:
            int(cigString[index])
            index+=1
            end = getNumEnd(cigString,index)
        except ValueError:
            return index
        return end
    
    cigList = []
    i=0
    while i < len(cigarString):
        num = None
        let = None
        
        try:
            int(cigarString[i])
            try:
                end = getNumEnd(cigarString,i)
                num = int(cigarString[i:end])
                i = end
            except ValueError:
                raise Exception('parseEnsemblCigar.getNumEnd: failed check code.')
        except ValueError:
            let = cigarString[i]
            i+=1
            
        if num:
            let = cigarString[i]
            i+=1
            cigList.append((let,num))
        elif let:
            cigList.append((let,1))
        else:
            raise InvalidOptionError('parseEnsemblCigar: passed through a cycle without assigning let OR num. Check for valid cigarString: %s' % (cigarString))
    return cigList
            

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
    cigInfo = ' '.join(line[10:])  # actual CIGAR encoding converted back to string for parser.
    
    # Format cigar info for easy use:
    cigar = parseCigarString(cigInfo,type='exonerate')
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
    
    