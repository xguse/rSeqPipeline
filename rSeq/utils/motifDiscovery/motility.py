import os
from rSeq.utils.errors import *

        
def mtlty_formatHits(seqName,seq,motifObj,ref='start',threshFrac=0.75):
    """Given SeqName, seq, and motility obj, return formated string of results referenced from
    the 'start' or the 'end' of the string.  Start position is 1-based and always
    most 5' bp with regard to the source sequence.  HitSequences are always
    displayed in same orientation as the motif used for the search and those that
    occur on the revComp strand at denoted with a '-1'."""
    
    goodRefs = ['start','end']
    if ref not in goodRefs:
        raise InvalidOptionError('ref argument must be in: %s' % (goodRefs))
    
    rawHits = motifObj.find(seq,(motifObj.max_score()*threshFrac))
    
    hits = []
    
    for h in rawHits:
        if ref == 'start':
            lEdge = h[0]+1
        else:
            lEdge = -(len(seq)-h[0])
        ori   = h[2]
        hSeq  = h[3]
        score = motifObj.calc_score(hSeq)/motifObj.max_score()
        hits.append((lEdge,ori,score,hSeq))

    hits.sort(key=lambda x: x[0])
    
    rStr = seqName
    for h in hits:
        rStr += ';(%s,%s,%.5g,%s)' % (h[0],h[1],h[2],h[3])
        
    return rStr
        
    