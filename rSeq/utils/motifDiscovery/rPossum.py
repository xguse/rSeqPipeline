from collections import namedtuple
from rSeq.utils.errors import *
from rSeq.utils.files import tableFile2namedTuple

# +++++ Useful Constants +++++
possumHeaders = ['id','ac','de','group_id',
                 'group_pos','match_start',
                 'match_len','ori','thresh_used',
                 'match_score','min_pssm_score',
                 'max_pssm_score','pvalue',
                 'evalue','mss','seq_no',
                 'seq_de','match_str']


# +++++  Helper Defs +++++
def getPossumProfileACs(path):
    """Returns a tuple of possum accessions from a 
    possum profile file."""
    pssmData = open(path,'rU').readlines()
    ACs      = []
    for line in pssmData:
        if line.startswith('AC '):
            ACs.append(line.strip('\n').split()[1])
    return ACs

def getPossumHitTable(path,headers=possumHeaders):
    """Return a tuple of named tuples representing the results
    of a PoSSuMSearch run."""
    return tuple(tableFile2namedTuple(path,sep='\t',headers=headers))

def getPossumHitDict(possumHitTable,allSeqNames,pAccessions):
    """Returns hitDict compatible with motifs.motifHyprGeoEnrichment()
    and a dict as a key to which motifs are represented in which indexes.
    The success values are 0 (absent), and >=1 (present) at least once.
    
    mNum = number of motifs searched."""
    #sanitize pAccessions
    for i in range(len(pAccessions)):
        pAccessions[i] = pAccessions[i].replace('.','_')
        
    MotifTup = namedtuple('MotifTup','%s'%(' '.join(pAccessions)))
    
    hDict = {}

    for seqName in allSeqNames:
        hDict[seqName] = [0]*len(pAccessions)
        
    
    for row in possumHitTable:
        seqName = row.seq_de.split('.')[0]
        if seqName in hDict:
            gPos = int(row.group_pos)
            hDict[seqName][gPos] += 1
        else:
            hDict[seqName] = [0]*mNum
            gPos = int(row.group_pos)
            hDict[seqName][gPos] += 1

    for seqName in hDict:
        hDict[seqName] = MotifTup._make(hDict[seqName])

    return hDict
    