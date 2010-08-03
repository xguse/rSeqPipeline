#import pyRXP
#import xmlutils
from xml.dom import minidom
from collections import namedtuple
import exceptions
import os
import hashlib
from motility import IUPAC,PWM
import motility
from cogent import DNA,SequenceCollection,Alignment
from rSeq.utils.files import ParseFastA
from rSeq.utils.stats import cumHypergeoP,seqStats
from rSeq.utils.errors import *

# ++++++++ useful constants ++++++++ 
ltr2wrd = {'A':'adenine',
           'C':'cytosine',
           'G':'guanine',
           'T':'thymine'}



# ++++++++ classes ++++++++ 
class Motif(object,IUPAC,PWM):
    """Represents the information attributed to a single motif matrix.

    """
    
    def __init__(self,motifData, mType='xms',rank=None, source=None):
        """Fill in soon.
        """
        # define some constants:
        self.vaild_mTypes = {'scope':self._initSCOPE,
                             'xms':self._initXMS,
                             'IUPAC':self._initIUPAC,
                             'PWM':self._initMotPWM,
                             'jaspar':self._initJASPAR}

        self.bases = ('A','C','G','T')

        # Basic instance attributes:
        #     Other possible attr: pVal,eVal,sigVal,genes,rank
        self.id          = None
        self.rank        = rank
        self.name        = None
        self.accession   = None
        self.consensus   = None
        self.source      = source
        self.barcode     = None  # hash of pwm data to fingerprint this matrix.
        self.instances   = ()    # tuple of namedtuples keys=(sequence, strand<-/+>,start[negVersion],end[negVersion],geneName))
        self.pwm         = {}    # key(nulceotide),value(tuple of freqs at each position)
        self.threshPval  = None
        self.threshPval_target = None
        self.threshScore = None
        self.kmersOverThresh =None

        # init instance attributes with correct function
        self.vaild_mTypes[mType](motifData)


    # +++++ general instance methods +++++
    def _setID(self):
        """Sets self.id after consensus and barcode have been defined."""
        if self.source == 'scope':
            self.id = '%s_%' % (self.consensus,int(float(self.sigvalue)))
        elif self.source == 'jaspar':
            self.id = '%s_%s' % (self.accession.replace('.','v'),self.name)
        else:
            self.id = '%s_%s' % (self.consensus,self.barcode[:6])
        
    def _getMotilityMatrix(self,pwm):
        """Return matrix expected by motility from
        self.pwm."""
        moMat = []
        for i in range(len(pwm['A'])):
            pos_i = []
            for b in sorted(self.bases):
                pos_i.append(pwm[b][i])
            moMat.append(tuple(pos_i))
        return tuple(moMat)

    def _getFreq(self,col,base):
        """get fractional weight for a base at given column."""
        counts = [float(self.pwm['A'][col]),
                  float(self.pwm['C'][col]),
                  float(self.pwm['G'][col]),
                  float(self.pwm['T'][col])]
        tot = sum(counts)
        return float(self.pwm[base][col])/tot
    
    def _generateConsensus(self):
        """Trick motility into making our consensus for us if need be.
        You must have already called a motility __init__ to do this.
        Specifically, self.matrix must exist."""
        return motility.make_iupac_motif(self.generate_sites_over(self.max_score()*0.99))
        
    def _initJASPAR(self,motifData):
        """Sets Motif's attributes using a JASPAR motif record."""
        
        def _buildPWM(dataString):
            """Convert jaspars's matrix data to my format and
            return the PWM/barcode. 
            dataString = 'A  [0  3]\nC  [94 75]\nG  [1  0]\nT  [2 19]"""
            
            dataString = dataString.replace('[',' ')
            dataString = dataString.replace(']',' ')
            data       = dataString.split('\n')
            data       = [x.split() for x in data]

            pwm = {}
            for d in data:
                pwm[d[0]] = tuple([int(x) for x in d[1:]])
            
            bc = hashlib.md5(str(pwm)).hexdigest()
            return pwm,bc
        
        acc,mName = motifData[0].lstrip('>').rstrip('\n').split()
        self.accession = acc
        self.name      = mName
        self.pwm,self.barcode = _buildPWM(motifData[1]) 
        PWM.__init__(self,self._getMotilityMatrix(self.pwm)) # call motility's PWM's init
        self.consensus = self._generateConsensus()
        self.source = 'jaspar'
        self._setID()
    
    def _initMotPWM(self,motilityStylePWM):
        """Sets Motif's attributes using a motilityStylePWM."""
        raise NotImplementedError()
        
    def _initIUPAC(self,consensus):
        """Sets Motif's attributes using an IUPAC motif string."""
        IUPAC.__init__(self,consensus)
        
        def _motil2myMat():
            """Convert motility's matrix data to my format and
            return the PWM/barcode."""
            pwm = {}
            bases = sorted(self.bases)
            for i in range(4):
                baseData = []
                for pos in self.matrix:
                    baseData.append(pos[i])
                pwm[bases[i]] = tuple(baseData)
            bc = hashlib.md5(str(pwm)).hexdigest()
            return pwm,bc
        
        self.consensus = self.motif
        self.pwm,self.barcode = _motil2myMat()
        self._setID()
        self.accession = self.id
        
    def _initXMS(self,motifData):
        """Sets Motif's attributes using a single 'motif' node from a
        minidom object representing the XMS file."""

        def _xmsBaseWeights(base,motifData):
            """Return list of base counts/freqs for all columns
            of motif."""
            baseWeights = []
            weights  = motifData.getElementsByTagName('weight')
            for w in weights:
                if w.attributes.getNamedItem('symbol').value.upper().startswith(base):
                    baseWeights.append(float(w.childNodes[0].data))
            return tuple(baseWeights)

        def _xmsProp(motifData,propKey):
            """Returns stored property values if they exist from XMS motif."""
            props     = motifData.getElementsByTagName('prop')
            propVal = None
            for p in props:
                if p.getElementsByTagName('key')[0].childNodes[0].data.lower() == propKey.lower():
                    propVal = p.getElementsByTagName('value')[0].childNodes[0].data
                    break
            return propVal

        def _xmsPWM(motifData):
            """Returns pwm and barcode data from xms motif."""
            pwm = {}
            for base in self.bases:
                pwm[base] = _xmsBaseWeights(base,motifData)
            bc = hashlib.md5(str(pwm)).hexdigest()
            return pwm,bc

        # init instance attributes from xms motifData
        self.name       = _xmsProp(motifData,'name')
        self.accession  = _xmsProp(motifData,'accession')
        self.consensus  = _xmsProp(motifData,'consensus')
        self.sigvalue   = _xmsProp(motifData,'sigvalue')
        self.rank       = _xmsProp(motifData,'rank')
        self.algorithm  = _xmsProp(motifData,'algorithm')

        self.pwm,self.barcode = _xmsPWM(motifData)
        self._setID()
        PWM.__init__(self,self._getMotilityMatrix(self.pwm))

    def _initSCOPE(self,motifData):
        """Sets Motif's attributes using a single SCOPE 'motif' data set."""

        def _scopeBaseWeights(base,motifData):
            """Return list of base counts/freqs for all columns
            of motif."""
            baseWeights = []
            bases  = motifData.getElementsByTagName('base')
            for b in bases:
                if b.attributes.getNamedItem('name').value.upper().startswith(base):
                    for w in b.getElementsByTagName('weight'):
                        baseWeights.append(float(w.childNodes[0].data))
                        break
            return tuple(baseWeights)

        def _scopePWM(motifData):
            """Returns pwm and barcode data from scope motif."""
            pwm = {}
            for base in self.bases:
                pwm[base] = _scopeBaseWeights(base,motifData)
            bc = hashlib.md5(str(pwm)).hexdigest()
            return pwm,bc

        def _scopeInstances(motifData):
            """Returns tuple of namedtuples:
            keys=(sequence, strand<-/+>,start[negVersion],end[negVersion],geneName))"""
            # Build instances namedtuples:
            Instance = namedtuple('instance','sequence strand begin end gene')
            instList = [Instance(x.getAttribute('sequence'),
                                 x.getAttribute('strand'),
                                 x.getAttribute('begin'),
                                 x.getAttribute('end'),
                                 x.getAttribute('gene')) for x in motifData.getElementsByTagName('instance')]
            return tuple(instList)

        self.consensus  = motifData.getAttribute('sequence')
        self.sigvalue   = motifData.getAttribute('sigvalue')
        self.algorithm  = motifData.getAttribute('algorithm')
        self.genes      = tuple(sorted(list(set([x.getAttribute('gene') for x in motifData.getElementsByTagName('instance')])))) # (sorted tuple)
        self.instances  = _scopeInstances(motifData) # tuple of namedtuples keys=(sequence, strand<-/+>,start[negVersion],end[negVersion],geneName))
        self.name       = '%s_%s' % (self.sigvalue,self.consensus)
        self.accession  = '%s_%s' % (self.consensus,int(float(self.sigvalue)))
        self.pwm,self.barcode = _scopePWM(motifData) # key(nulceotide),value(tuple of freqs at each position)
        self._setID()
        PWM.__init__(self,self._getMotilityMatrix(self.pwm))

    def toPossumMotif(self):
        """Return String representing the PossumSearch form of this motif.
        """
        # ++ initialize possum motif string ++
        if type(self.pwm['A'][0]) == type(1):
            pType = 'INT'
        elif type(self.pwm['A'][0]) == type(1.0):
            pType = 'FLOAT'

        pMtf  = 'BEGIN %s\n' % (pType)
        pMtf += 'ID %s\n'    % (self.name)
        pMtf += 'AC %s\n' % (self.accession)
        pMtf += 'DE %s\n'    % (self.consensus)
        pMtf += 'AL ACGT\n'
        pMtf += 'LE %s\n'    % (len(self.pwm['A']))

        # ++ create and add each column's data ++
        for i in range(len(self.pwm['A'])):
            pMtf += 'MA '
            for base in sorted(self.pwm.keys()):
                pMtf += '%.6g ' % (self.pwm[base][i])
            pMtf += '\n'                    
        pMtf += 'END\n'
        return pMtf


    def toXMSmotif(self):
        """Return string representing the XMS motif element for
        this motif.
        """
        # ++ initialize xms motif string ++
        xMtf = '<motif>\n\t<name>%s_%.3f</name>\n\t\t<weightmatrix alphabet="DNA" columns="%s">\n' %\
             (self.consensus,float(self.sigvalue),len(self.pwm['A']))

        # ++ create and add each column's data ++                  
        for i in range(len(self.pwm['A'])):
            xMtf += '\t\t\t<column pos="%s">\n' % (i)
            for base in sorted(self.pwm.keys()):
                freq = self._getFreq(i,base)
                xMtf += '\t\t\t\t<weight symbol="%s">%.4f</weight>\n' % \
                     (ltr2wrd[base],freq)
            xMtf += "\t\t\t</column>\n"
        xMtf += '\t\t</weightmatrix>\n'
        xMtf += '\t\t<threshold>0</threshold>\n' # for now we'll just use null threshold
        xMtf += '\t\t<prop>\n\t\t\t<key>consensus</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.consensus)
        xMtf += '\t\t<prop>\n\t\t\t<key>sigvalue</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.sigvalue)
        xMtf += '\t\t<prop>\n\t\t\t<key>rank</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.rank)
        xMtf += '\t\t<prop>\n\t\t\t<key>algorithm</key>\n\t\t\t<value>%s</value></prop>\n' % \
             (self.algorithm)
        xMtf += '</motif>\n'

        return xMtf


    def toTAMO(self,bkFreqs=None):
        """Return TAMO motif Object from self.pwm
        """
        raise exceptions.NotImplementedError()

    def toTRANSFACstring(self):
        """Returns printable motif representation in TRANSFAC format sutible for use in MAST.
        """
        raise exceptions.NotImplementedError()

    def scan(self,dnaSeq,pValThresh=0.01,halfAT=0.25,halfGC=0.25,log=False):
        """Return motility hits tuple.  You should use this instead of 'find'."""
        def _setThresh(newPval):
            """Return threshold that is approximately
            exquivilent to a p-value of 'pValThresh'"""
            maxScore = self.max_score()
            minScore = self.min_score() # to be used for lower bound (but later.. :( )
            scoreSteps = [maxScore*(x/100.0) for x in range(1,101)]
            scoreSteps.reverse()
            thresh = maxScore
            for score in scoreSteps:
                pVal = self.weight_sites_over(score,AT_bias=halfAT, GC_bias=halfGC)
                if pVal <= newPval:
                    thresh = score
                elif pVal > newPval:
                    break
            threshPval = self.weight_sites_over(thresh,AT_bias=halfAT, GC_bias=halfGC)
            if threshPval <= newPval:
                return thresh,threshPval
            else:
                return None

        if (self.threshPval == pValThresh) or (self.threshPval_target == pValThresh):
            #print 'No change in Thresh.'
            pass
        else:
            self.threshPval_target = pValThresh
            newThreshs = _setThresh(pValThresh)
            if not newThreshs:
                return False
            else:
                self.threshScore,self.threshPval = _setThresh(pValThresh)
                self.kmersOverThresh = len(self.generate_sites_over(self.threshScore))
                #print '%s new score: %s pVal: %s  kmersOver: %s' % (self.id,
                                                                    #self.threshScore,
                                                                    #self.threshPval,
                                                                    #self.kmersOverThresh)
        if self.kmersOverThresh < 1:
            return False
        else:
            return self.find(dnaSeq,self.threshScore)
        



# ++++++++ meta functions ++++++++ 

# --- read and write motif files ---
def parseJASPARfile(filePath):
    """Returns list of Motif Objects derived from filePath."""

    mList = []
    # +++ JASPAR files are pseudo-FastA style +++
    parser = ParseFastA(filePath,joinWith='\n',key=lambda x:x)
    for rec in parser:
        mList.append(Motif(rec,mType='jaspar'))
    
    return mList
        

def parseSCOPEfile(filePath):
    """Returns list of Motif Objects derived from filePath."""
    raise exceptions.NotImplementedError()

def parseXMSfile(filePath):
    """Returns list of Motif Objects derived from filePath."""

    mList = []

    # Get XMS data into dom obj
    xms = minidom.parse(filePath)
    motifs = xms.getElementsByTagName('motif')

    for m in motifs:
        mList.append(Motif(m, mType='xms'))

    return mList

def writeXMSfile(motifList,outPath):
    """Given a list of Motif objs, write out the contents in XMS format."""
    oFile = open(outPath, 'w')
    oFile.write('<motifset>\n')
    for m in motifList:
        oFile.write('%s' % (m.toXMSmotif()))
    oFile.write('</motifset>\n')
    oFile.close()

def writePoSSuMfile(motifList,outPath):
    """Given a list of Motif objs, write out the contents in PoSSuM format."""
    oFile = open(outPath, 'w')
    oFile.write('BEGIN GROUP\n')
    for m in motifList:
        oFile.write('%s' % (m.toPossumMotif()))
    oFile.write('END')
    oFile.close()
    
    
    
# --- facilitate analyses with motifs ---
def toTAMOmotifs(motifList,seqData=None):
    """Returns a list of TAMO motif Objects when given
    a list of the native motif Objs. 
    If seqData is a dict, it is assumed to be
    nucFreqs (seqData['A']==.26),
    else, it is assumed to be a path to a seqPopulation
    that the motifs will be used on for stats to be
    calculated against."""
    
    from TAMO.MotifTools import Motif_from_counts
    if not seqData:
        if type(seqData) == type({}):
            try:
                aFreq = seqData['A']
                cFreq = seqData['C']
                gFreq = seqData['G']
                tFreq = seqData['T']
            except KeyError:
                raise InvalidOptionError("toTAMOmotifs: Unrecognized key in seqData as dict.")
        elif type(seqData) == type(''):
            sDict = ParseFastA(seqData).toDict()
            seqData = seqStats(sDict)
            del(sDict)
            tot = float(seqData['nonNs'])
            aFreq = seqData['aCnt']/tot
            cFreq = seqData['cCnt']/tot
            gFreq = seqData['gCnt']/tot
            tFreq = seqData['tCnt']/tot
    else:
        aFreq = 0.25
        cFreq = 0.25
        gFreq = 0.25
        tFreq = 0.25
    
    tList = []
    for m in motifList:
        counts = []
        for i in range(len(m)):
            counts.append({'A':m.pwm['A'][i],'C':m.pwm['C'][i],'G':m.pwm['G'][i],'T':m.pwm['T'][i]})
        t = Motif_from_counts(counts,beta=0.01,bg={'A':aFreq,'C':cFreq,'G':gFreq,'T':tFreq})
        tList.append(t)
    return tList
    

def getEvalHitDict(motifList,seqDict,pThresh=0.01,halfAT=0.25,halfGC=0.25):
    """Return a dict of which motifs return more hits at given threshold than expected by
    the pval*kmers searched (E-val). 
    
    Keys=geneNames, Values=namedTuple of results named by each motif's ID"""
    
    hDict = {}
    
    # define named tuple
    MotifTup = namedtuple('MotifTup','%s'%(' '.join([x.id for x in motifList])))
    for seq in seqDict:
        motifData = []
        for m in motifList:
            hits = m.scan(seqDict[seq],pValThresh=pThresh,halfAT=halfAT,halfGC=halfGC)
            if hits == False:
                motifData.append(hits)
            else:
                seqLen      = len(seqDict[seq])
                kmersTested = (seqLen-len(m)+1)*2 # we search on both strands
                Eval        = pThresh*kmersTested
                if len(hits) > Eval:
                    motifData.append(len(hits))
                    #print 'Yay! %s' % (m.id)
                else:
                    motifData.append(None)
        hDict[seq] = MotifTup._make(motifData)
    return hDict
                
            
def motifHyprGeoEnrichment(motifList,hitDict,foregroundSeqs,expect=0):
    """Calculate the hypergeometric enrichment p-Value of each motif in motifList.
    n = # of positives in population
    i = # of positives in sample
    m = # of negatives in population
    N = sample size

    P(i) = sum([as i->N] (choose(n,i)choose(m,N-i))/choose(n+m,N))"""
    
    def _get_n_and_i(motifID,expect):
        """Number of positives in population and sample."""
        n = 0
        i = 0
        for item in hitDict.iteritems():
            if item[1].__getattribute__(motifID) > expect:
                if item[0] in foregroundSeqs:
                    n += 1
                    i += 1
                else:
                    n += 1
        return n,i

    def _get_m():
        """Number of negatives in population"""
        return len(hitDict)-n
        
    def _get_N():
        """sample size"""
        return len(foregroundSeqs)
    
    foregroundSeqs = set(foregroundSeqs)
    pVals = []
    for mtf in motifList:
        n,i = _get_n_and_i(mtf.id,expect)
        m = _get_m()
        N = _get_N()
        pVal = cumHypergeoP(n,i,m,N)
        pVals.append(tuple([mtf.id,cumHypergeoP(n,i,m,N)]))
        
    return tuple(pVals)
    

    
