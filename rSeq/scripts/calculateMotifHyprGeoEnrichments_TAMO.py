#!/usr/bin/env python
import sys
import os
import optparse
import csv
import random as r
from time import time
import math
from matplotlib import pylab as pl
from TAMO.MotifMetrics import ProbeSet
from TAMO.MotifTools import load as tamoLoad
from rSeq.utils.motifDiscovery.motifs import getEvalHitDict,parseSCOPEfile,parseXMSfile,motifHyprGeoEnrichmentTAMO,toTAMOmotifs
from rSeq.utils.motifDiscovery.rPossum import *
from rSeq.utils.stats import seqStats
from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp
from rSeq.utils.misc import Bag
from rSeq.utils.files import ParseFastA


valideMotifTypes = {'scope':parseSCOPEfile,
                    'xms'  :parseXMSfile,
                    'tamo' :tamoLoad}

def getProbSet(fastaPath,factor):
    """Returns ProbeSet object with seqDict and methods."""
    return ProbeSet(fastaPath,factor)

def getSeqs(pSet):
    """Returns ref to SeqDict in ProbSet obj."""
    sDict  = pSet.probes
    return sDict

def getMotifs(pathToMotifFile,nucFreqs,fileType='scope'):
    """Return list of motif objects."""
    if fileType not in valideMotifTypes:
        raise InvalidOptionError()
    if fileType != "tamo":
        nativeMotifs = valideMotifTypes[fileType](pathToMotifFile)
        tMotifs      = toTAMOmotifs(nativeMotifs,seqData=nucFreqs)
    else:
        tMotifs      = valideMotifTypes[fileType](pathToMotifFile)
        if not type(nucFreqs) == type(dict):
                raise InvalidOptionError("**ERROR: in getMotifs(), nucFreqs must be type(dict) when motif type is TAMO**")
        for i in range(len(tMotifs)):
            tMotifs[i].background = nucFreqs
            tMotifs[i]._compute_ll()
            tMotifs[i].id = "%s_%s" % (tMotifs[i].__repr__()[:-4],tMotifs[i].barcode[:6])
            
    return tMotifs

def filterMotifs(motifList,key=None):
    """Return only motifs that meet the 'key' filter function.
    The default function is "lambda x: float(x.sigvalue) >= 5" """
    if not key:
        key = lambda x: (float(x.sigvalue) >= 5)
    mList = []
    
    for m in motifList:
        if key(m):
            #print key(m),' ',type(m.sigvalue),' ',m.sigvalue
            mList.append(m)
    return mList

def getSeqFreqs(seqDict):
    """"Returns DICT with correct values for TAMO
    AT_bias,GC_bias keywords as calculated from seqDict."""
    data = seqStats(seqDict)
    seqFreqs = {}
    seqFreqs['A'] = float(data['aCnt'])/data['nonNs']
    seqFreqs['C'] = float(data['cCnt'])/data['nonNs']
    seqFreqs['G'] = float(data['gCnt'])/data['nonNs']
    seqFreqs['T'] = float(data['tCnt'])/data['nonNs']

    return seqFreqs

def getForeground(geneList):
    """Return forground list"""
    return geneList.split(',')

def appendData(newData,outData):
    """Add new results to outData."""
    for m in newData:
        ID   = m[0]
        pVal = m[1]
        if ID in outData:
            outData[ID].append(pVal)
        else:
            outData[ID] = []
            outData[ID].append(pVal)
            
def getRandomForground(foregroundSeqs,seqDict):
    """Return list of randomly chosen geneNames of same
    length as foregroundSeqs."""
    return r.sample(seqDict.keys(),len(foregroundSeqs))

def getOutBaseStr(outPath,jobName):
    """Returns base string to use for writing out data."""
    if not outPath.startswith('/'):
        raise InvalidOptionError('--out: must be full path (%s)' (outPath))
    outBaseStr = '%s/%s' % (outPath.rstrip('/'),jobName)
    return outBaseStr

def writeDataTable(outBaseStr,outData,motifList):
    """Write out data table to file."""
    outFile = '%s.table.txt' % (outBaseStr)
    outFile = open(outFile, 'w')

    # --- write data ---
    ctrlLen = len(outData[outData.keys()[0]][1:])
    dataWriter = csv.writer(outFile, delimiter='\t')
    dataWriter.writerow(['#motifs','real']+['ctrl_%s'%(x) for x in range(ctrlLen)])
    for m in [x.id for x in motifList]:
        dataWriter.writerow([m]+outData[m])
    outFile.close()
    
def plotHist(outBaseStr,outData,motifID):
    """Plot histograms for each motif results."""
    def _chooseBins():
        """Uses one of a number of heuristics to 
        choose the bins for each hist.
        
        More discussed here:
        http://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width"""
        # sqrtChoice
        # k = int(sqrt(n))+1
        return int(math.sqrt(len(ctrls)))*3 

        
    outFile = '%s_%s_hist.pdf' % (outBaseStr,motifID)
    ctrls   = [x for x in outData[motifID][1:]]
    yLab    = 'Controls [n=%s]' % (len(ctrls))
    xLab    = 'Enrichment [p-value]'
    
    fig  = pl.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel(yLab)
    ax.set_xlabel(xLab)
    ax.set_xscale('log')
    
    bins = _chooseBins()
    hData = ax.hist(ctrls, bins=bins,histtype='stepfilled',color='grey',cumulative=1)
    ax.axvline(x=outData[motifID][0] ,linewidth=2, color='r')
    pl.savefig(outFile)


def makeMotifListFromPossum(pAccessions):
    """Returns a list of DummyPlug Motif objects derived from
    possum data that are compatible with how THIS script
    expects an object in 'motifList' to act."""
    mList = []
    for ac in pAccessions:
        mList.append(Bag({'accession':ac, 'id':ac}))
    return mList


# +-+-+-+-+ MAIN DEF +-+-+-+-+
def main():
        #+++++++++++ File Parseing Etc +++++++++++
    desc = """Calls the folowing funcs: 'add this'"""
    
    usage = """python %prog args"""
    parser = optparse.OptionParser(usage=usage, description=desc)
    
    parser.add_option('--motifs', type='str',default=None,
                      help="""Path to motif file (default=%default).""")
    parser.add_option('--motif-type', type='str',default='scope',
                      help="""Format of motif file (default=%default).""")
    parser.add_option('--motif-filter', type='str',default='lambda x: x',
                      help="""Filter fuction to allow filtering of motifs in motif file.  Its a lambda function.  The default returns all motifs. If you don't understand this, please leave it alone. (default=%default).""")
    parser.add_option('--thresh', type='float',default=0.75,
                      help="""Fractional score threshold for motif score cut-off (default=%default).""")
    parser.add_option('--promoters', type='str',default=None,
                      help="""Path to fasta file with promoter population (default=%default).""")
    parser.add_option('--plen', type='int',default=None,
                      help="""Max promoter length to use -- starting from 3'-end!! (default=%default).""")
    parser.add_option('--genes', type='string',default=None,
                      help="""Unbroken string of gene/Tx names representing the true forground set, sep=','. Exp: 'gene,gene,gene' (default=%default).""")
    parser.add_option('--job', type='string',default='int(time())',
                      help="""String to identify this run (default=%default).""")
    parser.add_option('--out', type='string',default=None,
                      help="""Path to results dir -- must use full path (default=%default).""")
    parser.add_option('--plot-fdrs', dest="fdrs", type='int',default=0,
                      help="""How many random sets of genes equal in length to --genes to run for FDR estimation. (default=%default).""")
    ##parser.add_option('--from-possum',default=False,
                      ##help="""Path to possumSearch outFile, skips motif finding step. (default=%default)""")    
    parser.add_option('--verbose',action='store_true',default=False,
                      help="""Include if stdOut/stdErr is desired. (default=%default).""")
    parser.add_option('--check-seqs',action='store_true',default=False,
                      help="""Print info about promoter sequences and exit. (default=%default).""")
    ##parser.add_option('--expect',action='store_true',default=False,
                      ##help="""Use median seqLen to set success threshold at greater than the estimated expected number of occurences in each promoter [pValThresh*searchesPerSeq]. (default=%default).""")
    

    
    (opts, args) = parser.parse_args()
    
    # +++++ Argument Validations +++++
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    if not opts.out:
        raise InvalidOptionError("--out argument is required.")
    if not opts.genes:
        raise InvalidOptionError('--genes argument is required.')
    if opts.job == 'int(time())':
        opts.job = int(time())
    if not (opts.motifs or opts.motif_type or opts.thresh or opts.promoters):
        raise InvalidOptionError("""Unless --from-possum, the following argument are ALL required:
        --motif-type
        --thresh""")
    if not opts.motif_filter.startswith('lambda x:'):
        raise InvalidOptionError("**ERROR: the --motif-function option must begin with 'lambda x:'**")
    else:
        opts.motif_filter = eval(opts.motif_filter)
    if opts.plen:
        try:
            opts.plen = int(opts.plen)
        except ValueError:
            raise InvalidOptionError("**ERROR: the --plen option must be a number**")
    
    # +++++ Lets Begin +++++
    if opts.verbose: sys.stdout.write('\n%s\n\n' % (' '.join(sys.argv)))
    
    if opts.verbose: sys.stdout.write('preparing out directory...\n')
    outBaseStr = getOutBaseStr(opts.out,opts.job)
    mkdirp(opts.out)
    
    if opts.verbose: sys.stdout.write('building seqDict...\n')
    probeset = getProbSet(opts.promoters,opts.thresh)
    #print len(probeset.probes.items()[0][1])
    #print probeset.probes.items()[0][1]
    if opts.plen:
        if opts.verbose: sys.stdout.write("adjusting promoter lengths to no more than %s and conserving the 3' ends...\n" % (opts.plen))
        for s in probeset.probes.iteritems():
            probeset.probes[s[0]]= s[1][-opts.plen:]
        #print len(probeset.probes.items()[0][1])
        #print probeset.probes.items()[0][1]
            
    seqDict  = getSeqs(probeset)
    if opts.check_seqs:
        seqStats(seqDict,show=True)
        exit(0)
    
    if opts.verbose: sys.stdout.write('getting nucFreqs...\n')
    nucFreqs  = getSeqFreqs(seqDict)    

    if opts.verbose: sys.stdout.write('building filtered motifList...\n')
    motifList = getMotifs(opts.motifs,nucFreqs,opts.motif_type)
    preFilt   = len(motifList)
    motifList = filterMotifs(motifList,key=opts.motif_filter)
    postFilt  = len(motifList)
    if opts.verbose: sys.stdout.write('using %s of %s motifs...\n' % (postFilt,preFilt))
    
    if opts.verbose: sys.stdout.write('getting forgroundSeqs...\n')
    foregroundSeqs = getForeground(opts.genes)
    outData        = {}
    
    if opts.verbose: sys.stdout.write('calculating real data p-values...\n')
    realData = motifHyprGeoEnrichmentTAMO(motifList,probeset,foregroundSeqs,factor=opts.thresh,bestFactor=False)
    appendData(realData,outData)
    
    
    for i in range(opts.fdrs):
        if opts.verbose: sys.stdout.write('ctrl_%s:\n' % (i))
        if opts.verbose: sys.stdout.write('\tgetting random forground...\n')
        rForground = getRandomForground(foregroundSeqs,seqDict)
        
        if opts.verbose: sys.stdout.write('\tcalculating p-values...\n')
        ctrlData   = motifHyprGeoEnrichmentTAMO(motifList,probeset,rForground,factor=opts.thresh,bestFactor=False)
        appendData(ctrlData,outData)
    
    if opts.verbose: sys.stdout.write('writting outData...\n')
    writeDataTable(outBaseStr,outData,motifList)
    
    
    if opts.verbose: sys.stdout.write('plotting histograms...\n')
    for m in motifList:
        plotHist(outBaseStr,outData,m.id)
    
if __name__ == "__main__":

    main()

        
        
