import sys
import os
import optparse
import csv
import random as r
from time import time
import math
from cogent import LoadSeqs
from matplotlib import pylab as pl
from rSeq.utils.motifDiscovery.motifs import getEvalHitDict,parseSCOPEfile,parseXMSfile,motifHyprGeoEnrichment
from rSeq.utils.motifDiscovery.rPossum import *
from rSeq.utils.stats import seqStats
from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp
from rSeq.utils.misc import Bag


valideMotifTypes = {'scope':parseSCOPEfile,
                    'xms':parseXMSfile}
    

def getSeqs(fastaPath,pLen):
    """Returns seqDict with correct lengths."""
    cogSeqs = LoadSeqs(fastaPath,aligned=False, label_to_name=lambda x: x.split()[0])
    sDict = cogSeqs.todict()
    
    for s in sDict.iteritems():
        sDict[s[0]]= s[1][-pLen:]
        
    return sDict

def getMotifs(pathToMotifFile,fileType='scope'):
    """Return list of motif objects."""
    if fileType not in valideMotifTypes:
        raise InvalidOptionError()
    return valideMotifTypes[fileType](pathToMotifFile)

def getSeqFreqs(seqDict):
    """"Returns tuple with correct values for motility
    AT_bias,GC_bias keywords as calculated from seqDict."""
    
    # These look weird bc motility expects these info in a weird way.
    # I have confirmed from Titus Brown (motility author) that this is correct.
    GC_bias = seqStats(seqDict)['percentGC']/2
    AT_bias = 0.5 - GC_bias
    return AT_bias,GC_bias

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
    

    
if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    desc = """Calls the folowing funcs: 'add this'"""
    
    usage = """python %prog args"""
    parser = optparse.OptionParser(usage=usage, description=desc)
    
    parser.add_option('--motifs', type='str',default=None,
                      help="""Path to motif file (default=%default).""")
    parser.add_option('--motif-type', type='str',default='scope',
                      help="""Format of motif file (default=%default).""")
    parser.add_option('--thresh', type='float',default=0.005,
                      help="""P-value threshold for motif score cut-off (default=%default).""")
    parser.add_option('--promoters', type='str',default=None,
                      help="""Path to fasta file with promoter population (default=%default).""")
    parser.add_option('--plen', type='int',default=1000,
                      help="""Max promoter length to use -- starting from 3'-end!! (default=%default).""")
    parser.add_option('--genes', type='string',default=None,
                      help="""Unbroken string of gene/Tx names representing the true forground set, sep=','. Exp: 'gene,gene,gene' (default=%default).""")
    parser.add_option('--job', type='string',default='int(time())',
                      help="""String to identify this run (default=%default).""")
    parser.add_option('--out', type='string',default=None,
                      help="""Path to results dir -- must use full path (default=%default).""")
    parser.add_option('--plot-fdrs', dest="fdrs", type='int',default=0,
                      help="""How many random sets of genes equal in length to --genes to run for FDR estimation. (default=%default).""")
    parser.add_option('--from-possum',default=False,
                      help="""Path to possumSearch outFile, skips motif finding step. (default=%default)""")    
    parser.add_option('--verbose',action='store_true',default=False,
                      help="""Include if stdOut/stdErr is desired. (default=%default).""")
    parser.add_option('--check-seqs',action='store_true',default=False,
                      help="""Print info about promoter sequences and exit. (default=%default).""")
    parser.add_option('--expect',default=False,
                      help="""Use median seqLen to set success threshold at greater than the estimated expected number of occurences in each promoter [pValThresh*searchesPerSeq]. (default=%default).""")
    

    
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
    if not opts.from_possum:
        if not opts.motifs or opts.motif_type or opts.thresh or opts.promoters or opts.plen:
            raise InvalidOptionError("""Unless --from-possum, the following argument are ALL required:
            --motif-type
            --thresh
            --plen""")
    else:
        if not opts.motifs:
            raise InvalidOptionError("When using --from-possum, --motifs should be the PSSM file used to generate this particular hit set.")
    
    # +++++ Lets Begin +++++
    if opts.verbose: sys.stdout.write('\n%s\n\n' % (' '.join(sys.argv)))
    
    if opts.verbose: sys.stdout.write('preparing out directory...\n')
    outBaseStr = getOutBaseStr(opts.out,opts.job)
    mkdirp(opts.out)
    
    if opts.verbose: sys.stdout.write('building seqDict...\n')
    seqDict = getSeqs(opts.promoters,opts.plen)
    if opts.check_seqs:
        seqStats(seqDict,show=True)
        exit(0)
    else:
        seqData = seqStats(seqDict,show=False)
    if opts.expect:
        opts.expect = seqData['medLen']*opts.thresh*2
    else:
        opts.expect = 0
        
    # --- am i doing the searching myself? ---
    if not opts.from_possum:
        # -- yes --
    
        if opts.verbose: sys.stdout.write('building motifList...\n')
        motifList      = getMotifs(opts.motifs,opts.motif_type)
    
        if opts.verbose: sys.stdout.write('getting nucFreqs...\n')
        halfAT,halfGC  = getSeqFreqs(seqDict)
    
        if opts.verbose: sys.stdout.write('building hitDict...\n')
        motifHits      = getEvalHitDict(motifList,seqDict,pThresh=opts.thresh,halfAT=halfAT,halfGC=halfGC)
    else:
        # -- Oh thank god, no!  All I have to do is some parseing! --
        if opts.verbose: sys.stdout.write('skipping to building hitDict step...\n')
        pACs        = getPossumProfileACs(opts.motifs)
        possumTable = getPossumHitTable(opts.from_possum,headers=possumHeaders)
        motifHits   = getPossumHitDict(possumTable,seqDict.keys(),pACs)
        motifList   = makeMotifListFromPossum(pACs) # create list of DummyPlug MotifObjs for compatibility
    
    if opts.verbose: sys.stdout.write('getting forgroundSeqs...\n')
    foregroundSeqs = getForeground(opts.genes)
    outData        = {}
    
    if opts.verbose: sys.stdout.write('calculating real data p-values...\n')
    realData = motifHyprGeoEnrichment(motifList,motifHits,foregroundSeqs,opts.expect)
    appendData(realData,outData)
    
    
    for i in range(opts.fdrs):
        if opts.verbose: sys.stdout.write('ctrl_%s:\n' % (i))
        if opts.verbose: sys.stdout.write('\tgetting random forground...\n')
        rForground = getRandomForground(foregroundSeqs,seqDict)
        
        if opts.verbose: sys.stdout.write('\tcalculating p-values...\n')
        ctrlData   = motifHyprGeoEnrichment(motifList,motifHits,rForground,opts.expect)
        appendData(ctrlData,outData)
    
    if opts.verbose: sys.stdout.write('writting outData...\n')
    writeDataTable(outBaseStr,outData,motifList)
    
    
    if opts.verbose: sys.stdout.write('plotting histograms...\n')
    for m in motifList:
        plotHist(outBaseStr,outData,m.id)
        
        
