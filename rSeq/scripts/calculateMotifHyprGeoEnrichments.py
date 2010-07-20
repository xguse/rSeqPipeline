import sys
import os
import optparse
import csv
import random as r
import math
from cogent import LoadSeqs
from matplotlib import pylab as pl
from rSeq.utils.motifDiscovery.motifs import getHitDict,parseSCOPEfile,parseXMSfile,motifHyprGeoEnrichment
from rSeq.utils.stats import seqStats
from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp


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
    dataWriter = csv.writer(outFile, delimiter='\t')
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
        return int(math.sqrt(len(ctrls)))*2 

        
    outFile = '%s_%s_hist.pdf' % (outBaseStr,motifID)
    ctrls   = [-math.log10(x) for x in outData[motifID][1:]]
    yLab    = 'Enrichment [-log10(pVal)]'
    xLab    = 'Controls [n=%s]' % (len(ctrls))
    
    fig  = pl.figure()
    ax = fig.add_subplot(111)
    
    bins = _chooseBins()
    hData = ax.hist(ctrls, bins=bins, range=None,
                    normed=False, weights=None,
                    cumulative=False, bottom=None,
                    histtype='bar', align='mid',
                    orientation='vertical',
                    rwidth=None, log=False, hold=None)
    ax.axvline(x=-math.log10(outData[motifID][0]) ,linewidth=2, color='r')
    pl.savefig(outFile)

    
    

if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    desc = """Calls the folowing funcs: 'add this'"""
    
    usage = """python %prog [add this]"""
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
    parser.add_option('--verbose',action='store_true',default=False,
                      help="""Include if stdOut/stdErr is desired. (default=%default).""")
    parser.add_option('--plot-fdrs', dest="fdrs", type='int',default=0,
                      help="""How many random sets of genes equal in length to --genes to run for FDR estimation. (default=%default).""")
    

    
    (opts, args) = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    outBaseStr = getOutBaseStr(opts.out,opts.job)
    if opts.verbose:
        sys.stdout.write('\n%s\n\n' % (' '.join(sys.argv)))
    
        
    mkdirp(opts.out)
    if opts.verbose: sys.stdout.write('building seqDict...\n')
    seqDict        = getSeqs(opts.promoters,opts.plen)
    
    if opts.verbose: sys.stdout.write('building motifList...\n')
    motifList      = getMotifs(opts.motifs,opts.motif_type)
    
    if opts.verbose: sys.stdout.write('getting nucFreqs...\n')
    halfAT,halfGC  = getSeqFreqs(seqDict)
    
    if opts.verbose: sys.stdout.write('building hitDict...\n')
    motifHits      = getHitDict(motifList,seqDict,pThresh=opts.thresh,halfAT=halfAT,halfGC=halfGC)
    
    if opts.verbose: sys.stdout.write('getting forgroundSeqs...\n')
    foregroundSeqs = getForeground(opts.genes)
    outData        = {}
    
    if opts.verbose: sys.stdout.write('calculating real data p-values...\n')
    realData = motifHyprGeoEnrichment(motifList,hitDict,foregroundSeqs)
    appendData(realData,outData)
    
    
    for i in range(opts.fdrs):
        if opts.verbose: sys.stdout.write('ctrl_%s:\n' % (i))
        if opts.verbose: sys.stdout.write('\tgetting random forground...')
        rForground = getRandomForground(foregroundSeqs,seqDict)
        
        if opts.verbose: sys.stdout.write('\tcalculating p-values...\n')
        ctrlData   = motifHyprGeoEnrichment(motifList,hitDict,rForground)
        appendData(ctrlData,outData)
        
        if opts.verbose: sys.stdout.write('\tplotting histograms...\n')
        plotHist(outBaseStr,outData,motifID)
