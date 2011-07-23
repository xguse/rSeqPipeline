import sys
import optparse
import csv
import scipy.stats as stats
import matplotlib as mpl
#mpl.use('TkAgg')
from matplotlib import pylab as pl
from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp

def main():
    """Inputs:
    -- DEGseq output table file
    -- bed file of DEGseq input A
    -- bed file of DEGseq input B
    -- plotting options
    -- output directory
    Outputs:
    -- pdf plot of normalized A reads vs B reads and the Pearson stats."""
    
    desc  = """This script creates a scatterplot of normalized transcript
read-counts from column 2 on the x-axis and column 3 on the right axis and reports the Pearson Correlation.  
It uses bedFileA and bedFileB to determine the correct normalization factor.  
It writes a pdf file of the plot to the supplied output directory.
bedFileA is bedfile whose data is in DEGseqFile column 2.
bedFileB is bedfile whose data is in DEGseqFile column 3."""
    
    usage = """python %prog [options] degSeqFile bedFileA bedFileB"""
    parser = optparse.OptionParser(usage=usage, description=desc)
    
    parser.add_option('--name-a',dest="name_a", type="string", default="Column 2",
                      help="""Name of data in column 2 of DEGseq file. (default=%default)""")
    parser.add_option('--name-b',dest="name_b", type="string", default="Column 3",
                      help="""Quoted string to use as name of data in column 3 of DEGseq file. (default=%default)""")
    parser.add_option('--dir',dest="dir", type="string", default='.',
                      help="""Directory to write plot to. (default=%default)""")
    parser.add_option('--norm',dest="norm",type='string',default=False,
                      help="""Provide a normalizer set and skip the bed counting step.  Example: "0.46,1" or "1,0.36". (default=%default)""")
    parser.add_option('--log',dest="log",action="store_true",default=False,
                      help="""Use log scale. (default=%default)""")
    parser.add_option('--show',dest="show",action="store_true",default=False,
                      help="""Show plot(s) in window. (default=%default)""")


    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    if opts.norm:
        try:
            n = map(float,tuple(opts.norm.split(',')))
            assert 1 in n, "** ERROR: At least one number in --norm must be '1'."
            opts.norm = n
        except:
            raise InvalidOptionError("Error in --norm option.  Please review.")
    if (len(args) != 3) and not (opts.norm):
        parser.print_help()
        print "\n\n** ERROR: Please supply exactly one of each: degSeqFile, bedFileA, bedFileB. **"
        exit(1)
    if not ((len(args) >= 1) and (opts.norm)):
        parser.print_help()
        print "\n\n** ERROR: Please supply degSeqFile. **"
        exit(1)
    
    
    opts.dir = opts.dir.rstrip('/')+'/'
    
    degPath = args[0]
    bedA    = args[1]
    bedB    = args[2]
    
    txCnts       = getTxCounts(degPath)
    
    if not opts.norm:
        nFactor = getNormFactor(bedA,bedB)
        print "NormFactor: %s" % (str(nFactor))
    else:
        nFactor = opts.norm
    
    normdTxCnts  = normTxCounts(txCnts,nFactor)
    pearsonStats = stats.pearsonr(normdTxCnts[0],normdTxCnts[1])
    plotScatter(pearsonStats,normdTxCnts,opts)
    print "Good bye!"
    
def getTxCounts(DEGseqPath):
    """Take DEGseq file path, return list of vectors:
    [[TxCntsA],[TxCntsB]]"""
    
    reader = csv.reader(open(DEGseqPath),delimiter='\t')
    header = reader.next() # we dont actually want this
    
    txCnts = [[],[]]
    for line in reader:
        try:
            txCnts[0].append(int(line[1]))
        except ValueError:
            txCnts[0].append(0)
        try:
            txCnts[1].append(int(line[2]))
        except ValueError:
            txCnts[1].append(0)

    return txCnts
    
    
def getNormFactor(bedApath,bedBpath):
    """Takes paths to bed_A and bed_B. Returns a tuple of fractions
    to be applied to TxCntsA and TxCntsB respectively in order to
    normalize the set with more reads down to the level of the smaller
    set."""
    
    bedACnt = 0
    bedBCnt = 0
    
    for line in open(bedApath,'rU'):
        if not line.startswith('track name='):
            bedACnt += 1
    for line in open(bedBpath,'rU'):
        if not line.startswith('track name='):
            bedBCnt += 1
    
    if float(bedACnt)/bedBCnt > 1:
        return (float(bedBCnt)/bedACnt,1)
    else:
        return (1,float(bedACnt)/bedBCnt)
    
def normTxCounts(txCntsList,normTup):
    """Takes output from getTxCounts and getNormFactor.
    Applies normilization, and returns normalized [[TxCntsA],[TxCntsB]]"""
    
    for i in range(len(txCntsList[0])):
        txCntsList[0][i] = txCntsList[0][i] * normTup[0]
        txCntsList[1][i] = txCntsList[1][i] * normTup[1]
    
    return txCntsList
    
def plotScatter(pearsonStats,normedTxCntsList,opts):
    """"""
    fig = pl.figure()
    ax  = fig.add_subplot(111)
    if opts.log:
        ax.set_xscale('log')
        ax.set_yscale('log')
    
    
    ax.scatter(normedTxCntsList[0],normedTxCntsList[1], s=15, c='b', marker='o', alpha=1)
    if not opts.log:
        ax.set_autoscale_on(False)
    ax.set_xlabel(opts.name_a)
    ax.set_ylabel(opts.name_b)
    upperLim = max(normedTxCntsList[0]+normedTxCntsList[1])
    

    
    m,b  = pl.polyfit(normedTxCntsList[0],normedTxCntsList[1],1)
    bfYs = pl.polyval([m,b], [1,max(normedTxCntsList[0])])
    
    ax.plot([1,max(normedTxCntsList[0])],bfYs,'r-')
    
    pl.text(0.01,0.99,'Pearson: %.4f, %s\nBest Fit: y=%.3f*x+%.3f' % (pearsonStats[0],pearsonStats[1],m,b),
            bbox=dict(facecolor='#87AACD', alpha=1),
            horizontalalignment='left',
            verticalalignment='top',
            transform = ax.transAxes)
    
    mkdirp(opts.dir)
    pl.savefig('%s%s_vs_%s.png' % (opts.dir,opts.name_a,opts.name_b))
    print 'Show?  %s' % (opts.show)
    if opts.show:
        pl.show()

        
if __name__ == "__main__":
    main()