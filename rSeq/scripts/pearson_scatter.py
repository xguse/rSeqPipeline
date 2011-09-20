import sys,os
import argparse
import csv

import scipy.stats as stats
import matplotlib as mpl
from matplotlib import pylab as pl

from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp
from rSeq.utils.files import tableFile2namedTuple

def main():
    """Inputs:
    -- Data tables containing at least 1 column of data tied to a column of data symbols
    -- Column names/positions to ID data symbols and data
    -- Plotting options
    -- Outfile path
    Outputs:
    -- Image file containing the scatterplot, Pearson stats, other useful info."""
    
    desc  = """This script creates a scatterplot of FPKM transcript read-counts and reports the Pearson Correlation."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('table1', type=str,
                        help="""Path to first data-table.""")
    parser.add_argument('data1', type=str,
                        help="""Column header name or 0-based column number where data for table1 lives.""")
    parser.add_argument('ids1', type=str,
                        help="""Column header name or 0-based column number where the data symbols for table1 lives.""")
    parser.add_argument('--label1', type=str, default=' ',
                        help="""Axes name for data1. (Default:%(default)s)""")
    parser.add_argument('table2', type=str,
                        help="""Path to second data-table.""")
    parser.add_argument('data2', type=str,
                        help="""Column header name or 0-based column number where data for table2 lives.""")
    parser.add_argument('ids2', type=str,
                        help="""Column header name or 0-based column number where the data symbols for table2 lives.""")
    parser.add_argument('--label2', type=str, default=' ',
                        help="""Axes name for data2. (Default:%(default)s)""")
    parser.add_argument('--log', action='store_true',
                        help="""Plot the points on a log:log scale. (Default: %(default)s)""")
    parser.add_argument('--show', action='store_true',
                        help="""Plot the image for interactive manipulation, otherwise just write the file. (Default: %(default)s)""")
    parser.add_argument('--pdf', action='store_true',
                        help="""Plot the image as a pdf: png otherwise. Png is preferable when data size is large. (Default: %(default)s)""")
    parser.add_argument('--galaxy', action='store_true',
                        help="""Use symplified namings suitable for use with Galaxy tool template. (Default: %(default)s)""")
    parser.add_argument('--out', type=str, default='',
                        help="""Base path for output [superseded by --galaxy]. (Default: current working directory)""")




    # print the called command:
    print " ".join(sys.argv)
    
    args = parser.parse_args()
    
    
    
    try:
        args.table1 = tableFile2namedTuple(args.table1)
        args.table2 = tableFile2namedTuple(args.table2)
    except:
        raise
    
    # deal with names/header posibilities
    
    data = collectData(args) 
    pearsonStats = stats.pearsonr(data[0],data[1])
    plotScatter(pearsonStats,data,args)
    print "Good bye!"
    
def collectData(args):
    """
    Culls data from the data-tables, associates it with the data-symbols,
    links the two sets and outputs [[data1],[data2]] for plotting and pearson.
    """
    tmp = {}
    for row in args.table1:
        symbol = row.__getattribute__(args.ids1)
        datum  = row.__getattribute__(args.data1)
        try:
            tmp[symbol][0].append(datum)
        except KeyError:
            tmp[symbol]    = [[],[]]
            tmp[symbol][0] = datum
    
    for row in args.table2:
        symbol = row.__getattribute__(args.ids2)
        datum  = row.__getattribute__(args.data2)
        try:
            tmp[symbol][1] = datum
        except KeyError:
            tmp[symbol]    = [[],[]]
            tmp[symbol][1] = datum
    
    data = [[],[]]
    for each in tmp:
        d = tmp[each]
        data[0].append(float(d[0]))
        data[1].append(float(d[1]))
    del(tmp)
    return data
        
def plotScatter(pearsonStats,data,args):
    """"""
    fig = pl.figure()
    ax  = fig.add_subplot(111)
    if args.log:
        ax.set_xscale('log')
        ax.set_yscale('log')
    
    
    ax.scatter(data[0],data[1], s=15, c='b', marker='o', alpha=1)
    if not args.log:
        ax.set_autoscale_on(False)
    ax.set_xlabel(args.label1)
    ax.set_ylabel(args.label2)
    upperLim = max(data[0]+data[1])
    

    
    m,b  = pl.polyfit(data[0],data[1],1)
    bfYs = pl.polyval([m,b], [1,max(data[0])])
    
    ax.plot([1,max(data[0])],bfYs,'r-')
    
    pl.text(0.01,0.99,'Pearson: %.4f, %s\nBest Fit: y=%.3f*x+%.3f' % (pearsonStats[0],pearsonStats[1],m,b),
            bbox=dict(facecolor='#87AACD', alpha=1),
            horizontalalignment='left',
            verticalalignment='top',
            transform = ax.transAxes)
    
    if args.pdf:
        pdf_or_png = 'pdf'
    else:
        pdf_or_png = 'png'
        
    # construct outfile name
    if not args.log:
        outName = '%s_%s_vs_%s.%s' % (args.out,args.label1,args.label2,pdf_or_png)
    else:
        outName = '%s_%s_vs_%s.log.%s' % (args.out,args.label1,args.label2,pdf_or_png)
        
    pl.savefig(outName)
    
    if args.galaxy:
        os.rename(outName,args.out)
        
    print 'Show?  %s' % (args.show)
    if args.show:
        pl.show()


        
if __name__ == "__main__":
    main()