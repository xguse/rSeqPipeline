"""Plot a transcript/gene's expression profile alone or in small groups using
table headers to coordinate data vs confidence/errs.  Extracts provided Tx/gene
data rows from full expression files.

Tx_symbol	gene_id	gene_short_name	status	2hr	2h_conf_lo	2hr_conf_hi	4hr	4hr_conf_lo	4hr_conf_hi
AGAP004678-RA	XLOC_000001	AGAP004678	OK	82.7505	28.0147	137.486	54.9937	19.2801	90.7074
AGAP004679-RB	XLOC_000002	AGAP004679	OK	8.98055	2.94927	15.0118	16.4392	5.58947	27.2889
AGAP004679-RA	XLOC_000002	AGAP004679	OK	26.4958	9.23621	43.7554	16.0875	5.47647	26.6985

Uses headers to determine which data to plot.

"""

import sys
import os
import optparse

import numpy as np
from matplotlib import rc
from matplotlib import pylab as pl
rc('text', usetex=True)

from rSeq.utils.externals import mkdirp,runExternalApp
from rSeq.utils.files import tableFile2namedTuple
from rSeq.utils.misc import Bag
from rSeq.utils.plotting import setTickSizes


def main():
    desc = "Plot a transcript/gene's expression profile alone or in small groups using table headers to coordinate data vs confidence/errs.  Extracts provided Tx/gene data rows from full expression files."
    usage = """python %prog inFile [options] """
    parser = optparse.OptionParser(usage,epilog=desc)
    parser.add_option('-i',dest="include",type='string',default=False,
                      help="""Quoted comma delim'd list of tx/gene symbols to plot. (default=%default)""")
    parser.add_option('--names-col',dest="names_col",type='string',default=False,
                      help="""Header name where tx/gene names live. (default=%default)""")
    parser.add_option('--dir',dest="dir",type='string',default='.',
                      help="""Directory to serve as the root of the gallery. (default=%default)""")
    parser.add_option('-p',dest="prefix",type='string',default=False,
                      help="""Save figs using file-name prefix provided. (default=inFile)""")
    parser.add_option('--conditions',dest="conditions",type="string", default=False,
                      help="""Quoted, ordered, comma delim'd list of headers where condition/time-point expression values live. (default=%default)""")
    parser.add_option('--sym-bounds',dest="sym_bounds",type="string", default=False,
                      help="""Quoted, ordered (matching --conditions), comma delim'd list of headers where symetrical confidence/err values live. (default=%default)""")
    parser.add_option('--hi-bounds',dest="hi_bounds",type="string", default=False,
                      help="""Quoted, ordered (matching --conditions), comma delim'd list of headers where upper confidence/err values live. (default=%default)""")
    parser.add_option('--lo-bounds',dest="lo_bounds",type="string", default=False,
                      help="""Quoted, ordered (matching --conditions), comma delim'd list of headers where lower confidence/err values live. (default=%default)""")
    parser.add_option('--xlabel',dest="xlabel",type="string", default=False,
                      help="""Quoted string [TeX allowed]: "the xlabel text". (default=%default)""")
    parser.add_option('--ylabel',dest="ylabel",type="string", default=False,
                      help="""Quoted string [TeX allowed]: "the ylabel text". (default=%default)""")
    parser.add_option('--xticks',dest="xticks",type="string", default=False,
                      help="""Quoted comma delim'd string for x-axis tick text [TeX allowed]: "4 hrs,6 hrs,...". (default=%default)""")
    parser.add_option('--ylog',dest="ylog",action='store_true', default=False,
                      help="""Plot as y-axis as log. (default=%default)""")
    parser.add_option('--overlay',dest="overlay",action='store_true', default=False,
                      help="""Plot all in same file. (default=%default)""")
    parser.add_option('--color',dest="color",type='string', default=None,
                  help="""Manually set marker color (html codes or names are legal) (default=%default)""")
    parser.add_option('--no-show', dest="no_show",action="store_true", default=False,
                      help="""Suppress figure displays. (default=%default)""")
    parser.add_option('--galaxy', dest="galaxy",action="store_true", default=False,
                      help="""Used to make suitable for use with Galaxy tool template. (default: %default)""")
    figTypes = ['pdf','png']
    parser.add_option('--type',dest="type",type='string',default='pdf',
                      help="""Save figs as: """ + str(figTypes) + """. (default=%default)""")
    
    (opts, args) = parser.parse_args()
    print opts
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    
    
    if not opts.type in figTypes:
        raise ValueError(" Option --types must be one of %s, given: %s " % (figTypes,opts.type))
    
    if len(args) != 1:
        parser.print_help()
        raise ValueError(" Please supply exactly one inFile. ")
    if not opts.prefix:
        opts.prefix = args[0]
        
    if opts.overlay != False:
        raise NotImplementedError(" --overlay has not been implemented yet (contact wadunn83@gmail.com). ")
    
    if opts.names_col == False:
        raise ValueError(" You must provide the header name where tx/gene names live. ")
    if opts.include == False:
        raise ValueError(" You must provide at least ONE tx/gene symbol. ")
    else:
        opts.include = opts.include.split(',')
    
    if opts.conditions == False:
        raise ValueError(" You must provide at least ONE condition. ")
    else:
        opts.conditions = opts.conditions.split(',')
    
    # Option translation kluge here to deal with galaxy
    if opts.sym_bounds == '':
        opts.sym_bounds = False
    if opts.hi_bounds == '':
        opts.hi_bounds = False
    if opts.lo_bounds == '':
        opts.lo_bounds = False
        
    if opts.sym_bounds != False:
        opts.sym_bounds = opts.sym_bounds.split(',')
    if ((opts.sym_bounds != False) and ((opts.hi_bounds != False) or (opts.lo_bounds != False))):
        raise ValueError(" If --sym-bounds: NEITHER --hi-bounds or --lo-bounds should be used. ")
    if ((opts.hi_bounds != False) and (opts.lo_bounds == False) or (opts.hi_bounds == False) and (opts.lo_bounds != False)):
        raise ValueError(" BOTH or NEITHER --hi-bounds and --lo-bounds must be used. ")
    elif ((opts.sym_bounds == False) and (opts.hi_bounds == False) and (opts.lo_bounds == False)):
        pass
    else:
        opts.hi_bounds = opts.hi_bounds.split(',')
        opts.lo_bounds = opts.lo_bounds.split(',')
    
    if opts.xticks:
        opts.xticks = opts.xticks.split(',')
    
    
    
      
    # --- Parse inFile ---
    table = tableFile2namedTuple(tablePath=args[0],sep='\t',headers=None)
    data  = []
    for row in table:
        if row.__getattribute__(opts.names_col) in opts.include:
            data.append(row)
            
    # --- Create out directory ---
    if not opts.galaxy:
        mkdirp(opts.dir)
    else:
        opts.prefix = opts.prefix.rstrip('/').split('/')[-1]
        opts.dir = opts.dir.rstrip('/') + '_files'
        mkdirp('%s/' % (opts.dir))
        
        
    
    # ------- Build Figs -------
    for gene in data:
        fig = pl.figure()
        ax  = fig.add_subplot(111)
        ax.set_title(r"{\Huge %s}" % (gene.__getattribute__(opts.names_col)))
        if opts.xlabel:
            ax.set_xlabel(r"{\Large %s}" % (opts.xlabel))
        if opts.ylabel:
            ax.set_ylabel(r"{\Large %s}" % (opts.ylabel))
        
        x    = range(len(opts.conditions))
        y    = np.array([float(gene.__getattribute__(i)) for i in opts.conditions])
        
        if opts.sym_bounds != False:
            yerr = np.array([float(gene.__getattribute__(i)) for i in opts.sym_bounds])
            yerr = abs(y-yerr)
        elif opts.hi_bounds != False:
            yHi  = np.array([float(gene.__getattribute__(i)) for i in opts.hi_bounds])
            yLo  = np.array([float(gene.__getattribute__(i)) for i in opts.lo_bounds])
            yerr = np.array([abs(y-yLo),abs(y-yHi)])
        else:
            yerr = None
            

        pl.errorbar(x, y, yerr=yerr, xerr=None, fmt='s', ecolor='0.55',
                    elinewidth=None, capsize=7, barsabove=False, lolims=False,
                    uplims=False, xlolims=False, xuplims=False, lw=3, ms=10, color=opts.color)
        ax.set_xlim(-1,len(opts.conditions))
        
        if opts.xticks:
            locs, labels = pl.xticks(x, [i.replace('_','\_') for i in opts.xticks])
        else:
            locs, labels = pl.xticks(x, [i.replace('_','\_') for i in opts.conditions])
        #pl.setp(labels, 'rotation', 'vertical')

        
        setTickSizes(axObj=ax,fontSize=16)
        pl.savefig('%s/%s.%s.expPlot.%s' % (opts.dir.rstrip('/'),opts.prefix,gene.__getattribute__(opts.names_col),opts.type))

        if not opts.no_show:
            pl.show()
        del(fig,ax)
        
    # ---- Build HTML Gallery ----
    curResult = runExternalApp('curator','%s' % (opts.dir))
    
    # ---- Create a link file to the sortindex.html file ----
    lnSrc  = '%s/sortindex.html' % (opts.dir)
    lnDest = '%s/%s' % ('/'.join(opts.dir.split('/')[:-1]),opts.prefix)
    os.symlink(lnSrc,lnDest)

        

    
if __name__ == "__main__":
    main()