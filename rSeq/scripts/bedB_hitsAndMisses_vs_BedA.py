from rSeq.utils.bed import divByWindow
from rSeq.utils.errors import *
import optparse
import sys

#bedA = '/home/dunnw/data/genomes/aaegypti.SUPERCONTIGS-Liverpool.AaegL1/aaegypti.GENES-AaegL1.2.bed'
#bedB = '/home/dunnw/data/solexa/bowtie_out/LSxLBx_bowtie.bed'
#outDir = '/home/dunnw/data/solexa/bowtie_out/LSxLBx_bowtie_results'




if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    epilog = """Calls the folowing funcs:
**divByWindow(bedA_Path,bedB_Path,win,cols,side,outDir)**:
Create files separating features in bedB by those alling within the area defined
by <win> and those outside this area in bedA.  If A.bed is stranded, the area is defined
by win[0] upstrm and win[1] dwnstrm on the FEATURE's strand.  Otherwise its 
win[0] upstrm and win[1] dwnstrm on the CONTIG/CHROM's plus strand.  Files ouput
to outDir.

**extractFromDoubleSidedBedtoolOut(filePath,cols,side,outDir)**:
Creates new file from filePath using only the bedInfo from the
left/right (based on 'side') side of a BEDtools outFile with double-
sided output (side=[3,6]). 'cols' must be a list with length of columns
in each 'side' of the double output.  'side' = keep the 'right' or
'left' side of the output line."""
    
    usage = """python %prog [options]  bedPathA bedPathB"""
    parser = optparse.OptionParser(usage=usage, epilog=epilog)
    
    parser.add_option('--out', dest="out", type='string',default='',
                      help="""Directory path for output (default=%default)""")
    parser.add_option('--win',dest="win",type="string", default="2000,500", 
                      help="""Comma separated string of left and right window sizes.  (default=%default)""")
    parser.add_option('--cols', dest="cols", type='string',default='6,6',
                      help="""Comma separated string of left and right window sizes **See "extractFromDoubleSidedBedtoolOut()" (default=%default)""")
    parser.add_option('--side', dest="side", type='string',default='right',
                      help="""**See "extractFromDoubleSidedBedtoolOut()". (default=%default)""")
    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    if len(opts.win.split(',')) != 2:
        raise InvalidOptionError('Malformed --win option (%s).' % (opts.win))
    if len(opts.cols.split(',')) != 2:
        raise InvalidOptionError('Malformed --cols option (%s).' % (opts.cols))
    if opts.side not in ['left','right']:
        raise InvalidOptionError('Illegal --side option (%s). Legals are %s.' % (opts.side,['left','right']))



print divByWindow(args[0],
                  args[1],
                  win=[int(x) for x in opts.win.split(',')],
                  cols=[int(x) for x in opts.cols.split(',')],
                  side=opts.side,
                  outDir=opts.out)
