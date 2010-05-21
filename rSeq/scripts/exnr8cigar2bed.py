from rSeq.utils.align import exonerateCigar2BEDline
from rSeq.utils.externals import mkdirp
from rSeq.utils.errors import *
import optparse
import sys




if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    epilog = """Calls the folowing funcs:
**exonerateCigar2BEDline(cigarLine,rgb="0,0,0"):
Takes a cigar line from an exonerate output file and returns
a BED formated line representing the alignment blocks of the query
as aligned to the target. """
    
    usage = """python %prog [options]  exonerateOut"""
    parser = optparse.OptionParser(usage=usage, epilog=epilog)
    
    parser.add_option('--out', dest="out", type='string',default='.',
                      help="""Directory path for output (default=%default).""")
    parser.add_option('--track', dest="track", type='string',default='untitled',
                      help="""Unbroken string to use for track name (default=%default).""")
    parser.add_option('--description', dest="description", type='string',default='none given',
                      help="""Quoted string for description of track (default=%default).""")
    parser.add_option('--rgb', dest="rgb", type='string',default='0,0,0',
                      help="""Unbroken comma separated string to denote color of track (default=%default).""")

    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    if opts.out.endswith('/'):
        opts.out = opts.out[:-1]
    if not len(opts.rgb.split(',')):
        raise InvalidOptionError('Malformed --rgb option: %s.' % (opts.rgb))
    
    # Prepare outdir 
    mkdirp(opts.out)
    outFile = open('%s/%s' % (opts.out,args[0].split('/')[-1].replace('.txt', '.bed')))
    outFile.write('track name=%s description="%s" useScore=0\n' % (opts.track,opts.description))
    
    # For every line starting "cigar:" convert to bed and write out
    for line in open(args[0],'rU'):
        if line.startswith('cigar:'):
            bedLine = exonerateCigar2BEDline(line,opts.rgb)
            outFile.write(bedLine)
    outFile.close()
