# modified from http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_protocol#Python_APIs_.28Pysam.29
import os, sys, re, optparse

import pysam

from rSeq.utils.convert import bam2bed

       

if __name__ == "__main__":
    # setup command line parser
    epilog = """DESCRIPTION: Opens each inFile.bam with pysam to access its contents then converts each read's info to bed format and appends to OUT_FILE."""
    usage = "python %prog [options] inFile1.bam [inFile2.bam ..]"
    parser = optparse.OptionParser(usage=usage, epilog=epilog)
    parser.add_option("-r", dest="region", type="string", default=None,
                      help="samtools region string [default=%default].")
    parser.add_option("-o", dest="out_file", type="string", default=None,
                      help="Path to outFile.bed.  If None: uses sys.stdout. [default=%default].")

    (opts, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        exit(0)
    if opts.out_file == None:
        opts.out_file = sys.stdout
    else:
        opts.out_file = open(opts.out_file,'a')
    
    # lets do this:
    for bamPath in args:
        bam2bed(bamPath = bamPath,
                bedFile = opts.out_file,
                region  = opts.region)