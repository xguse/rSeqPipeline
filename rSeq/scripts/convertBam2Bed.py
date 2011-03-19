# modified from http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_protocol#Python_APIs_.28Pysam.29
import os, sys, re, optparse
import pysam

def main( argv = None ):
    if not argv: argv = sys.argv

    # setup command line parser
    usage = "python %prog inFile.bam > outFile.bed"
    parser = optparse.OptionParser( version = "%prog version: $Id$", 
                                    usage = usage )
    parser.add_option("-r", "--region", dest="region", type="string",
                      help="samtools region string [default=%default]."  )
    parser.set_defaults( region = None, )
    (options, args) = parser.parse_args( argv[1:] )
    if len(args) != 1:
        parser.print_help()
        exit(0)

    # open BAM and get the iterator
    samfile = pysam.Samfile( args[0], "rb" )
    if options.region != None:
        it = samfile.fetch( region=options.region )
    else:
        it = samfile.fetch()

    # calculate the end position and print out BED
    take = (0, 2, 3) # CIGAR operation (M/match, D/del, N/ref_skip)
    outfile = sys.stdout
    for read in it:
        if read.is_unmapped: continue
        # compute total length on reference
        t = sum([ l for op,l in read.cigar if op in take ])
        if read.is_reverse: strand = "-"
        else: strand = "+"
        outfile.write("%s\t%d\t%d\t%s\t%d\t%c\n" %\
                          ( samfile.getrname( read.rname ),
                            read.pos, read.pos+t, read.qname,
                            read.mapq, strand) )            

if __name__ == "__main__":
    sys.exit( main( sys.argv) )