# modified from http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_protocol#Python_APIs_.28Pysam.29
import os, sys, re, optparse

from rSeq.utils.errors import *
from rSeq.utils.convert import bam2bed
from rSeq.utils.misc import RseqHelpFormatter

       

if __name__ == "__main__":
    # setup command line parser
    desc = """DESCRIPTION:
Opens bams with pysam to access the contents then converts each read's info to bed format and writes to beds.

BAMS:\tpath string or quoted list of comma-delimited path strings
BEDS:\tpath string or quoted list of comma-delimited path strings
region:\tsamtools region string or quoted list of comma-delimited region strings
ncpus:\tnum of processors to use if parallel job is run (autodetect -> all avail)

If multiple BAMS and one BEDS path, all bams written to same bed.

If multi BAMS AND BEDS, they must be of equal length, and BAM[i] is written to BED[i].

If multi REGIONS, they must be of same length as BAMs and BAM[i] will use REGION[i].

If single REGIONS, it will be used for all BAMs.

When avail, the jobs will be run in parallel.

\nThis is modified from http://sourceforge.net/apps/mediawiki/samtools/index.php?title=SAM_protocol#Python_APIs_.28Pysam.29
"""
    
    usage = "python %prog [options] --bams 'bam1,bam2..' --beds 'bed1,bed2..' -r 'region1,region2..'"
    parser = optparse.OptionParser(usage=usage, epilog=desc, formatter=RseqHelpFormatter())
    
    parser.add_option("--bams", dest="bams", type="string", default=None,
                      help="Path(s) to BAM files. REQUIRED [default=%default].")
    parser.add_option("--beds", dest="beds", type="string", default=None,
                      help="Path(s) to BEDS[i].bed.  If None: uses sys.stdout. [default=%default].")
    parser.add_option("--ncpus", dest="ncpus", type="string", default='autodetect',
                      help="Number of processors to use. autodetect uses all avail. [default=%default].")
    parser.add_option("-r", dest="regions", type="string", default=None,
                      help="samtools region string [default=%default].")
    
    (opts, args) = parser.parse_args()
    
    
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        exit(0)
    
    # process bams
    if not opts.bams:
        raise MissingArgumentError("ERROR: BAMS option is required.")
    else:
        opts.bams = opts.bams.split(',')
    
    # process beds
    if opts.beds:
        opts.beds = opts.beds.split(',')
        
    # process regions
    if opts.regions:
        opts.regions = opts.regions.split(',')
    
    # Validate ncpus
    if opts.ncpus != 'autodetect':
        try:
            opts.ncpus = int(opts.ncpus)
        except ValueError:
            raise InvalidOptionError('ERROR: NCPUS must be an integer.')
    
    # lets do this:
    bam2bed(bamPath=opts.bams,
            bedFile=opts.beds,
            region=opts.regions,
            ncpus=opts.ncpus)
    