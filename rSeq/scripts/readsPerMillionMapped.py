import sys
import optparse
import csv

import pysam
import scipy.stats as stats

from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp

def main():
    """Inputs:
    -- Expression output table file
    -- Transcript/gene symbol containing column number
    -- Reads per Tx/gene containing column number
    -- Bam file
    -- List of transcript/gene symbols
    -- Name of condition/sample
    Outputs:
    -- table of transcript/gene symbols, condition/sample name, and read counts normalized to each million mapped reads."""
    
    
    usage = """python %prog [options]"""
    parser = optparse.OptionParser(usage=usage)
    
    parser.add_option('-t', type="string", default=False,
                      help="""Name of a table file containing Tx/Gene symbols and respective read counts for at least one condition/sample. (default=%default)""")
    parser.add_option('-c', type='string', default=False,
                      help="""Name of the condition/sample. (default=%default)""")
    parser.add_option('--tx-col', dest="tx_col", type='string', default=False,
                      help="""The column number containing Tx/Gene symbols. (default=%default)""")
    parser.add_option('--rd-col', dest="rd_col", type='string', default=False,
                      help="""The column number containing read counts. (default=%default)""")
    parser.add_option('-b', type='string', default=False,
                      help="""Path to the bam file representing the aligned reads for the desired condition/sample. (default=%default)""")
    parser.add_option('-l', type='string', default="All",
                      help="""List of transcript/gene symbols. (default=%default)""")



    
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    
    
    opts.dir = opts.dir.rstrip('/')+'/'
    
    # get million mapped reads
    milMapped = float(pysam.flagstat(opts.b)[3].split()[0])/1000000
    
    # parse list of tx/gene symbols into dict for tracking successes
    symbols = {}
    
    # open expression table into rows
    rows = csv.reader(open(opts.t),delimiter='\t')
    
    # for each tx/gene in expFile: output million mapped reads (MMR)
    #   if row[opts.tx-col].startswith(<any of the requested symbols>)
    #   and update symbols with True if a hit is found.
    for row in rows:
        pass
    
    # report whether/which requested symbols were not found or had errors??


if __name__ == "__main__":
    main()