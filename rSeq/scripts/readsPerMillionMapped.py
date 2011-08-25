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
    -- table of transcript/gene symbols, condition/sample name,original read count, and read counts normalized to each million mapped reads."""
    
    
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
                      help="""Quoted, comma-delim'd list of transcript/gene symbols. (default=%default)""")


    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)
    
    if not opts.l == "All":
        opts.l = opts.l.split(',')
    else:
        opts.l = ''

    opts.tx_col = int(opts.tx_col)
    opts.rd_col = int(opts.rd_col)
    
    # get million mapped reads
    milMapped = float(pysam.flagstat(opts.b)[3].split()[0])/1000000
    #milMapped = 7.5 # for quick debugging to bybass the bam filtering
    
    # parse list of tx/gene symbols into dict for tracking successes
    ##symbols = {}
    ##if not opts.l == "All":
        ##for sym in opts.l:
            ##symbols[sym] = []
    
        
    
    # open expression table into rows
    rows = csv.reader(open(opts.t),delimiter='\t')
    
    # for each tx/gene in expFile: output million mapped reads (MMR)
    #   if row[opts.tx_col].startswith(<any of the requested symbols>)
    #   and update symbols with True if a hit is found.
    print "Tx_symbol\tCondition\tOriginal_reads\tReads_per_million_mapped"
    for row in rows:
        if row[opts.tx_col].startswith(tuple(opts.l)):
            try:
                print "%s\t%s\t%s\t%s" % (row[opts.tx_col],
                                      opts.c,
                                      row[opts.rd_col],
                                      float(row[opts.rd_col])/milMapped)
            except:
                print "failed:  %s" % (';'.join(row))
    
    # report whether/which requested symbols were not found or had errors??


if __name__ == "__main__":
    main()