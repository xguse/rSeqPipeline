from rSeq.utils.motifDiscovery.scope import *
from rSeq.utils.errors import *
import random
from time import time
import optparse
import sys

if __name__ == "__main__":

    
    #+++++++++++ File Parseing Etc +++++++++++
    desc = """Calls the folowing funcs: 'runSCOPE(pLen,genes,jobName,scopeDir,outDir,paramName,jMem='2000',verbose=False)'"""
    
    usage = """python %prog [args]"""
    parser = optparse.OptionParser(usage=usage, description=desc)
    
    parser.add_option('--plen', dest="len", type='int',default=1000,
                      help="""Promoter length to use (default=%default).""")
    parser.add_option('--genes', dest="genes", type='string',default=None,
                      help="""Unbroken string of gene/Tx names, sep=','. Exp: 'gene,gene,gene' (default=%default).""")
    parser.add_option('--job', dest="job", type='string',default='int(time())',
                      help="""String to identify this run (default=%default).""")
    parser.add_option('--scope', dest="scope", type='string',default=None,
                      help="""Full path to SCOPE's home dir (default=%default).""")
    parser.add_option('--out', dest="out", type='string',default=None,
                      help="""Path to results dir (default=%default).""")
    parser.add_option('--param', dest="param", type='string',default=None,
                      help="""Name of parameter file to use. (default=%default).""")
    parser.add_option('--list-params', dest="list", action='store_true',default=False,
                      help="""List names of param files for 'scopeDir'. (default=%default).""")
    parser.add_option('--mem', dest="mem", type='int',default=2000,
                      help="""Int representing how many mega-bytes of mem to give to java. (default=%default).""")
    parser.add_option('--verbose', dest="verb", action='store_true',default=False,
                      help="""Include if stdOut/stdErr is desired. (default=%default).""")
    parser.add_option('--random', dest="random", type='int',default=False,
                      help="""How many random sets of genes equal in length to --genes to run for FDR estimation. (default=%default).""")
    

    
    (opts, args) = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    if opts.list:
        if not opts.scope:
            raise InvalidOptionError('runSCOPE.py: --list-params requires --scope.')
        listParams(opts.scope)
        exit(0)

        
        
    # runSCOPE(pLen,genes,jobName,scopeDir,outDir,paramName,jMem='2000',verbose=False) 
    print "setting up scope run (real data)..."
    opts.genes = opts.genes.replace(',',';')
    realRes = runSCOPE(opts.len,
                       opts.genes,
                       opts.job,
                       opts.scope,
                       opts.out,
                       opts.param,
                       opts.mem,
                       opts.verb)
    if opts.random:
        # get correct number of random gene sets and run/save results.
        ctrlResList = []
        print "setting up scope run (random data)..."
        for i in range(int(opts.random)):
            randGenes = getRandomGeneSet(opts.scope,opts.param,len(opts.genes.split(';')))
            randGenes = ';'.join(randGenes)
            ctrlRes = runSCOPE(opts.len,
                               '%s' % (randGenes),
                               '%s.ctrl%s' % (opts.job,i),
                               opts.scope,
                               opts.out,
                               opts.param,
                               opts.mem,
                               opts.verb)
            ctrlResList.append(ctrlRes)
    