# Run the requested operations to initialize standard indexes and files
# related to the common operations run on genomic data with a new set of
# genomic source files.

import sys
import os
import optparse

from iniparse import SafeConfigParser

from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp,runExternalApp
from rSeq.utils.align import bowtie_index
from rSeq.utils.sitRep import start_sitrep

def init_dir_structure(spcName,versionID,baseDir,isCurrent=False):
    """Ensure that the correct directory structure exists to accept
    the genome and index data files.  Create anything that does not already exist.
    If isCurrent != False: soft link the versionID dir as "current"
    
    Example args:
    spcName = 'aedes_aegypti'
    versionID   = 'release_7'
    baseDir     = '/home/data'
    """
    
    spcVerDir   = '%s/genomes/%s/%s' % (baseDir,spcName,versionID)
    fastas      = '%s/fasta/' % (spcVerDir)
    annotations = '%s/annotations/' % (spcVerDir)
    mysql       = '%s/mysql/' % (spcVerDir)
    indexes     = '%s/indexes/' % (baseDir)
    
    print "creating dir: %s" % (fastas)
    mkdirp(fastas)
    print "creating dir: %s" % (annotations)
    mkdirp(annotations)
    print "creating dir: %s" % (mysql)
    mkdirp(mysql)
    print "creating dir: %s" % (indexes)
    mkdirp(indexes)
    
    
    if isCurrent:
        print "creating sym link from %s to 'current' because <isCurrent> is set to True" % (spcVerDir)
        os.symlink(spcVerDir, '%s/genomes/%s/current' % (baseDir,spcName))
    
    

def blast_formatdb():
    """"""
    
def build_esa():
    """"""

def init_sqlDB():
    """"""


if __name__ == "__main__":
    ##os.environ['BOWTIE_INDEXES'] ="/home/augustine/work/data/genomes/bowtie_indexes"
    # Default values
    defaultXXXXOpts = '--opts x'
    
    #+++++++++++ File Parseing Etc +++++++++++
    epilog = """DESCRIPTION: xxxx"""
    
    usage = """python %prog xxxxx"""
    parser = optparse.OptionParser(usage=usage, epilog=epilog)
    
    parser.add_option('--xxx', dest="xxx", type='string',default=None,
                      help="""xxxx. (default=%default)""")
    parser.add_option('--yyy', dest="yyy", action='store_true',default=False,
                      help="""yyyy. (default=%default)""")
    parser.add_option('--cent-log', dest="cent_log", action='store_true',default=False,
                      help="""Include to redirect stdout and stderr to a central log file in out_dir. [useful for reproducibility and debugging] (default=%default)""")
    
    (opts, args) = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    if opts.out_dir:
        mkdirp(opts.out_dir)
        opts.out_dir = opts.out_dir.rstrip('/')
    else:
        opts.out_dir = os.getcwd()
    if opts.cent_log:
        start_sitrep(opts.out_dir)
    if not 'BOWTIE_INDEXES' in os.environ:
        raise Exception('ERROR: please set the BOWTIE_INDEXES environment variable!')
    try:
        opts.bt_index
        opts.bam_base
        opts.fastqs
    except KeyError:
        raise MissingArgumentError('missing at least one of the required command line arguments: --bt-index, --bam-base, --fastqs')        
    if not opts.override:
        opts.bt_opts = "%s %s" % (defaultBtOpts,opts.bt_opts)
    

    
