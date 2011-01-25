import sys
import os
import optparse

import pysam

from rSeq.utils.errors import *
from rSeq.utils.externals import mkdirp,runExternalApp
from rSeq.utils.align import bowtie_align,samtools_sort,samtools_view,samtools_index
from rSeq.utils.sitRep import start_sitrep






if __name__ == "__main__":
    ##os.environ['BOWTIE_INDEXES'] ="/home/augustine/work/data/genomes/bowtie_indexes"
    # Default values
    defaultBtOpts = '--solexa1.3-quals -v 2 -m 1 -S'
    
    #+++++++++++ File Parseing Etc +++++++++++
    epilog = """DESCRIPTION: Runs bowtie, converts to BAM, sorts and indexes the BAM file to output directory.  Requires that BOWTIE_INDEXES environment variable be set to the home of your precompiled bowtie indexes.  Stared (*) arguments are required.  """
    
    usage = """python %prog --bt-index <idx_name> --bam-base <bam_name> --fastqs <fastq_files> --bt-opts <bowtie_options>"""
    parser = optparse.OptionParser(usage=usage, epilog=epilog)
    
    parser.add_option('--bt-index', dest="bt_index", type='string',default=None,
                      help="""*Base name of pre-compied bowtie index. (default=%default)""")
    parser.add_option('--bam-base',dest="bam_base",type="string", default=None, 
                      help="""*Base name of resulting alignment files. Do not include file extentions OR a full path. (default=%default)""")
    parser.add_option('--fastqs', dest="fastqs", type='string',default=None,
                      help="""*Appropriate quoted string representing which fastq files to use (see "bowtie -h"). (default=%default)""")
    parser.add_option('--bt-opts', dest="bt_opts", type='string',default=defaultBtOpts,
                      help="""Quoted string to pass as arguments to bowtie that have not already been provided.  Defualt mimics Eland. Unless --override is used, the options will be appended to the default. (default=%default)""")
    parser.add_option('--out-dir', dest="out_dir", type='string',default=None,
                      help="""Central directory to deposit output files. (default=%default)""")
    parser.add_option('--cent-log', dest="cent_log", action='store_true',default=False,
                      help="""Include to redirect stdout and stderr to a central log file in out_dir. [useful for reproducibility and debugging] (default=%default)""")
    parser.add_option('--override', dest="override", action='store_true',default=False,
                      help="""A switch that causes the default --bt-opts to be replaced by those provided as arguments to --bt-opts instead of adding to them. (default=%default)""")
    
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
    
    # Change to out_dir to make sure that any tempfiles are not droped in random
    #   places in the file system in case of fatal errors
    os.chdir(opts.out_dir)
        
    # Construct the bowtie comand
    ##btResults = bowtie_align(ebwt=opts.bt_index,
                             ##readsString=opts.fastqs,
                             ##hit='%s.sam' % (opts.bam_base),
                             ##runDir=opts.out_dir,
                             ##options=opts.bt_opts)
    
    # +++ Convert SAM to BAM +++ #
    samViewIn  = '%s/%s.sam' % (opts.out_dir,opts.bam_base)
    samViewOut = '%s/%s.bam' % (opts.out_dir,opts.bam_base)
    
    samViewArgs = ['-b',         # output should be BAM
                   '-S',         # input is SAM
                   samViewIn,
                   '-o',         # write output to file vs stdout
                   samViewOut]
    
    #viewResults = pysam.view(*samViewArgs)
    viewResults = samtools_view(samViewArgs)
    
    
    # +++ Sort BAM +++ #
    samSortIn   = samViewOut
    samSortOut  = '%s/%s.srt' % (opts.out_dir,opts.bam_base)
    samSortArgs = [samSortIn,samSortOut]
    
    sortResults = samtools_sort(samSortArgs)
    
    
    # +++ Index BAM +++ #
    samIndexArgs = [samSortOut+'.bam']
    
    indexResults = samtools_index(samIndexArgs)
    
