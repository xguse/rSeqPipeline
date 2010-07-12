import sys
import optparse
from rSeq.utils.files import tableFile2namedTuple
from rSeq.utils.align import parseCigarString


strandReps = {'+':'+',
              '-':'-',
              '1':'+',
              '-1':'-',}

def groupFeatureAlignments(features,opts):
    """Returns Dict with key,val mappings == (EST_ID,ChrmID),[rowsFrom_features]"""
    featureDict = {}
    for row in features:
        k = (row.__getattribute__(opts.featName),row.__getattribute__(opts.chrm))
        if k in featureDict:
            featureDict[k].append(row)
        else:
            featureDict[k] = []
            featureDict[k].append(row)
    return featureDict



def printBEDline(listOfRowsByFeature,opts):
    """Takes a list of block info rows grouped by a single feature.
    Prints the BED line calculated from this information."""
    
    # +++++ func specific Defs +++++
    def getBlockSizes(feat):
        blkSzList = []
        for blk in feat:
            blkSzList.append(str(int(blk.__getattribute__(opts.blkChmEnd))-int(blk.__getattribute__(opts.blkChmStrt))+1))
        return ','.join(blkSzList)
    
    def getBlockStarts(feat,chrmStart):
        blkStrtList = []
        for blk in feat:
            blkStrtList.append(str(int(blk.__getattribute__(opts.blkChmStrt))-1-int(chrmStart)))
        return ','.join(blkStrtList)
    
    # +++ Do some data validation +++
    # --- Check Chrom info ---
    chrms = set([x.__getattribute__(opts.chrm) for x in listOfRowsByFeature])
    if len(chrms) != 1:
        print >> sys.stderr, 'ERROR: the %s feature has blocks on more than one chromosome/contig %s (%s column).  SKIPPING' % \
              (listOfRowsByFeature[0].__getattribute__(opts.featName),
               list(chrms),
               opts.chrm)
        return None
        #raise Exception, \
              #'ERROR: the %s feature seems to have blocks on more than one chromosome/contig\n%s (%s column)' % \
              #(listOfRowsByFeature[0].__getattribute__(opts.featName),
               #list(chrms),
               #opts.chrm)
    
    # --- Check strand info ---
    if opts.strand:
        strands = set([x.__getattribute__(opts.strand) for x in listOfRowsByFeature])
        
        # ensure all blocks on same strand
        if len(strands) != 1: 
            print >> sys.stderr, 'ERROR: the %s feature has more than one strand type: %s (%s column).  SKIPPING' % \
                  (listOfRowsByFeature[0].__getattribute__(opts.featName),
                   list(strands),
                   opts.strand)
            return None
            #raise Exception, \
                  #'ERROR: the %s feature has more than one strand type: %s (%s column)' % \
                  #(listOfRowsByFeature[0].__getattribute__(opts.featName),
                   #list(strands),
                   #opts.strand)
    
        # ensure the strand representation is recognized
        if not list(strands)[0] in strandReps: 
            raise Exception, \
                  "ERROR: the %s feature's strand representation is not recognized: %s (%s column)" % \
                  (listOfRowsByFeature[0].__getattribute__(opts.featName),
                   list(strands)[0],
                   opts.strand)
        
    # --- Check ThickEdge info ---
    # ensure all thkStarts agree
    if opts.thkStrt:
        thkStrts = set([x.__getattribute__(opts.thkStrt) for x in listOfRowsByFeature])
        if len(thkStrts) != 1: 
            raise Exception, \
                  'ERROR: the %s feature has more than one thickStart: %s (%s column)' % \
                  (listOfRowsByFeature[0].__getattribute__(opts.featName),
                   list(thkStrts),
                   opts.thkStrt)
    # ensure all thkEnds agree
    if opts.thkEnd:
        thkEnds = set([x.__getattribute__(opts.thkEnd) for x in listOfRowsByFeature])
        if len(thkEnds) != 1: 
            raise Exception, \
                  'ERROR: the %s feature has more than one thickEnd: %s (%s column)' % \
                  (listOfRowsByFeature[0].__getattribute__(opts.featName),
                   list(thkEnds),
                   opts.thkEnd)
    
    # +++ Sort based on blockStarts +++
    listOfRowsByFeature.sort(key=lambda x: int(x.__getattribute__(opts.blkChmStrt)))
    
    # +++ Gather the required pieces +++    
    chrm      = listOfRowsByFeature[0].__getattribute__(opts.chrm)
    chrmStart = str(int(listOfRowsByFeature[0].__getattribute__(opts.blkChmStrt))-1)
    chrmEnd   = listOfRowsByFeature[-1].__getattribute__(opts.blkChmEnd)
    name      = listOfRowsByFeature[0].__getattribute__(opts.featName)
    score     = '0'
    rgb       = opts.rgb
    blkCount  = str(len(listOfRowsByFeature))
    blkSizes  = getBlockSizes(listOfRowsByFeature)
    blkStarts = getBlockStarts(listOfRowsByFeature,chrmStart)
    # --- Initialize pieces that depend on user options ---
    if opts.strand:
        strand = strandReps[listOfRowsByFeature[0].__getattribute__(opts.strand)]
    else:
        strand = '+' 
    if opts.thkStrt:
        thkStart  = listOfRowsByFeature[0].__getattribute__(opts.thkStrt)
    else:
        thkStart  = chrmStart
    if opts.thkEnd:
        thkEnd    = listOfRowsByFeature[0].__getattribute__(opts.thkEnd)
    else:
        thkEnd    = chrmEnd
    
    print '%s' % ('\t'.join([chrm,    
                             chrmStart,
                             chrmEnd,
                             name,
                             score,
                             strand,
                             thkStart,
                             thkEnd,
                             rgb,
                             blkCount,
                             blkSizes,
                             blkStarts]))


if __name__ == "__main__":
    
    
    #+++++++++++ File Parseing Etc +++++++++++    
    usage = """python %prog inFile [options]"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-t',dest="track_name",type="str", default="untitled", 
                      help="""Space-less string for the display name for this track. (default=%default)""")
    parser.add_option('-d',dest="description",type="str", default="no description", 
                      help="""Quoted string: Exp: "Clone Paired Reads". (default=%default)""")
    parser.add_option('--sep',dest="sep",type="str", default="\t", 
                      help="""String used as delimiter. (default='\\t')""")
    parser.add_option('--featName',dest="featName",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the names of the features. Exp: transcript_name (default=%default)""")
    parser.add_option('--blkChmStrt',dest="blkChmStrt",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the chrom start positions of the blocks. Exp: seq_region_start (default=%default)""")
    parser.add_option('--blkChmEnd',dest="blkChmEnd",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the chrom end positions of the blocks. Exp: seq_region_end (default=%default)""")
    parser.add_option('--chrm',dest="chrm",type="str", default=None, 
                      help="""REQUIRED - Exact Title of Column holding the name of the chrom on which each block resides. Exp: name -or- Chromosome/plasmid (default=%default)""")
    parser.add_option('--strand',dest="strand",type="str", default=None, 
                      help="""Exact Title of Column holding the strand indicator. Exp: strand. If None: all will be assigned '+' strand. (default=%default)""")
    parser.add_option('--thkStrt',dest="thkStrt",type="str", default=None, 
                      help="""Exact Title of Column holding the position of the thickStart. Exp: start_of_coding (default=%default)""")
    parser.add_option('--thkEnd',dest="thkEnd",type="str", default=None, 
                      help="""Exact Title of Column holding the position of the thickEnd. Exp: end_of_coding (default=%default)""")
    parser.add_option('--rgb',dest="rgb",type="str", default="0,0,0", 
                      help="""Space-less string to assign the color for this track. Exp: 0,0,0=black; 255,0,0=red; 0,255,0=green; 0,0,255=blue (default=%default)""")
    parser.add_option('--cigars',dest="cigars",type="str", default=False, 
                      help="""Exact Title of Column holding the cigar strings. Exp: cigar_string (default=%default)""")
    cigTypes = ['ensembl','exonerate']
    parser.add_option('--cigar-type',dest="cigar_type",type="str", default=False, 
                      help="""Type of cigar string.  REQUIRED when using '--cigars'.  Options: %s (default=%default)""" % (cigTypes))

    
    (opts, args) = parser.parse_args()
    
    if len(args) != 1:
        parser.print_help()
        exit()
    
    
    features   = tableFile2namedTuple(args[0],sep=opts.sep)
    rowsByAlgn = groupFeatureAlignments(features,opts)
    
    print """track name=%s description="%s" useScore=0""" % (opts.track_name, opts.description)
    
    for alnmnt in rowsByAlgn:
        printBEDline(rowsByAlgn[alnmnt],opts)
    