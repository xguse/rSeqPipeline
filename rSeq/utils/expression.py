from rSeq.utils.errors import *
import cogent
import pysam

def get_feature_reads(pyCgFeat,bams):
    """Given pyCogent feature, return reads (alignments optional) that map to that feature.

    pyCgFeat : pycogent feature object (exon,gene,etc)
    bams     : list of paths to BAM files (union of files will be used will be use)
    """
    if 1==1:
        raise NotImplementedError("ERROR: this function is under development.")
    
    # pysam-ize BAM files
    bams = [pysam.Samfile( x, "rb" ) for x in bams]
    
    # Record feature info
    
    #