import sys
import os

import scipy.stats as stats

from rSeq.utils.errors import *
import cogent
import pysam


def pearsonExpnFilter(modelVector, targetVectors, filterFunc=None):
    """
    modelVector:   ordered list of expression values per condition of model gene
    targetVectors: dict of ordered expression vectors for each gene with geneName as key
    filterFunc:    func that returns True if gene should be kept based on r value, False otherwise (usually a lambda func)
    
    Returns a dict of gene names and expression vectors whose pearson r values
    passed the filterFunc.
    """
    if filterFunc == None:
        filterFunc = lambda r: r >= 0.7
    rDict = {}
    
    for t in targetVectors:
        rStats = stats.pearsonr(modelVector,targetVectors[t])
        if filterFunc(rStats[0]):
            rDict[rStats + (t,)] = targetVectors[t]
    return rDict


#def get_feature_reads(pyCgFeat,bams):
    #"""Given pyCogent feature, return reads (alignments optional) that map to that feature.

    #pyCgFeat : pycogent feature object (exon,gene,etc)
    #bams     : list of paths to BAM files (union of files will be used will be use)
    #"""
    #if 1==1:
        #raise NotImplementedError("ERROR: this function is under development.")
    
    ## pysam-ize BAM files
    #bams = [pysam.Samfile( x, "rb" ) for x in bams]
    
    ## Record feature info
    
    ##