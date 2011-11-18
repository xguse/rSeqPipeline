import sys
import os

import cogent
import pysam
import scipy.stats as stats

from rSeq.utils.errors import *
from rSeq.utils.files import tableFile2namedTuple




def mangle_expn_vectors(expnPath,txNameHeader,condHeaders,manualHeaders=False):
    """
    GIVEN:
    1) expnPath: expn file path.
    2) txNameHeader: header name where txName lives.
    3) condHeaders: list of header names of conditions to be
       included (order should be meaningful).
    4) manualHeaders: manually set headers for the expn table file
       (provide header name for EVERY column in file ONLY if no headers are already present).
    
    DO:
    1) Extract the txName and condition expn levels from expnPath
       and store in vectDict with txName as key and expn vector(list) as value.
       
    RETURN:
    1) vectDict: see "DO:".
    """
    # Sanity checks
    if type('') != type(txNameHeader):
        raise SanityCheckError("txNameHeader must be type(''); you gave: %s." % (txNameHeader))
    if type([]) != type(condHeaders):
        raise SanityCheckError("condHeaders must be type([]); you gave: %s." % (condHeaders))
    if manualHeaders:
        if type([]) != type(manualHeaders):
            raise SanityCheckError("manualHeaders must be type([]); you gave: %s." % (manualHeaders))
    
    # lets go    
    expnTable = tableFile2namedTuple(tablePath=expnPath,sep='\t',headers=manualHeaders)
    
    vectDict = {}
    for row in expnTable:
        vectDict[row.__getattribute__(txNameHeader)] = [float(row.__getattribute__(x)) for x in condHeaders]
    
    return vectDict

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
        rStats = stats.pearsonr([float(x) for x in modelVector],[float(x) for x in targetVectors[t]])
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