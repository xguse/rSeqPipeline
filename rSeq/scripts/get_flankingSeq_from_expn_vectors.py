"""
1: Collect Tx from one or more species that are within at least some r value of similarity to
   a provided example Tx or a submitted hypothetical expression vector.
2: Use GTFs, BEDtools, and genome FASTAs to extract the upstream flanking sequences into a new FASTA
   for use in motif discovery.
"""

import sys,os
import argparse
from tempfile import NamedTemporaryFile

import gtf_to_genes
import pybedtools

from rSeq.utils.errors import *
from rSeq.utils.misc import Bag
from rSeq.utils.files import fastaRec_length_indexer,tableFile2namedTuple,mv_file_obj
from rSeq.utils.expression import pearsonExpnFilter

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
    if manualHeaders != False:
        if type([]) != type(manualHeaders):
            raise SanityCheckError("manualHeaders must be type([]); you gave: %s." % (manualHeaders))
    
    # lets go    
    expnTable = tableFile2namedTuple(tablePath=expnPath,sep='\t',headers=manualHeaders)
    
    vectDict = {}
    for row in expnTable:
        vectDict[row.__getattribute__(txNameHeader)] = [row.__getattribute__(x) for x in condHeaders]
    
    return vectDict
    

def filter_GTF_4_Tx(txList,g2gObj):
    """
    GIVEN:
    1) txList: list of requested Tx names
    2) g2gObj: a gtf_to_genes obj
    
    RETURN:
    1) a dict of Tx g2g dicts matching supplied Tx list
    """
    keptTx = {}
    
    txIdxDict = gtf_to_genes.index_transcripts(g2gObj, by_prot_id=False)
    
    for t in txList:
        if t in txIndex:
            keptTx[t] = txIdxDict[t]
        else:
            # TODO: !*!*! add logging code here !*!*!
            pass
    return keptTx
    
def convert_2_bed(txDict):
    """
    GIVEN:
    1) txDict: a dict of Tx g2g dicts
    
    DO:
    1) construct a bed entry for every Tx in list with start,stop,chrom,strand
    2) write this bed to a tmp file sorted by chrom:minCoord
    
    RETURN:
    1) path to out file
    """
    # NOTE: because gtf_to_genes uses the same coord scheme there is no need for coord
    # mangling in a simple BED representation (we dont care about exon blocks here)
    
    tmpFile = NamedTemporaryFile(mode='w+t',suffix=".bed",delete=False)
    
    bedLines = []
    
    for tx in txDict:
        chrom      = str(tx.gene.contig)
        chromStart = str(tx.beg)
        chromEnd   = str(tx.end)
        name       = str(tx.cdna_id)
        score      = str(999)
        strand     = str(tx.gene.strand)
        
        bedLines.append([chrom,chromStart,chromEnd,name,score,strand])
        
    sortedLines = sorted(sorted(bedLines,key=lambda x: min(int(x[1]),int(x[2]))),key=lambda x: x[0])
    
    tmpFile.writelines(['\t'.join(l)+'\n' for l in sortedLines])
    tmpFile.close()
    
    return tmpFile
    
    
    
    
def write_fastas(txBed,genomeFasta,lenIndex,lenFlanks,tmpFileDict):
    """
    GIVEN:
    1) txBed: path to TxBED tmp file
    2) genomeFasta: species genome fasta path
    3) lenIndex: chrom length index path
    4) lenFlanks: length of flank regions
    5) tmpFileDict: tmp_files dict so that the tmp beds can be recovered
    
    DO:
    1) create new bed from TxBED of flanking regions
    2) get and write seqs defined in flankBED to fasta tmp file
    
    RETURN:
    1) BedTool obj containing flankBed coords file (flankBed.fn) and fastaSeqFilePath (flankBed.seqfn)
    """
    
    # create flankBed
    txBed = pybedtools.BedTool(txBed)
    flankBed = txBed.flank(genome=lenIndex,l=lenFlanks,r=0,s=True).sequence(fi=genomeFasta,name=True,s=True)
    tmpFileDict['flankBed'] = flankBed
    return flankBed

def main(args):
    """
    1: Collect Tx from one or more species that are within at least some r value of similarity to
       a provided example Tx or a submitted hypothetical expression vector.
    2: Use GTFs, BEDtools, and genome FASTAs to extract the upstream flanking sequences into a new FASTA
       for use in motif discovery.
    """
    
    desc = """(1) Collect Tx from one or more species that are within 
at least some r value of similarity to a provided example Tx or a 
submitted hypothetical expression vector. (2) Use GTFs, BEDtools, and 
genome FASTAs to extract the upstream flanking sequences into a new 
FASTA for use in motif discovery."""
    
    parser = argparse.ArgumentParser(description=desc)
    
    parser.add_argument('--tx-list', dest='txList', type=str, required=True, nargs='+',
                        help="""help text goes here %(default)s""")
    parser.add_argument('--gtf', dest='gtf', type=str, required=True,
                        help="""help text goes here %(default)s""")
    parser.add_argument('--genome-fastas', dest='genomeFastas', type=str, required=True, nargs='+',
                        help="""(allow cat of files and/or dir contents)help text goes here %(default)s""")
    parser.add_argument('--flank-len', dest='lenFlanks', type=int, required=True, nargs='+',
                        help="""help text goes here %(default)s""")
    
    # tmp files will be stored here
    tmp_files = Bag()
    
    # 1: Use a correlation filter to pull out any Tx that is sufficiently similar to the model Tx
    vectDict = mangle_expn_vectors(expnPath,txNameHeader,condHeaders,manualHeaders=args.toSet)
    
    filterDict = pearsonExpnFilter(modelVector=vectDict[args.toSet], targetVectors=vectDict, filterFunc=None)
    
    # remove vectors whose r's pVal is not significant (<=0.05)
    matchVectors = {}
    for key in filterDict:
        if key[1] <= args.rpval:
            matchVectors[key] = filterDict[key]
    
    # Sort txList so that the highest r values are at the top
    # and save vectors and this info out to file
    txList = sorted(matchVectors.keys(),key=lambda x: x[0], reverse=True)
    sortedTxListFile = NamedTemporaryFile(mode='w+t',prefix='txExpnVectFilteredBy_r.',suffix=".tsv",delete=False)
    for row in txList:
        sortedTxListFile.write('%s\t$s\n' % ('\t'.join(row), '\t'.join(txDict[row])))
    tmp_files['sortedTxListFile'] = sortedTxListFile
    sortedTxListFile.close()
    
    txDict = filter_GTF_4_Tx(txList=[x[2] for x in txList],g2gObj=args.toSet)
    tmp_files['txBedFile'] = convert_2_bed(txDict=txDict)
    
    # 2: Use GTFs, BEDtools, and genome FASTAs to extract the upstream flanking sequences into a new FASTA
    fastaRecLengths = fastaRec_length_indexer(fastaFiles=args.toSet)
    tmpFastaRecLengthFile = NamedTemporaryFile(mode='w+b',prefix='tmpFastaRecLengthFile.',suffix=".fas")
    for seqRec in fastaRecLengths:
        tmpFastaRecLengthFile.write("%s\t%s\n" % (seqRec,fastaRecLengths[seqRec]))
    tmpFastaRecLengthFile.flush()
        
    flankBed = write_fastas(txBed=tmp_files.txBedFile.name,genomeFasta=tmp_files.megaFastaFile.name,lenIndex=tmpFastaRecLengthFile.name,lenFlanks=args.toSet)
    # 
    
    # CLEAN UP:
    # TODO: Close all tmp_files, and move to args.outDir
    for f in tmp_files:
        tmp_files[f] = mv_file_obj(tmp_files[f],args.outDir,0755)
        tmp_files[f].close()
        
    
if __name__ == "__main__":

    main(args=sys.argv)
