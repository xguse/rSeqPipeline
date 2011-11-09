"""
1: Collect Tx from one or more species that are within at least some r value of similarity to
   a provided example Tx or a submitted hypothetical expression vector.
2: Use GTFs, BEDtools, and genome FASTAs to extract the upstream flanking sequences into a new FASTA
   for use in motif discovery.
"""

import sys,os
import argparse
import logging
from tempfile import NamedTemporaryFile

import gtf_to_genes
import pybedtools

from rSeq.utils.errors import *
from rSeq.utils.misc import Bag
from rSeq.utils.files import fastaRec_length_indexer,tableFile2namedTuple,mv_file_obj
from rSeq.utils.expression import pearsonExpnFilter
from rSeq.utils.externals import mkdirp

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
        if t in txIdxDict:
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

def main():
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
    FileType = argparse.FileType
    
    logger = logging.getLogger(sys.argv[0].split('/')[-1])
    
    parser.add_argument('--expn-path', type=str, required=True,
                        help="""Path to expression table file. \n(default: %(default)s)""")
    parser.add_argument('--tx-name', type=str, required=True,
                        help="""Name of the Tx you want to use as a model. (default: %(default)s)""")
    parser.add_argument('--pearson-filter-type', type=str, default='>=', choices=['>=','<='],
                        help="""Use >= to find similar expn profiles or <= to find opposite profiles. (default: %(default)s)""")
    parser.add_argument('--pearson-filter-thresh', type=float, default=0.7,
                        help="""Set the threshold of the Pearson r value for the filter. (default: %(default)s)""")
    parser.add_argument('--tx-name-header', type=str, required=True,
                        help="""The text of the header in the expn table where tx names are stored. (default: %(default)s)""")
    parser.add_argument('--cond-headers', type=str, required=True, nargs='+',
                        help="""A list of the text of the headers in the expn table where the values for each condition are stored (--cond-headers cond1 cond2 ...). (default: %(default)s)""")
    parser.add_argument('--manual-headers', type=str, required=False, nargs='?',
                        help="""If the expn table does not have headers, provide a list of ordered names for them here. (default: %(default)s)""")
    parser.add_argument('--gtf', type=str, required=True,
                        help="""The path to the gtf file that you want to use for your annotation. (default: %(default)s)""")
    parser.add_argument('--gtf-index', type=str, required=True,
                        help="""The path to the gtf index file generated from "gtf_to_genes". (default: %(default)s)""")
    parser.add_argument('--genome-fastas', type=str, required=True, nargs='+',
                        help="""A list of paths to genomic fasta files or directories where they are stored. (default: %(default)s)""")
    parser.add_argument('--flank-len', type=int, default=2000,
                        help="""The length in bp that should be harvested from the 5' end of the tx. (default: %(default)s)""")
    parser.add_argument('--out-dir', type=str, default='.',
                        help="""A path to a directory where you would like the output files to be stored. (default: %(default)s)""")
    
    args = parser.parse_args()
    
    # tmp files will be stored here
    tmp_files = Bag()
    
    # 1: Use a correlation filter to pull out any Tx that is sufficiently similar to the model Tx
    vectDict = mangle_expn_vectors(expnPath=args.expn_path,txNameHeader=args.tx_name_header,condHeaders=args.cond_headers,manualHeaders=args.manual_headers)
    
    filterFunc = eval("lambda x: x %s %f" % (args.pearson_filter_type, args.pearson_filter_thresh))
    filterDict = pearsonExpnFilter(modelVector=vectDict[args.tx_name], targetVectors=vectDict, filterFunc=filterFunc)
    
    # remove vectors whose r's pVal is not significant (<=0.05)
    matchVectors = {}
    for key in filterDict:
        if key[1] <= 0.05:
            matchVectors[key] = filterDict[key]
    
    # Sort txList so that the highest r values are at the top
    # and save vectors and this info out to file
    txList = sorted(matchVectors.keys(),key=lambda x: x[0], reverse=True)
    sortedTxListFile = NamedTemporaryFile(mode='w+t',prefix='txExpnVectFilteredBy_r.',suffix=".tsv",delete=False)
    for row in txList:
        sortedTxListFile.write('%s\n' % ('\t'.join(map(str,row))))
    tmp_files['sortedTxListFile'] = sortedTxListFile
    sortedTxListFile.close()
    
    g2gObj = gtf_to_genes.get_indexed_genes_matching_gtf_file_name(index_file_name=args.gtf_index, logger=logger, regex_str=args.gtf)[-1]
    txDict = filter_GTF_4_Tx(txList=[x[2] for x in txList],g2gObj=g2gObj)
    tmp_files['txBedFile'] = convert_2_bed(txDict=txDict)
    
    # 2: Use GTFs, BEDtools, and genome FASTAs to extract the upstream flanking sequences into a new FASTA
    fastaRecLengths = fastaRec_length_indexer(fastaFiles=args.genome_fastas)
    tmpFastaRecLengthFile = NamedTemporaryFile(mode='w+b',prefix='tmpFastaRecLengthFile.',suffix=".fas")
    for seqRec in fastaRecLengths:
        tmpFastaRecLengthFile.write("%s\t%s\n" % (seqRec,fastaRecLengths[seqRec]))
    tmpFastaRecLengthFile.flush()
        
    flankBed = write_fastas(txBed=tmp_files.txBedFile.name,genomeFasta=tmp_files.megaFastaFile.name,lenIndex=tmpFastaRecLengthFile.name,lenFlanks=args.flank_len)
    # 
    
    # CLEAN UP:
    # TODO: Close all tmp_files, and move to args.outDir
    mkdirp(args.out_dir)
    for f in tmp_files:
        tmp_files[f] = mv_file_obj(tmp_files[f],args.outDir,0755)
        tmp_files[f].close()
        
    
if __name__ == "__main__":

    main()
