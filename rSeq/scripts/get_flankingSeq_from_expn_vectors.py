"""
1: Collect Tx from one or more species that are within at least some r value of similarity to
   a provided example Tx or a submitted hypothetical expression vector.
2: Use GTFs, BEDtools, and genome FASTAs to extract the upstream flanking sequences into a new FASTA
   for use in motif discovery.
"""

import sys,os
import argparse
import logging
import shutil
from tempfile import NamedTemporaryFile

import gtf_to_genes
import pybedtools
import numpy as np

from rSeq.utils.errors import *
from rSeq.utils.misc import Bag
from rSeq.utils.stats import basic_bootstrap_est
from rSeq.utils.files import fastaRec_length_indexer,tableFile2namedTuple,mv_file_obj
from rSeq.utils.expression import pearsonExpnFilter,mangle_expn_vectors
from rSeq.utils.externals import mkdirp


    

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
        chrom      = str(txDict[tx].gene.contig)
        chromStart = str(txDict[tx].beg)
        chromEnd   = str(txDict[tx].end)
        name       = str(txDict[tx].cdna_id)
        score      = str(999)
        if txDict[tx].gene.strand == True:
            strand = '+'
        else:
            strand = '-'
        
        bedLines.append([chrom,chromStart,chromEnd,name,score,strand])
        
    sortedLines = sorted(sorted(bedLines,key=lambda x: min(int(x[1]),int(x[2]))),key=lambda x: x[0])
    
    tmpFile.writelines(['\t'.join(l)+'\n' for l in sortedLines])
    tmpFile.close()
    
    return tmpFile
    
    
    
    
def get_fastas(txBed,genomeFasta,lenIndex,lenFlanks):
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
    flankBed = txBed.flank(g=lenIndex,l=lenFlanks,r=0,s=True)
    flankBed.sequence(fi=genomeFasta,name=True,s=True)
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
    parser.add_argument('--pval-filter-thresh', type=float, default=0.05,
                            help="""Set the upper threshold for the p-value of the Pearson r values to keep. (default: %(default)s)""")    
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
    parser.add_argument('--dump-megafasta', action='store_true',
                        help="""Save concatonated fasta file for debugging. (default: %(default)s)""")
    parser.add_argument('--dump-stats', action='store_true',
                            help="""Print a list of Tx/gene names and the r- p-values that passed the filter and exit without getting fastas. (default: %(default)s)""")    
    
    args = parser.parse_args()
    
    # tmp files will be stored here
    tmp_files = Bag()
    
    # 1: Use a correlation filter to pull out any Tx that is sufficiently similar to the model Tx
    vectDict = mangle_expn_vectors(expnPath=args.expn_path,txNameHeader=args.tx_name_header,condHeaders=args.cond_headers,manualHeaders=args.manual_headers)
    
    filterFunc = eval("lambda x: x %s %f" % (args.pearson_filter_type, args.pearson_filter_thresh))
    filterDict = pearsonExpnFilter(modelVector=vectDict[args.tx_name], targetVectors=vectDict, filterFunc=filterFunc)
    
    # remove vectors whose r's pVal is not significant (<=0.05)
    sigVectors = {}
    for key in filterDict:
        if key[1] <= args.pval_filter_thresh:
            sigVectors[key] = filterDict[key]
    
    # Impose a distance filter to further refine the gene set
    # incorperating magnitudes of the absolute levels of gene expression
    
    # set the boundries of acceptable deviation for the target gene mean expression
    # mangitude by bootstrapping.  The metric for comparison will be the average of
    # the differences of each point in remaining vectors against the target
    # vector.
    
    # 1) calc the metrics for each remaining gene's vector
    #    PS: numpy rocks.
    avgDists = {}
    for key in sigVectors:
        avgDist_i = np.mean(np.subtract(vectDict[args.tx_name],
                                           sigVectors[key]))
        avgDists[key] = avgDist_i
        
    # 2) bootstrap that bitch and give me a stdErr!
    medianEst,stdErrEst,lo95,hi95 = basic_bootstrap_est(avgDists.values())
    
    # 3) recover keys that fall within +/- 1 SE
    matchVectors = {}
    for key in avgDists:
        avgDist = avgDists[key]
        if (avgDist >= -stdErrEst) and (avgDist <= stdErrEst):
            matchVectors[key] = sigVectors[key]
    
        
    
    # Sort txList so that the highest r values are at the top
    # and save vectors and this info out to file
    txList = sorted(matchVectors.keys(),key=lambda x: x[0], reverse=True)
    sortedTxListFile = NamedTemporaryFile(mode='w+t',prefix='txExpnVectFilteredBy_r.',suffix=".tsv",delete=False)
    for row in txList:
        sortedTxListFile.write('%s\t%s\n' % ('\t'.join(map(str,row)),'\t'.join(matchVectors[row])))
    tmp_files['sortedTxListFile'] = sortedTxListFile
    sortedTxListFile.close()
    
    g2gObj = gtf_to_genes.get_indexed_genes_matching_gtf_file_name(index_file_name=args.gtf_index, logger=logger, regex_str=args.gtf)[-1]
    txDict = filter_GTF_4_Tx(txList=[x[2] for x in txList],g2gObj=g2gObj)
    tmp_files['txBedFile'] = convert_2_bed(txDict=txDict)
    
    # 2: Use GTFs, BEDtools, and genome FASTAs to extract the upstream flanking sequences into a new FASTA
    fastaRecLengths,fastaSeqs = fastaRec_length_indexer(fastaFiles=args.genome_fastas)
    tmpFastaRecLengthFile = NamedTemporaryFile(mode='w+b',prefix='tmpFastaRecLengthFile.',suffix=".txt")
    for seqRec in fastaRecLengths:
        tmpFastaRecLengthFile.write("%s\t%s\n" % (seqRec,fastaRecLengths[seqRec]))
    tmpFastaRecLengthFile.flush()

    # TODO: concatonate fasta files
    megaFastaFile = NamedTemporaryFile(mode='w+b',prefix='tmpMegaFastaFile.',suffix=".fas")
    for fasta in fastaSeqs:
        megaFastaFile.write('>%s\n%s\n' % (fasta,fastaSeqs[fasta]))
    megaFastaFile.flush()
        
    tmp_files['flankBed'] = get_fastas(txBed=tmp_files.txBedFile.name,genomeFasta=megaFastaFile.name,lenIndex=tmpFastaRecLengthFile.name,lenFlanks=args.flank_len)
    
    
    # CLEAN UP:
    # TODO: Close all tmp_files, and move to args.outDir
    mkdirp(args.out_dir)
    for f in tmp_files:
        try:
            tmp_files[f].delete = False
        except AttributeError:
            pass
        try:
            tmp_files[f].close()
        except AttributeError:
            pass
    # ['sortedTxListFile', 'flankBed', 'txBedFile', 'flankFasta']
    sortedTxListFile = "%s/sortedTxList.tsv" % (args.out_dir)
    flankBed         = "%s/flankBed.bed" % (args.out_dir)
    txBedFile        = "%s/txBed.bed" % (args.out_dir)
    flankFasta       = "%s/flankFasta.fas" % (args.out_dir)
    
    
    shutil.move(tmp_files.sortedTxListFile.name, sortedTxListFile)
    os.chmod(sortedTxListFile,0775)
    
    tmp_files.flankBed.saveas(flankBed)
    os.chmod(flankBed,0775)
    
    shutil.move(tmp_files.txBedFile.name, txBedFile)
    os.chmod(txBedFile,0775)
    
    shutil.move(tmp_files.flankBed.seqfn, flankFasta)
    os.chmod(flankFasta,0775)
    
    if args.dump_megafasta:
        megaFasta = "%s/megaFasta.fas" % (args.out_dir)
        megaFastaFile.delete = False
        megaFastaFile.close()
        shutil.move(megaFastaFile.name, megaFasta)
        os.chmod(megaFasta,0775)

    
if __name__ == "__main__":

    main()
